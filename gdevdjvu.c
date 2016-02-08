/* -*- Mode: c; c-basic-offset: 4; -*- 
   
   -------------------------------------------------------------------------- 
   DjVu Device for Ghostscript 
   -- Copyright (C) 2000 AT&T Corp.
   -- Copyright (C) 2002-2007 Leon Bottou.
   ------------------------------------------------------------------------ 

   This file is derived from the gsdjvu files released in June 2005 
   by AT&T Corp. under the terms of the Common Public Licence. 
   See the files COPYING and COPYING.CPL for more information.

   This file includes further contributions authored by Leon Bottou,
   including (i) the gathering of all gsdjvu code into a single 
   self-contained file, (ii) various bug fixes and compatibility fixes, 
   and (iii) support for text extraction.  These contributions 
   are governed by the Common Public License.

   Alternatively, Leon Bottou also releases these contributions under 
   the terms of the GNU Public License.  Should AT&T re-release 
   the gsdjvu code under the GNU Public Licence, the present file 
   would become distributable under the terms of the GNU Public License.
   As of June 2005, this is not the case.

   ------------------------------------------------------------------------ 

   This driver implements two devices:

   * Device "djvumask" outputs the foreground mask as PBM files.

   * Device "djvusep" outputs separation files containing one or several
     pages. Each page contains the foreground encoded as a CRLE image
     (cf. comments on Color RLE image below), possibly followed by a
     subsampled PPM image representing the background, possibly followed by an
     arbitrary number of comment lines starting with character '#'.
     These files should be piped through the back-end encoder "csepdjvu".
   
   The following ghostscript command line options are supported:

   * "-sOutputFile=<string>"  Selects the name of the output file.
                              The usual Ghostscript conventions apply.

   * "-dThreshold=<0..100>"   Selects the separation threshold.  Larger
                              values put more things into the foreground
                              layer.  Default value is 80

   * "-dBgSubsample=<1-12>"   Selects the background subsampling factor.
                              Default value is 3.

   * "-dFgColors=<1-4000>"    Selects the maximal number of colors in the
                              foreground layer. Default value is 256.

   * "-dFgImgColors=<0-4000>" Images with less than the specified number of
                              colors might be moved into the foreground.  
                              Default value is 256.

   * "-dMaxBitmap=<num>"      Specifies the maximal amount of memory used
                              by the banding code when rendering the background
                              layer.  Default is 10000000.

   * "-dAutoHires"            This option enables an algorithm that lets
                              device "djvusep" ignore -dBgSubsample and 
			      generate a full resolution background layer
			      when the foreground is empty.

   * "-dExtractText"          This option causes the generation of sep file 
                              comments describing the textual information.

   * "-dReopenPerPage"        Specifies that the output file must be reopened
                              for each page.  Default is to reopen if the 
			      filename contains a "%d" specification for 
			      the page number.

   ------------------------------------------------------------------------ */   
   
#include "stdpre.h"
#include "arch.h"
#include "ctype_.h"
#include "stdio_.h"
#include "memory_.h"
#include "string_.h"
#include "math_.h"
#include <stdlib.h>   /* for qsort */
#include "gserrors.h"
#include "gsmemory.h"
#include "gsdevice.h"
#include "gstypes.h"
#include "gscdefs.h"
#include "gsrect.h"
#include "gscoord.h"
#include "gsstruct.h"
#include "gsmalloc.h"
#include "gsparam.h"
#include "gspath.h"
#include "gx.h"
#include "gxstdio.h"
#include "gxdevice.h"
#include "gxdevmem.h"
#include "gxbitmap.h"
#include "gxfont.h"
#include "gxgetbit.h"
#include "gximage.h"
#include "gxlum.h"
#include "gxalloc.h"
#include "gxcindex.h"
#include "gxstdio.h"
#include "gdevprn.h"

#ifndef GS_VERSION
# define GS_VERSION 0
#endif

/* Debugging characters for gs option -z */
#define DEBUG_CHAR_DLIST    0
#define DEBUG_CHAR_COLOR    0
#define DEBUG_CHAR_P2MEM    0
#define DEBUG_CHAR_CHROME   0
#define DEBUG_CHAR_LOCOLOR  0
#define DEBUG_CHAR_BG       0

/* For ghostscript-8.61 and above */
#ifndef private
# define private static
#endif


/* ====================================== 
              U T I L I T I E S
   ======================================  */

/* Return memory error code */
#define return_VMerror \
    return_error(gs_error_VMerror)

/* Check for internal errors. */
#define ASSERT(cond) \
  while (!(cond)) lbassertfail(__FILE__,__LINE__);
#ifdef __GNUC__
private void lbassertfail(const char *, int) __attribute__ ((noreturn));
extern void abort(void) __attribute__ ((noreturn));
#endif

/* Forward declarations */
typedef struct gx_device_djvu_s gx_device_djvu;

/* STDOUT - for regular printout */
#ifdef gs_stdout
# define STDOUT gs_stdout
#else
# define STDOUT (((gx_device*)cdev)->memory->gs_lib_ctx->fstdout)
#endif

/* Called by macro ASSERT defined in the H file */
private void 
lbassertfail(const char *file, int line)
{
    FILE *serr = fdopen(2,"a");
    fprintf(serr,"Internal error at %s:%d\n", file, line);
    fclose(serr);
    abort();
}


/* ======================================
              P 2 M E M
   ====================================== */

/* Runmaps were initially putting a lot of stress on the Ghostscript standard
   memory allocator, causing fragmentation and creating lots of long-lived
   small objects.  The fragmentation issue was the major concern in gs6_23,
   sometimes accounting for more than 95% of the computation time.  Allocator
   improvements by L. Peter Deutsch basically eliminated that problem.  Yet
   the Ghostscript garbage collector was still spending a lot of time scanning
   these long lists of small objects, never finding memory to reclaim.  The
   following memory allocator was included in order to ``stay out of the
   way''.  Memory is obtained from the low level GS allocator in chunks.  Some
   chunks merely contain a single large object (large blocks).  The other
   chunks contain a number of small objects (small blocks).  Small block sizes
   are always constrained to be a power-of-two, minus two marker bytes located
   before and after the block. The marker byte simply contains the log2 of the
   block size (3 for 8 bytes, 4 for 16 bytes, etc.).  A marker byte of 255
   indicates a large block allocated alone in its chunk. */


/* Required alignment for p2mem allocated objects */
#define p2mem_log2_align  log2_obj_align_mod
#define p2mem_align       obj_align_mod

/* Parent allocator */
#define p2mem_parent_alloc(size)           malloc(size)
#define p2mem_parent_resize(data, newsize) realloc(data, newsize)
#define p2mem_parent_free(data)            free(data)

/* Internal: minimal and maximal size of small objects */
#ifndef arch_log2_sizeof_ptr
#if arch_sizeof_ptr == 4
#define arch_log2_sizeof_ptr 2
#else
#if arch_sizeof_ptr == 8
#define arch_log2_sizeof_ptr 3
#endif
#endif
#endif
#if p2mem_log2_align > arch_log2_sizeof_ptr
#define p2mem_log2_min    p2mem_log2_align
#define p2mem_log2_max    (11)
#else
#define p2mem_log2_min    (arch_log2_sizeof_ptr+1)
#define p2mem_log2_max    (11)
#endif
#define p2mem_min         (1<<p2mem_log2_min)
#define p2mem_max         (1<<p2mem_log2_max)
#define p2mem_chunksize   (31)

/* Memory allocator */
typedef struct p2mem_s {
    byte *cc;  /* slist of chunks containing small blocks. */
    byte *cb;  /* dlist of chunks containing a single large block. */
    uint cc_used; 
    uint cc_total;
    byte *freelist[p2mem_log2_max+1];
} p2mem;

/* Create an allocator. */
private p2mem *
p2mem_create(void) 
{
    p2mem *mem = p2mem_parent_alloc(sizeof(p2mem));
    if (mem) memset(mem, 0, sizeof(p2mem));
    return mem;
}

/* Free everything in an allocator. */
private void
p2mem_freeall(p2mem *mem) 
{
    if (mem) {
        while (mem->cc) {
            byte *cc = mem->cc;
            mem->cc = *(byte**)cc;
            p2mem_parent_free(cc);
        }
        while (mem->cb) {
            byte *cb = mem->cb;
            mem->cb = *(byte**)cb;
            p2mem_parent_free(cb);
        }
    }
    memset(mem, 0, sizeof(p2mem));
}

/* Destroy an allocator. */
private void
p2mem_destroy(p2mem *mem) 
{
    if (mem) {
        p2mem_freeall(mem);
        p2mem_parent_free(mem);
    }
}

/* Internal: get a small block. */
private byte *
p2mem_getblock(p2mem *mem, uint log2)
{
    if (mem->freelist[log2]) {
        /* Get it from freelist */
        byte *data = mem->freelist[log2];
        mem->freelist[log2] = *(byte**)data;
        mem->cc_used += (1<<log2);
        return data;
    } else if (log2 < p2mem_log2_max) {
        /* Get it by splitting larger block */
        int size = (1<<log2);
        byte *data = p2mem_getblock(mem, log2+1);
        if (!data) return 0;
        data[-1] = data[size-2] = log2;
        data[size-1] = data[size+size-2] = log2;
        *(byte**)&data[size] = mem->freelist[log2];
        mem->freelist[log2] = &data[size];
        mem->cc_used -= size;
        return data;
    } else {
        /* Allocate new chunk */
        int overhead = sizeof(byte*) | ((1<<p2mem_log2_align)-1);
        int chunksize = overhead + p2mem_chunksize * p2mem_max;
        byte *cc = p2mem_parent_alloc(chunksize);
        byte *data = cc + chunksize + 1;
        int nblocks = p2mem_chunksize;
        if (!cc) return 0;
        *(byte**)cc = mem->cc;
        mem->cc = cc;
        mem->cc_total += chunksize;
        while (nblocks-- > 0) {
            data -= p2mem_max;
            data[-1] = data[p2mem_max-2] = p2mem_log2_max;
            *(byte**)data = mem->freelist[p2mem_log2_max];
            mem->freelist[p2mem_log2_max] = data;
        }
        /* Return block */
        data = mem->freelist[p2mem_log2_max];
        mem->freelist[p2mem_log2_max] = *(byte**)data;
        mem->cc_used += p2mem_max;
        return data;
    }
}

/* Allocate data. */
private void *
p2mem_alloc(p2mem *mem, uint size) 
{
    if (size > 0) {
        size = size + 2;
        if (size <= p2mem_max) {
            /* Compute suitable power of two */
            int log2 = p2mem_log2_min;
            int blocksize = (1<<log2);
            while (blocksize < size) {
                log2++;
                blocksize <<= 1;
            }
            /* Allocate */
            return p2mem_getblock(mem, log2);
        } else {
            /* Large block */
            int overhead = 2*sizeof(byte*) + ((1<<p2mem_log2_align)-1);
            byte *cb = p2mem_parent_alloc(size + overhead);
            if (!cb) return 0;
            if (mem->cb) { ((byte**)(mem->cb))[1] = cb; }
            ((byte**)cb)[0] = mem->cb;
            ((byte***)cb)[1] = &(mem->cb);
            mem->cb = cb;
            cb[overhead] = 0xff;
            cb[overhead+size-1] = 0xff;
            return cb + overhead + 1;
        }
    }
    return 0;
}

/* Free data block. */
private void
p2mem_free(p2mem *mem, void *vdata)
{
    if (vdata) {
        byte *data = vdata;
        int log2 = data[-1];
        if (log2 <= p2mem_log2_max) {
            /* Small blocks */
            int blocksize = (1<<log2);
            ASSERT(log2>=p2mem_log2_min);
            ASSERT(data[blocksize-2] == log2);
            *(byte**)data = mem->freelist[log2];
            mem->freelist[log2] = data;
            mem->cc_used -= blocksize;
        } else {
            /* Large blocks */
            int overhead = 2*sizeof(byte*) + ((1<<p2mem_log2_align)-1);
            byte *cb = data - overhead - 1;
            byte *cb2 = ((byte**)cb)[0];
            ASSERT(log2 == 0xff);
            if (cb2) { ((byte***)cb2)[1] = ((byte***)cb)[1]; }
            ((byte***)cb)[1][0] = ((byte**)cb)[0];
            p2mem_parent_free(cb);
        }
    }
}

/* Resize data block. */
private void *
p2mem_resize(p2mem *mem, void *vdata, uint newsize)
{
    if (newsize == 0) {
        p2mem_free(mem, vdata);
    } else if (vdata == 0) {
        return p2mem_alloc(mem, newsize);
    } else {
        byte *data = vdata;
        int log2 = data[-1];
        int blocksize = (1<<log2);
        newsize += 2;
        if (log2 <= p2mem_log2_max) {
            /* Small block reallocation */
            if (newsize <= blocksize) {
                if (newsize < p2mem_min)
                    newsize = p2mem_min;
                /* No change if blocksize is adequate. */
                if (newsize+newsize > blocksize)
                    return data;
                /* Otherwise fallback to copy.  Splitting the block until
                   reaching the right blocksize causes very inefficient
                   fragmentation. */
            }
        } else if (newsize >= p2mem_max) {
            /* Large block reallocation */
            int overhead = 2*sizeof(byte*) + ((1<<p2mem_log2_align)-1);
            byte *cb = data - overhead - 1;
            byte *cb2 = ((byte**)cb)[0];
            ASSERT(log2 == 0xff);
            if (cb2) { ((byte***)cb2)[1] = ((byte***)cb)[1]; }
            ((byte***)cb)[1][0] = ((byte**)cb)[0];
            cb2 = p2mem_parent_resize(cb, newsize+overhead);
            if (!cb2) return 0;
            if (mem->cb) { ((byte**)(mem->cb))[1] = cb2; }
            ((byte**)cb2)[0] = mem->cb;
            ((byte***)cb2)[1] = &(mem->cb);
            mem->cb = cb2;
            cb2[overhead] = 0xff;
            cb2[overhead+newsize-1] = 0xff;
            return cb2 + overhead + 1;
        }
        /* Fallback: copy data */
        data = p2mem_alloc(mem, newsize);
        if (!data) return 0;
        if (log2 > p2mem_log2_max) blocksize=newsize;
        memcpy(data, vdata, ((blocksize<newsize)?blocksize:newsize)-2);
        p2mem_free(mem, vdata);
        return data;
    }
    return 0;
}

#ifdef DEBUG
private void
p2mem_diag(p2mem *mem, gx_device_djvu *cdev)
{
    int i;
    byte *ptr, **pptr;
    int cbnum = 0;
    int cctotalfree = 0;
    fprintf(STDOUT,"cc: %d-%d=%d  (", 
            mem->cc_total, mem->cc_used, 
            mem->cc_total - mem->cc_used );
    for (i=p2mem_log2_min; i<p2mem_log2_max+1; i++) {
        int ccfree = 0;
        ptr = mem->freelist[i];
        while (ptr) { ccfree += 1; ptr = *(byte**)ptr; }
        cctotalfree += (1<<i) * ccfree;
        fprintf(STDOUT," %d", ccfree);
    }
    fprintf(STDOUT," ) = %d ", cctotalfree);
    pptr = &mem->cb;
    while ((ptr = *pptr)) {
        ASSERT( ((byte***)ptr)[1] == pptr );
        pptr = (byte**)ptr;
        cbnum += 1;
    }
    fprintf(STDOUT," cb: %d\n", cbnum);
}
#endif








/* ======================================
                R U N M A P
   ====================================== */

/* Memory efficient algorithms for performing boolean operations between run
   length encoded bitonal images.  Each row is represented by the lengths of
   alternating white and black pixels.  Each length can be encoded on one to
   three bytes.  Data for row <y> is located in array <data> between offsets
   <rows[y-ymin]> and <rows[y-ymin+1]>.  Empty runmaps are represented by null
   pointers */

typedef struct runmap_s {
    uint ymin, ymax, xmin, xmax;
    uint rlen, dlen;
    p2mem *mem;
    uint  *rows;
    byte  *data;
} runmap;

/* Deallocate runmap */
private void
runmap_free(runmap *a)
{
    if (a) { 
        p2mem *mem = a->mem;
        p2mem_free(mem, a->data);
        p2mem_free(mem, a->rows);
        p2mem_free(mem, a);
    }
}

/* Encode run length x at location data */
#define rl_encode(data, x) { \
  if ((uint)(x)<0xc0) { *data++ = (x); } \
  else if ((uint)(x)<0x2000) { \
         *data++ = (((x)>>8)|0xc0); *data++ = ((x)&0xff); } \
  else if ((uint)(x)<0x200000) { *data++ = (((x)>>16)|0xe0); \
         *data++ = ((x)&0xff00)>>8; *data++ = ((x)&0xff); } \
  else ASSERT(0); }

/* Decode run length at location data and increment x */
#define rl_decode(data, x) { \
  if (data[0] < 0xc0) { x += data[0]; data+=1; } \
  else if (data[0] < 0xe0) { x += (data[1]|(data[0]<<8))-0xc000; data+=2; } \
  else { x += (data[2]|(data[1]<<8)|(data[0]<<16))-0xe00000; data+=3; } }

/* Check that data is available and decode run length.*/
#define rl_check_decode(data, dataend, x, c, finished) { \
  if (data < dataend) { c = !c; rl_decode(data,x); } \
  else { x = ARCH_MAX_UINT; finished; } }

/* Accumulating small objects into a large runmap can be very inefficient
   because each call to runmap_or must copy the unchanged part of the large
   runmap.  The problem is easily circumvented by using a smaller intermediate
   runmap.  The small objects are accumulated into the intermediate runmap.
   The intermediate runmap is accumulated into the large runmap whenever its
   size exceeds a predefined threshold.  The following macro checks whether it
   is time to transfer the intermediate runmap into the large runmap. */
#define runmap_becoming_big(r)  ((r) && (r)->dlen > 2048)

/* Check whether runmaps a and b collide. */
private bool
runmap_hittest(runmap *a, runmap *b)
{
    int y;
    int ymin, ymax;
    /* determine ymin ymax */
    if (!a || !b) 
        return false;
    ymin = max(a->ymin, b->ymin);
    ymax = min(a->ymax, b->ymax);
    if (ymin > ymax || a->xmin > b->xmax || a->xmax < b->xmin) 
        return false;
    /* loop over lines */
    for (y=ymin; y<=ymax; y++) {
        /* access rows */
        byte *adata = a->data + a->rows[y - a->ymin];
        byte *adataend = a->data + a->rows[y - a->ymin + 1];
        byte *bdata = b->data + b->rows[y - b->ymin];
        byte *bdataend = b->data + b->rows[y - b->ymin + 1];
        /* process row */
        if (bdata < bdataend && adata < adataend) {
            /* something in both rows a and b */
            uint ax = 0;
            uint bx = 0;
            uint cx = 0;
            bool ac = true;
            bool bc = true;
            for(;;) {
                if (ax <= cx) 
		    rl_check_decode(adata, adataend, ax, ac, break);
                if (bx <= cx) 
		    rl_check_decode(bdata, bdataend, bx, bc, break);
		cx = min(ax,bx);
		if (ac && bc)
                    return true;
            }                                   
        }
    }
    /* finish */
    return false;
}

/* Internal: allocate runmap structure. */
private int 
runmap_alloc(p2mem *mem, int dlen, int rlen, runmap **out)
{
    uint *r;
    byte *d;
    runmap *a;
    dlen = max(1,dlen);
    rlen = max(1,rlen);
    a = p2mem_alloc(mem, sizeof(runmap));
    r = p2mem_alloc(mem, sizeof(uint)*rlen);
    d = p2mem_alloc(mem, dlen);
    if (! (a && d && r))
        return_VMerror;
    memset(a, 0, sizeof(runmap));
    a->mem = mem;
    memset(r, 0, sizeof(uint)*rlen);
    a->ymin = ARCH_MAX_UINT;
    a->xmin = ARCH_MAX_UINT;
    a->rlen = rlen;
    a->rows = r;
    a->dlen = dlen;
    a->data = d;
    *out = a;
    return 0;
}

/* Internal: cleanup after boolean operations. */
private int
runmap_finish(p2mem *mem, runmap *c, 
	      runmap *a, runmap *b, int del, runmap **out)
{
    int code;
    /* Check whether a is empty */
    if (c && c->ymin > c->ymax) {
        runmap_free(c);
        c = 0;
    }
    /* Definitely not empty */
    if (c) {
        runmap *from = 0;
        /* Check if we must make a copy */
        if (c == a) { if (del & 1) { a = 0; } else { from = a; } }
        if (c == b) { if (del & 2) { b = 0; } else { from = b; } }
        if (from) {
            /* Make the copy */
            code = runmap_alloc(mem, from->dlen, from->rlen, out);
            if (code < 0) return code;
            c = *out;
            c->ymin = from->ymin;
            c->ymax = from->ymax;
            c->xmin = from->xmin;
            c->xmax = from->xmax;
            memcpy(c->data, from->data, c->dlen);
            memcpy(c->rows, from->rows, sizeof(uint) * c->rlen);
        } else {
            /* Adjust size of runmap memory */
            p2mem *mem = c->mem;
            int rlen = c->ymax - c->ymin + 2;
            int dlen = c->rows[rlen-1];
            if (rlen < c->rlen) 
                c->rows = p2mem_resize(mem, c->rows, rlen*sizeof(uint));
            c->rlen = rlen;
            if (dlen < c->dlen)
                c->data = p2mem_resize(mem, c->data, dlen);
            c->dlen = dlen;
            if (! (c->rows && c->data))
                return_VMerror;
        }
    }
    /* Delete operand runmaps */
    *out = c;
    if (a && (del & 1))
        runmap_free(a);
    if (b && (del & 2))
        runmap_free(b);
    return 0;
}


/* Check that there are remaining triples and obtain a run length.  */
#define rp_check_decode(runs, eruns, y, x, c, finished) { \
  if (runs < eruns && runs[0][0] == y) { c = !c; \
     x = (c ? (runs++)[0][2]+1 : runs[0][1]); } \
  else { x = ARCH_MAX_UINT; finished; } }

/* Pointer to triples {y,x1,x2} representing runs. */
typedef uint (*runptr)[3];

/* Compute bitwise OR of a runmap and an ordered array of triples <{y,x1,x2}>.
   Runmap A will be freed if <(del&1)> is non zero.  Resulting runmap is
   returned through pointer <out>. */
private int
runmap_or_runs(p2mem *mem, runmap *a, runptr runs, int nruns, 
	       int del, runmap **out)
{
    int y;
    int ymin, ymax, dlen;
    runmap *c = 0;
    byte *cdata;
    runptr eruns = runs + nruns;
    int code;
    /* determine ymin ymax */
    ymin = runs[0][0];
    if (a && a->ymin < ymin)  
        ymin = a->ymin;
    ymax = runs[nruns-1][0];
    if (a && a->ymax > ymax)
        ymax = a->ymax;
    /* allocate runmap */
    dlen = 6*nruns + 4*(ymax-ymin) + (a ? a->dlen : 0);
    code = runmap_alloc(mem, dlen, ymax-ymin+2, &c);
    if (code < 0) return code;
    c->xmin = (a ? a->xmin : ARCH_MAX_UINT);
    c->xmax = (a ? a->xmax : 0);
    cdata = c->data;
    /* loop over lines */
    for (y=ymin; y<=ymax; y++) {
        uint xmin = ARCH_MAX_UINT;
        uint xmax = 0;
        const byte *adata = 0;
        const byte *adataend = 0;
        /* access rows */
        if (a && y >= a->ymin && y <= a->ymax) {
            adata = a->data + a->rows[y - a->ymin];
            adataend = a->data + a->rows[y - a->ymin + 1];
        }
        /* process row */
        if (! (runs<eruns && runs[0][0]==y)) {
            /* nothing in b */
            if (adata < adataend) {
                memcpy(cdata, adata, adataend-adata);
                cdata += adataend - adata;
		xmin = a->xmin;
		xmax = a->xmax;
            }
        } else {
            /* general case */
            uint ax = 0;
            uint bx = 0;
            uint cx = 0;
            uint cx2 = 0;
            bool ac = true;
            bool bc = true;
            uint cc = false;
            uint lastcx = 0;
            uint lastcc = false;
            for(;;) {
                cx = cx2;
                if (ax <= cx)
		    rl_check_decode(adata, adataend, ax, ac, ac=false);
                if (bx <= cx)
		    rp_check_decode(runs, eruns, y, bx, bc, bc=false);
                cx2 = min(ax,bx);
                cc = ac || bc;
                if (cx2 >= ARCH_MAX_UINT)       
                    break;
                if (cc != lastcc) {             
                    if (lastcc) {
                        xmin = min(xmin, lastcx);
                        xmax = max(xmax, cx-1);
                    }
                    rl_encode(cdata, cx-lastcx); /*ENCODEC*/
                    lastcc = !lastcc;
                    lastcx = cx;
                }                               
            }                                   
            if (lastcc && cx>lastcx) {          
                xmin = min(xmin, lastcx);
                xmax = max(xmax, cx-1);
                rl_encode(cdata, cx-lastcx); /*ENCODEC*/
            }                                   
        }
        /* update bounding box */
	if (xmin <= xmax) {
	    c->ymin = min(c->ymin, y);
	    c->ymax = max(c->ymax, y);
	    c->xmin = min(c->xmin, xmin);
	    c->xmax = max(c->xmax, xmax);
	}
        /* fill rows array */
	if (c->ymin <= c->ymax)
	    c->rows[y - c->ymin + 1] = cdata - c->data;
        ASSERT(cdata <= c->data + c->dlen);
    }    
    /* finish */
    return runmap_finish(mem, c, a, 0, del, out);
}

/* Compute bitwise OR of two runmaps.  Runmap <a> will be freed if <(del&1)>
   is non zero.  Runmap <b> will be freed if <(del&2)> is non zero.  Resulting
   runmap is returned through pointer <out>. */
private int
runmap_or(p2mem *mem, runmap *a, runmap *b, 
	  int del, runmap **out)
{
    int y;
    int ymin, ymax, dlen;
    runmap *c = 0;
    byte *cdata;
    int code;
    /* determine ymin ymax dlen */
    if (! a) return runmap_finish(mem, b, a, b, del, out);
    if (! b) return runmap_finish(mem, a, a, b, del, out);
    ymin = min(a->ymin, b->ymin);
    ymax = max(a->ymax, b->ymax);
    dlen = a->dlen + b->dlen + 4*(ymax-ymin);
    /* allocate runmap */
    code = runmap_alloc(mem, dlen, ymax-ymin+2, &c);
    if (code < 0) return code;
    c->xmin = min(a->xmin, b->xmin);
    c->xmax = max(a->xmax, b->xmax);
    cdata = c->data;
    /* loop over lines */
    for (y=ymin; y<=ymax; y++) {
        uint xmin = ARCH_MAX_UINT;
        uint xmax = 0;
        byte *adata = 0;
        byte *adataend = 0;
        byte *bdata = 0;
        byte *bdataend = 0;
        /* access rows */
        if (y >= a->ymin && y <= a->ymax) {
            adata = a->data + a->rows[y - a->ymin];
            adataend = a->data + a->rows[y - a->ymin + 1];
        }
        if (y >= b->ymin && y <= b->ymax) {
            bdata = b->data + b->rows[y - b->ymin];
            bdataend = b->data + b->rows[y - b->ymin + 1];
        }
        /* process row */
        if (bdata >= bdataend) {
            /* nothing in b */
            if (adata < adataend) {
                memcpy(cdata, adata, adataend-adata);
                cdata += adataend - adata;
		xmin = a->xmin;
		xmax = a->xmax;
            }
        } else if (adata >= adataend) {
            /* nothing in a */
            memcpy(cdata, bdata, bdataend-bdata);
            cdata += bdataend - bdata;
	    xmin = b->xmin;
	    xmax = b->xmax;
        } else {
            /* general case */
            uint ax, bx, cx, cx2, lastcx;
	    bool ac, bc, cc, lastcc;
            ax = bx = cx = cx2 = lastcx = 0;
	    ac = bc = cc = lastcc = false;
	    rl_decode(adata, ax);
	    rl_decode(bdata, bx);
	    /* loop */
            for(;;) {
                cx = cx2;
                if (ax <= cx) 
		    rl_check_decode(adata, adataend, ax, ac, ac=false);
                if (bx <= cx) 
		    rl_check_decode(bdata, bdataend, bx, bc, bc=false);
                cx2 = min(ax,bx);
                cc = ac || bc;
                if (cx2 >= ARCH_MAX_UINT)       
                    break;
                if (cc != lastcc) {             
                    if (lastcc) {
                        xmin = min(xmin, lastcx);
                        xmax = max(xmax, cx-1);
                    }
                    rl_encode(cdata, cx-lastcx); /*ENCODEC*/
                    lastcc = !lastcc;
                    lastcx = cx;
                }                               
            }                                   
            if (lastcc && cx>lastcx) {          
                xmin = min(xmin, lastcx);
                xmax = max(xmax, cx-1);
                rl_encode(cdata, cx-lastcx); /*ENCODEC*/
            }                                   
        }
        /* update bounding box */
	if (xmin <= xmax) {
	    c->ymin = min(c->ymin, y);
	    c->ymax = max(c->ymax, y);
	    c->xmin = min(c->xmin, xmin);
	    c->xmax = max(c->xmax, xmax);
	}
        /* fill rows array */
	if (c->ymin <= c->ymax)
	    c->rows[y - c->ymin + 1] = cdata - c->data;
        ASSERT(cdata <= c->data + c->dlen);
    }
    /* finish */
    return runmap_finish(mem, c, a, b, del, out);
}

/* Compute bitwise AND of two runmaps.  Runmap <a> will be freed if <(del&1)>
   is non zero.  Runmap B will be freed if <(del&2)> is non zero.  Resulting
   runmap is returned through pointer <out>. */
private int
runmap_and(p2mem *mem, runmap *a, runmap *b, 
	   int del, runmap **out)
{
    int y;
    int ymin, ymax, dlen;
    runmap *c = 0;
    byte *cdata;
    int code;
    /* determine ymin ymax dlen */
    if (!a || !b) 
        return runmap_finish(mem, 0, a, b, del, out);
    ymin = max(a->ymin, b->ymin);
    ymax = min(a->ymax, b->ymax);
    if (ymin > ymax || a->xmin > b->xmax || a->xmax < b->xmin) 
        return runmap_finish(mem, 0, a, b, del, out);        
    dlen = 4 * (ymax - ymin);
    dlen += a->rows[ymax + 1 - a->ymin] - a->rows[ymin - a->ymin];
    dlen += b->rows[ymax + 1 - b->ymin] - b->rows[ymin - b->ymin];
    /* allocate runmap */
    code = runmap_alloc(mem, dlen, ymax-ymin+2, &c);
    if (code < 0) return code;
    c->xmin = ARCH_MAX_UINT;
    c->xmax = 0;
    cdata = c->data;
    /* loop over lines */
    for (y=ymin; y<=ymax; y++) {
        uint xmin = ARCH_MAX_UINT;
        uint xmax = 0;
        /* access rows */
        byte *adata = a->data + a->rows[y - a->ymin];
        byte *adataend = a->data + a->rows[y - a->ymin + 1];
        byte *bdata = b->data + b->rows[y - b->ymin];
        byte *bdataend = b->data + b->rows[y - b->ymin + 1];
        /* process row */
        if (bdata < bdataend && adata < adataend) {
            /* something in both rows a and b */
            uint ax = 0;
            uint bx = 0;
            uint cx = 0;
            uint cx2 = 0;
            bool ac = true;
            bool bc = true;
            uint cc = false;
            uint lastcx = 0;
            uint lastcc = false;
            for(;;) {
                cx = cx2;
                if (ax <= cx) 
		    rl_check_decode(adata, adataend, ax, ac, ac=false);
                if (bx <= cx) 
		    rl_check_decode(bdata, bdataend, bx, bc, bc=false);
                cx2 = min(ax,bx);
                cc = ac && bc;
                if (ax >= ARCH_MAX_UINT || bx >= ARCH_MAX_UINT)       
                    break;
                if (cc != lastcc) {             
                    if (lastcc) {
                        xmin = min(xmin, lastcx);
                        xmax = max(xmax, cx-1);
                    }
                    rl_encode(cdata, cx-lastcx); /*ENCODEC*/
                    lastcc = !lastcc;
                    lastcx = cx;
                }                               
            }                                   
            if (lastcc && cx>lastcx) {          
                xmin = min(xmin, lastcx);
                xmax = max(xmax, cx-1);
                rl_encode(cdata, cx-lastcx); /*ENCODEC*/
            }                                   
        }
        /* update bounding box */
	if (xmin <= xmax) {
	    c->ymin = min(c->ymin, y);
	    c->ymax = max(c->ymax, y);
	    c->xmin = min(c->xmin, xmin);
	    c->xmax = max(c->xmax, xmax);
	}
        /* fill rows array */
	if (c->ymin <= c->ymax)
	    c->rows[y - c->ymin + 1] = cdata - c->data;
        ASSERT(cdata <= c->data + c->dlen);
    }
    /* finish */
    return runmap_finish(mem, c, a, b, del, out);
}

/* Compute bitwise SUBTRACT of two runmaps.  Runmap <a> will be freed if
   <(del&1)> is non zero.  Runmap <b> will be freed if <(del&2)> is non zero.
   Resulting runmap is returned through pointer <out>. */
private int
runmap_andnot(p2mem *mem, runmap *a, runmap *b, 
	      int del, runmap **out)
{
    int y;
    int ymin, ymax, dlen;
    runmap *c = 0;
    byte *cdata;
    int code;
    /* determine ymin ymax dlen */
    if (! a) return runmap_finish(mem, 0, a, b, del, out);
    if (! b) return runmap_finish(mem, a, a, b, del, out);
    ymin = max(a->ymin, b->ymin);
    ymax = min(a->ymax, b->ymax);
    if (ymin > ymax || b->xmin > a->xmax || b->xmax < a->xmin  )
        return runmap_finish(mem, a, a, b, del, out);        
    dlen = a->dlen + 4 * (ymax - ymin);
    dlen += b->rows[ymax + 1 - b->ymin] - b->rows[ymin - b->ymin];
    ymin = a->ymin;
    ymax = a->ymax;
    /* allocate runmap */
    code = runmap_alloc(mem, dlen, ymax-ymin+2, &c);
    if (code < 0) return code;
    c->xmin = ARCH_MAX_UINT;
    c->xmax = 0;
    cdata = c->data;
    /* loop over lines */
    for (y=ymin; y<=ymax; y++) {
        uint xmin = ARCH_MAX_UINT;
        uint xmax = 0;
        /* access rows */
        byte *adata = a->data + a->rows[y - a->ymin];
        byte *adataend = adataend = a->data + a->rows[y - a->ymin + 1];
        byte *bdata = 0;
        byte *bdataend = 0;
        if (y >= b->ymin && y <= b->ymax) {
            bdata = b->data + b->rows[y - b->ymin];
            bdataend = b->data + b->rows[y - b->ymin + 1];
        }
        /* process */
        if (bdata >= bdataend ) {
            /* nothing in b */
            if (adata < adataend) {
                uint ax = 0;
                memcpy(cdata, adata, adataend-adata);
                cdata += adataend - adata;
		rl_decode(adata, ax);
		xmin = xmax = ax;
		for (;;) {
		    if (adata>=adataend) break;
                    rl_decode(adata, ax);
		    xmax = ax-1;
		    if (adata>=adataend) break;
                    rl_decode(adata, ax);
		}
            }
        } else if (adata >= adataend) {
            /* nothing in a */
        } else {
            /* general case */
            uint ax = 0;
            uint bx = 0;
            uint cx = 0;
            uint cx2 = 0;
            bool ac = true;
            bool bc = true;
            uint cc = false;
            uint lastcx = 0;
            uint lastcc = false;
            for(;;) {
                cx = cx2;
                if (ax <= cx) 
		    rl_check_decode(adata, adataend, ax, ac, ac=false);
                if (bx <= cx) 
		    rl_check_decode(bdata, bdataend, bx, bc, bc=false);
                cx2 = min(ax,bx);
                cc = ac && !bc;
                if (ax >= ARCH_MAX_UINT)
                    break;
                if (cc != lastcc) {             
                    if (lastcc) {
                        xmin = min(xmin, lastcx);
                        xmax = max(xmax, cx-1);
                    }
                    rl_encode(cdata, cx-lastcx); /*ENCODEC*/
                    lastcc = !lastcc;
                    lastcx = cx;
                }                               
            }                                   
            if (lastcc && cx>lastcx) {          
                xmin = min(xmin, lastcx);
                xmax = max(xmax, cx-1);
                rl_encode(cdata, cx-lastcx); /*ENCODEC*/
            }                                   
        }
        /* update bounding box */
	if (xmin <= xmax) {
	    c->ymin = min(c->ymin, y);
	    c->ymax = max(c->ymax, y);
	    c->xmin = min(c->xmin, xmin);
	    c->xmax = max(c->xmax, xmax);
	}
        /* fill rows array */
	if (c->ymin <= c->ymax)
	    c->rows[y - c->ymin + 1] = cdata - c->data;
        ASSERT(cdata <= c->data + c->dlen);
    }
    /* finish */
    return runmap_finish(mem, c, a, b, del, out);
}

/* Compute area and perimeter of runmap. */
private void
runmap_info(runmap *r, int *parea, int *pperimeter)
{
    int y, nrun, area, canc;
    byte *pdata, *pdataend;
    byte *data, *dataend;
    /* check */
    nrun = area = canc = 0;
    if (r && r->rows && r->data && r->ymax>=r->ymin) {
        /* haffner's computation of perimeter */
        pdataend = r->data;
        dataend = r->data;
        for (y=r->ymin; y<=r->ymax; y++) {
            int x = 0;
            int x1 = 0;
            int px = 0;
            int px1 = 0;
            pdata = pdataend;
            pdataend = data = dataend;
            dataend = r->data + r->rows[y - r->ymin + 1];
            while (data < dataend) {
                rl_decode(data, x);
                if (data >= dataend) 
                    break;
                x1 = x;
                rl_decode(data, x);
                nrun += 1;
                area += x - x1;
                while (px1 < x)
                {
                    if (px > x1) { 
                        canc += min(x,px) - max(x1,px1);
                        if (px >= x)
                            break;
                    }
                    if (pdata >= pdataend) 
                        break;
                    rl_decode(pdata, px);
                    if (pdata >= pdataend) 
                        break;
                    px1 = px;
                    rl_decode(pdata, px);
                }
            }
        }
    }
    if (parea)
        *parea = area;
    if (pperimeter)
        *pperimeter = 2*(nrun+area-canc);
}

/* Save a runmap as PBM */
private int
runmap_save(runmap *r, FILE *f, 
            uint xmin, uint xmax, uint ymin, uint ymax, 
            const char *comment)
{
    uint w = xmax - xmin + 1;
    uint h = ymax - ymin + 1;
    uint rowsize;
    byte *buf;
    uint y;
    rowsize = (w+7)/8;
    if (w<=0 || h<=0) return_error(gs_error_limitcheck);
    buf = p2mem_parent_alloc(rowsize);
    if (! buf) return_VMerror;
    /* PBM files begin with a header line "P4 <width> <height>\n"
       indicating the bitmap size.  The rest of the file encodes 
       the rows as a packed bitmap.  Each row occupies an integral
       number of byte.  The most significant bit of a byte represents
       the leftmost pixel coded by that byte. */
    fprintf(f,"P4\n");
    if (comment) 
        fprintf(f, "# %s\n", comment);
    fprintf(f,"%d %d\n", w, h);
    for (y=ymin; y<=ymax; y++) {
        uint x = xmin;
        uint nx = 0;
        int z = 0;
        int m = 0x80;
        byte *ptr = buf;
        bool c = true;
        byte *data = 0;
        byte *dataend = 0;
        if (r && y>=r->ymin && y<=r->ymax) {
            data = r->data + r->rows[y - r->ymin];
            dataend = r->data + r->rows[y - r->ymin + 1];
        }
        while (nx<=xmax) {
            rl_check_decode(data, dataend, nx, c, c=false);
            nx = min(xmax+1, nx);
            while (x < nx) {
                if (c) {
                    while (x<nx && m) { x++; z|=m; m>>= 1; } 
                } else {
                    while (x<nx && m) { x++; m>>= 1; } 
                }
                if (! m) {
                    *ptr++ = z;
                    z = (c ? 0xff : 0x00);
                    while (x+8<=nx) { x+=8; *ptr++ = z; }
                    z = 0x00;
                    m = 0x80;
                }
            }
        }
        if (m < 0x80)
            *ptr++ = z;
        ASSERT(ptr == buf+rowsize);
        if (fwrite(buf, 1, rowsize, f) != rowsize) {
            p2mem_parent_free(buf);
            return_error(gs_error_ioerror);
        }
    }
    p2mem_parent_free(buf);
    return 0;
}

/* Draws a runmap with a specified color into a raster pixmap.  Raster covers
   rectangle <cwxch+cx+cy>. Pointer <base> indicates position <(cx,cy)>.
   Integer <sraster> represents the pointer difference from one line to the
   next. */
private void
runmap_play(runmap *r, gx_color_index color,
            byte *base, int sraster, uint cx, uint cy, uint cw, uint ch)
{
    byte c[3];
    uint cx2 = cx + cw;
    uint cy2 = cy + ch;
    c[0] = color >> 16;
    c[1] = color >> 8;
    c[2] = color;
    if (r) {
        uint y1 = max(cy, r->ymin);
        uint y2 = min(cy2, r->ymax + 1);
        byte *p = base + (y1 - cy) * sraster - cx;
        while (y1 < y2) {
            byte *data = r->data + r->rows[y1 - r->ymin];
            byte *dataend = r->data + r->rows[y1 - r->ymin + 1];
            if (data < dataend) {
                uint x = 0;
                while (data < dataend) {
                    rl_decode(data, x);
                    if (data < dataend && x < cx2) {
                        uint x1 = x;
                        rl_decode(data, x);
                        if (x > cx) {
                            byte *p1 = p + 3 * max(x1, cx);
                            byte *p2 = p + 3 * min(x, cx2);
                            while (p1 < p2) {
                                *p1++ = c[0];
                                *p1++ = c[1];
                                *p1++ = c[2];
                            }
                        }
                    }
                }
            }
            y1 += 1;
            p += sraster;
        }
    }    
}



/* ======================================
            R U N M A P B L D
   ====================================== */

/* Runmap builders allow the efficient construction of rmap using arbitrarily
   ordered black runs.  Each black run is added using runmapbld_write() and
   extends from column x1 to x2 (inclusive) on row y.  Runs are stored in a
   temporary buffer and ORed into the current runmap.  The final runmap is
   computed when runmapbld_close() is called. */

/* The runmapbld data structure */
typedef struct runmapbld_s {
    p2mem *mem;
    runmap *rmap;
    runptr runs;
    uint nruns;
    uint nrunsmax;
    bool needsort;
} runmapbld;

/* Recommended buffer size for a runmapbld */
#define runmapbld_buffer_size 5000

/* Create a runmap builder. */
private int
runmapbld_open(p2mem *mem, uint nrunsmax, runmapbld **out)
{
    runmapbld *c = p2mem_alloc(mem, sizeof(runmapbld));
    runptr r = p2mem_alloc(mem, sizeof(uint)*3*nrunsmax);
    if (! (c && r)) return_VMerror;
    c->mem = mem;
    c->rmap = 0;
    c->runs = r;
    c->nruns = 0;
    c->nrunsmax = nrunsmax;
    c->needsort = false;
    *out = c;
    return 0;
}

/* Internal: qsort callback for ordering runs. */
private int 
runmapbld_sortsub(const void *av, const void *bv)
{
    const int *a = (const int *)av;
    const int *b = (const int *)bv;
    return a[0]!=b[0] ? a[0]-b[0] : a[1]!=b[1] ? a[1]-b[1] : a[2]-b[2];
}

/* Compute up-to-date runmap. */
private int
runmapbld_flush(runmapbld *bld)
{
    int code = 0;
    if (bld->needsort && bld->nruns>1) {
        int i,j;
        int nruns = bld->nruns;
        runptr runs = bld->runs;
        qsort(runs, nruns, sizeof(int[3]), runmapbld_sortsub);
        for (i=1,j=0; i<nruns; i++) {
            if (runs[i][0]==runs[j][0] && runs[i][1]<=runs[j][2]+1) {
                runs[j][2] = max(runs[i][2],runs[j][2]);
            } else {
                j += 1; 
                runs[j][0]=runs[i][0]; 
                runs[j][1]=runs[i][1]; 
                runs[j][2]=runs[i][2]; 
            }
        }
        bld->nruns = j+1;
    }
    if (bld->nruns > 0)
        code = runmap_or_runs(bld->mem, bld->rmap, bld->runs, bld->nruns, 
			      1, &bld->rmap);
    bld->nruns = 0;
    bld->needsort = false;
    return code;
}

/* Add run {y,x1,x2} to the runmap builder. */
private int 
runmapbld_write(runmapbld *bld, int y, int x1, int x2)
{
    uint *lastrun;
    if (bld->nruns > 0) {
	lastrun = bld->runs[bld->nruns-1];
        if (lastrun[0] == y) {
            if (x1 < lastrun[1])
                bld->needsort = true;
            if (lastrun[1]<=x2+1 && lastrun[2]+1>=x1) {
                lastrun[1] = min(x1, lastrun[1]);
                lastrun[2] = max(x2, lastrun[2]);
                return 0;
            }
        } else if (y<lastrun[0])
            bld->needsort = true;
    }            
    if (bld->nruns >= bld->nrunsmax) {
        int code = runmapbld_flush(bld);
        if (code < 0) return code;
    }
    lastrun = bld->runs[bld->nruns++];
    lastrun[0] = y;
    lastrun[1] = x1;
    lastrun[2] = x2;
    return 0;
}

/* Destroy runmap builder, release memory, and return runmap. */
private int 
runmapbld_close(runmapbld *bld, runmap **out)
{
    if (bld) {
        if (out) {
            int code = runmapbld_flush(bld);
            if (code < 0) return code;
            *out = bld->rmap;
            bld->rmap = 0;
        }
        runmap_free(bld->rmap);
        p2mem_free(bld->mem, bld->runs);
        p2mem_free(bld->mem, bld);
        return 0;
    }
    if (out) 
        *out = 0;
    return 0;
}

/* Accumulate a rectangle into a runmapbld */
private int
runmapbld_put_rect(runmapbld *bld,
                   int xx, int yy, int w, int h)
{
    int y;
    int code = 0;
    for (y=0; y<h; y++)
        if ((code = runmapbld_write(bld, yy+y, xx, xx+w-1)) < 0)
            return code;
    return 0;
}

/* Accumulate a bitmap into a runmapbld.
   Arguments <base>, <sourcex>, <sraster>, <xx>, <yy>, <w> and <h>
   are as in ghostscript.  Setting <xor> to <0> copies all black pixels
   into the runmapbld.  Setting <xor> to <0xff> copies all white pixels
   and produces a reverse video effect. */
private int
runmapbld_put_bitmap(runmapbld *bld,
                     const byte *base, int sourcex, int sraster, int xor,
                     int xx, int yy, int w, int h)
{
    int y;
    int code = 0;
    int basemask = 0x80 >> (sourcex & 0x7);
    base = base + (sourcex >> 3);
    for (y=0; y<h; y++, base+=sraster) {
        const byte *row = base;
        int z = (*row++) ^ xor;
        int m = basemask;
        int x1 = 0;
        int x = 0;
        while (x < w) {
            while (x < w) {
                if (m) {
                    if (z & m)
                        break;
                    x += 1;
                    m >>= 1;
                } else {
                    z = (*row++) ^ xor;
                    while (z==0x00 && x+7<w) {
                        z = (*row++) ^ xor;
                        x += 8;
                    }
                    m = 0x80;
                }
            }
            x1 = x;
            while (x < w) {
                if (m) {
                    if (! (z & m))
                        break;
                    x += 1;
                    m >>= 1;
                } else {
                    z = (*row++) ^ xor;
                    while (z==0xff && x+7<w) {
                        z = (*row++) ^ xor;
                        x += 8;
                    }
                    m = 0x80;
                }
            }
            if (x > x1)
                if ((code = runmapbld_write(bld, yy+y, xx+x1, xx+x-1)) < 0)
                    return code;
        }
    }
    return 0;
}




/* ======================================
           H A S H   T A B L E
      A N D   A P P L I C A T I O N S
   ====================================== */

/* Generic hash table.  Each hash table entry starts with a user
   defined payload structure whose first component is the key.  In the
   simplest case, you just specify the size of the payload when 
   creating the hash table.  You might then override procedure pointers 
   for comparing keys, copying keys, and releasing extra memory associated 
   with the payload structure.  Macro htable_lookup() gives a pointer to 
   the payload structure associated with a particular key.  You can then 
   consult or change any additional field in the payload structure.  
   A null pointer is returned if the key is not found and argument 
   create is false.  */

typedef struct htable_s {
    p2mem *mem;     /* memory allocator */
    int psize;      /* payload size */
    int xsize;      /* payload size rounded for alignment */
    int nelems;     /* number of elements */
    int nbuckets;   /* number of hash buckets */
    void **buckets; /* pointer to entries in each bucket */
    /* Procedure pointer you can override */
    bool (*key_equal)(const void *key1, const void *key2);
    void (*key_copy)(p2mem *mem, void *dst, const void *src);
    void (*payload_release)(p2mem *mem, void *payload);
} htable;

/* Internal: hash table entry fields following the payload */
typedef struct htablelink_s {
    uint hash;
    void *next;
} htablelink;

/* Internal: change the number of buckets. */
private int
htable_reorganize(htable *h, int nbuckets)
{
    int b;
    int oldnbuckets = h->nbuckets;
    void **oldbuckets = h->buckets;
    void **buckets = p2mem_alloc(h->mem, nbuckets*sizeof(void*));
    if (! buckets) return_VMerror;
    memset(buckets, 0, nbuckets*sizeof(void*));
    h->nbuckets = nbuckets;
    h->buckets = buckets;
    for (b=0; b<oldnbuckets; b++) {
	void *e = oldbuckets[b];
	while (e) {
	    htablelink *f = (htablelink*)((char*)e + h->xsize);
	    void *enext = f->next;
	    int nb = f->hash % nbuckets;
	    f->next = buckets[nb];
	    buckets[nb] = e;
	    e = enext;
	}
    }
    p2mem_free(h->mem, oldbuckets);
    return 0;
}

/* Create a hash table. */
private htable *
htable_alloc(p2mem *mem, int payloadsize)
{
    int a;
    htable *h = p2mem_alloc(mem, sizeof(htable));
    if (h) {
	memset(h, 0, sizeof(htable));
	h->mem = mem;
	h->psize = payloadsize;
	/* compute alignment of htablelink */
	a = (int)((size_t)&(((struct{char pad; htablelink e;}*)0)->e));
	/* compute offset of htablelink w.r.t. payload */
	h->xsize = (int)((payloadsize + a - 1) / a) * a;
	/* setup initial bucket table */
	if (htable_reorganize(h, 17) >= 0)
	    return h;
    }
    /* allocation failed */
    p2mem_free(mem, h);
    return 0;
}

/* Destroy a hash table */
private void
htable_free(htable *h) 
{
    if (h) {
	int b;
	for (b=0; b<h->nbuckets; b++) {
	    void *e = h->buckets[b];
	    while (e) {
		htablelink *f = (htablelink*)((char*)e + h->xsize);
		void *enext = f->next;
		if (h->payload_release) 
		    h->payload_release(h->mem, e);
		p2mem_free(h->mem, e);
		e = enext;
	    }
	}
	p2mem_free(h->mem, h->buckets);
	p2mem_free(h->mem, h);
    }
}

/* Htable lookup function. Returns a payload pointer if a matching entry is
   found, otherwise creates a new entry or return a null pointer according
   to flag ``create''.  Note the the type of argument ``pkey'' *must* be
   ``pointer to the key type''.  Note that argument hash must be the hashcode
   for the key pointed to by ``pkey''. */
#define htable_lookup(htbl, pkey, hash, create) \
   (htable_lookup_sub(htbl, pkey, sizeof(*pkey), hash, create))

/* Internal: helper for htable_lookup() */
private void *
htable_new_entry(htable *h, const void *key, 
                 int keysize, uint hash, int bucket)
{
    int esize = h->xsize + sizeof(htablelink);
    htablelink *f;
    void *e;
    if (h->nelems*3 > h->nbuckets*2)
        if (htable_reorganize(h, 2 * h->nbuckets - 1) >= 0)  
            bucket = hash % h->nbuckets;
    if (! (e = p2mem_alloc(h->mem, esize)))
        return 0;
    memset(e, 0, esize);
    if (h->key_copy) 
        (*h->key_copy)(h->mem, e, key);
    else 
        memcpy(e, key, keysize);
    f = ((htablelink*)((char*)e + h->xsize));
    f->hash = hash;
    f->next = h->buckets[bucket];
    h->buckets[bucket] = e;
    h->nelems += 1;
    return e;
}

/* Internal: helper for htable_lookup() */
private void *
htable_lookup_sub(htable *h, const void *key, 
                  int keysize, uint hash, bool create)
{
    int bucket = hash % h->nbuckets;
    void *e = h->buckets[bucket];
    if (h->key_equal) {
	while (e && !h->key_equal(e, key))
	    e = ((htablelink*)((char*)e + h->xsize))->next;
    } else {
	while (e && memcmp(e, key, keysize))
	    e = ((htablelink*)((char*)e + h->xsize))->next;
    }
    if (!e && create)
        e = htable_new_entry(h, key, keysize, hash, bucket);
    return e;
}

/* Iterate on hash table entries. */
#define htable_begin_loop(pptr, h) {\
   int _b_; \
   for ( _b_ = 0; _b_ < (h)->nbuckets; _b_++ ) { \
      pptr = (h)->buckets[_b_]; \
      while (pptr) { 

#define htable_end_loop(pptr, h) \
        pptr = ((htablelink*)((char*)(pptr) + (h)->xsize))->next; } \
      } }





/* ======================================
    C O L O R   Q U A N T I Z A T I O N
   ====================================== */


/* Payload for color histogram hash table */
typedef struct colordata_s {
    gx_color_index color; 
    uint w; 
} colordata;

/* Hash function for colors. */
#define hash_rgb(rgb) \
   (((rgb) & 0xffffff)^(((rgb) & 0xffff)<<4))

/* Allocate a color hash table */
#define colorhash_alloc(mem) \
   (htable_alloc(mem, sizeof(colordata)))

/* Free color hash table */
#define colorhash_free(colorhash) \
   htable_free(colorhash)

/* Add contribution to the histogram */
private int
colorhash_add(htable *colorhash, gx_color_index color, uint w)
{
    colordata *cd;
    if (! (cd = htable_lookup(colorhash, &color, hash_rgb(color), true)))
        return_VMerror;
    cd->w += w;
    return 0;
}

/* Reference for the median color quantization algorithm: Paul
   Heckbert: "Color Image Quantization for Frame Buffer Display",
   SIGGRAPH '82 Proceedings, page 297.  Also used in ppmquant. */

/* Internal: Colorspace box for quantizer */
typedef struct colorbox_s { 
    struct colorbox_s *next;
    colordata *cdata;
    int ndata; 
    int w; 
} colorbox;

/* Internal: sorting direction mask */
private gx_color_index colorcomp_dir; 

/* Internal: sorting subroutine */
private int
colorcomp_sub(const void *pa, const void *pb)
{
    gx_color_index a = colorcomp_dir & *(const gx_color_index*)pa;
    gx_color_index b = colorcomp_dir & *(const gx_color_index*)pb;
    return a - b;
}

/* Reduce number of colors.  Input <colorhash> is a color hash table
   containing a color histogram.  This function will allocate a palette array
   <*ppalette> and return its size as <*psize>.  The weights in the color hash
   table will be replaced by the palette index.  */
private int 
color_quantize(p2mem *mem, htable *colorhash, uint maxcolors,
               gx_color_index **ppalette, uint *psize)
{ 
    /* Variables */
    int code;
    gx_color_index *palette = 0;
    colordata *cdata = 0;
    colorbox  *cbox = 0;
    colorbox *cb;
    colorbox **pcb;
    colordata *cd;
    int ncolors;
    int w;
    int i;
    /* Copy histogram in cdata array */
    cdata = p2mem_alloc(mem, sizeof(colordata)*colorhash->nelems);
    if (!cdata) goto xitmem;
    w = 0;
    ncolors = 0;
    htable_begin_loop(cd, colorhash) {
	cdata[ncolors++] = cd[0];
	w += cd->w;
    } htable_end_loop(cd, colorhash);
    /* Create first quantization box */
    cbox = p2mem_alloc(mem, sizeof(colorbox));
    if (!cbox) goto xitmem;
    cbox->next = 0;
    cbox->cdata = cdata;
    cbox->ndata = ncolors;
    cbox->w = w;
    ncolors = 1;
    /* Median color quantization */
    while (ncolors < maxcolors) {
	int rl,bl,gl;
	byte pmin[3], pmax[3];
	colorbox *ncb;
	/* Find suitable box */
	for (pcb=&cbox; *pcb; pcb=&((*pcb)->next))
	    if ((*pcb)->ndata>=2)
		break;
	/* Unlink box or exit loop */
	if (! (cb = *pcb)) break;
	*pcb = cb->next;
	/* Compute box boundaries */
	pmin[0]=pmin[1]=pmin[2]=0xff;
	pmax[0]=pmax[1]=pmax[2]=0x00;
	cd = cb->cdata;
	for (i=0; i<cb->ndata; i++) {
	    gx_color_index color = cb->cdata[i].color;
	    byte r=(color>>16), g=(color>>8), b=(color);
	    pmin[0] = min(pmin[0],r);
	    pmin[1] = min(pmin[1],g);
	    pmin[2] = min(pmin[2],b);
	    pmax[0] = max(pmax[0],r);
	    pmax[1] = max(pmax[1],g);
	    pmax[2] = max(pmax[2],b);
	}
	rl = pmax[0] - pmin[0];
	gl = pmax[1] - pmin[1];
	bl = pmax[2] - pmin[2];
	/* Determine split direction */
	if (gl>=bl && gl>=rl)
	    colorcomp_dir = 0x00ff00;
	else if (rl>=gl && rl>=bl)
	    colorcomp_dir = 0xff0000;
	else 
	    colorcomp_dir = 0x0000ff;
	/* Sort colordata and find median */
	qsort(cb->cdata, cb->ndata, sizeof(colordata), colorcomp_sub);
	for (i=w=0; i<cb->ndata-1 && w+w<cb->w; i++) 
	    w += cb->cdata[i].w;
	/* Prepare new boxes */
	ncolors += 1;
	ncb = p2mem_alloc(mem, sizeof(colorbox));
	if (!ncb) goto xitmem;
	ncb->cdata = cb->cdata + i;
	ncb->ndata = cb->ndata - i;
	ncb->w = cb->w - w;
	cb->ndata = i;
	cb->w = w;
	/* Insert cb at proper location */
	pcb = &cbox;
	while (*pcb && (*pcb)->w >= cb->w)
	    pcb = &((*pcb)->next);
	cb->next = *pcb;
	*pcb = cb;
	/* Insert ncb at proper location */
	pcb = &cbox;
	while (*pcb && (*pcb)->w >= ncb->w)
	    pcb = &((*pcb)->next);
	ncb->next = *pcb;
	*pcb = ncb;
    }
    /* Allocate and fill palette */
    palette = p2mem_alloc(mem, ncolors*sizeof(gx_color_index));
    if (! palette) goto xitmem;
    *ppalette = palette;
    *psize = ncolors;
    ncolors = 0;
    for (cb=cbox; cb; cb=cb->next) {
	/* Compute mean color */
	byte r, g, b;
	double rsum, gsum, bsum;
	rsum = gsum = bsum = 0;
	for (i=0; i<cb->ndata; i++) {
	    gx_color_index color = cb->cdata[i].color;
	    w = cb->cdata[i].w;
	    rsum += (double)w * (double)(0xff & (color >> 16));
	    gsum += (double)w * (double)(0xff & (color >> 8));
	    bsum += (double)w * (double)(0xff & color);
	    cd = htable_lookup(colorhash, &color, hash_rgb(color), false);
	    cd->w = ncolors; /* Color hash now gives palette indice */
	}
	r = max(0, min(255, (int)(rsum/cb->w)));
	g = max(0, min(255, (int)(gsum/cb->w)));
	b = max(0, min(255, (int)(bsum/cb->w)));
	palette[ncolors] = (r<<16)|(g<<8)|(b);
	ncolors += 1;
    }
    /* Success */
    code = 0;
    goto xit;
 xitmem:
    /* Memory error */
    code = gs_note_error(gs_error_VMerror);
 xit:
    /* Free everything and return */
    p2mem_free(mem, cdata);
    while (cbox) {
	cb = cbox->next;
	p2mem_free(mem, cbox);
	cbox = cb;
    }
    return code;
}



/* ======================================
               C H R O M E
   ====================================== */

/* Some image components do not have a single color.  These components are
   normally classified as background.  In order to generate the background
   image, we need to record additional information besides the shape of the
   object.  This is stored as a sequence of rendering opcodes that can be
   replayed at later times.  Chrome data is moslty accessed sequentially.  It
   is stored in a linked list of data blocks and accessed like a file.  Using
   the ghostscript clist infrastructure might have been smarter. */

/* Internal: length of a chrome block */
#define chrome_block_len (p2mem_max - 8)

/* Internal: chrome blocks */
typedef struct chrome_block_s {
    struct chrome_block_s *next;
    byte data[chrome_block_len];
} chrome_block;

/* A chrome object represent a collection of images represented by a sequence
   of rendering opcodes that can be replayed at later times.  Each image is
   identified by a chrome_pos object that defines its offset and length in the
   chrome object.  */
typedef struct chrome_s {
    p2mem *mem;            /* Block memory */
    chrome_block *first;   /* Block chain  */
    chrome_block *last;    /* Current block for writing */
    byte *wpos;            /* Current pointer for writing */
    byte *wend;            /* Pointer to end of current block */
    byte *wlastrect;       /* Position following last rectangle opcode */
    int x, sx, y, sy;      /* Opcode context */
    int w, h;              /* Opcode context */
    gx_color_index c0, c1; /* Opcode context */
} chrome;

/* Internal: reset opcode context */
#define chrome_reset_ctx(ctx) { \
    (ctx)->x = (ctx)->sx = (ctx)->y = (ctx)->sy = 0; \
    (ctx)->w = (ctx)->h = 1; (ctx)->c0 = (ctx)->c1 = 0; }

/* Allocate chrome */
private chrome *
chrome_create(void)
{
    p2mem *nmem = p2mem_create();
    chrome *cr = 0;
    if (nmem) {
        cr = p2mem_alloc(nmem, sizeof(chrome));
        if (cr) {
            memset(cr, 0, sizeof(chrome));
            cr->mem = nmem;
            chrome_reset_ctx(cr);
            return cr;
        }
        p2mem_destroy(nmem);
    }
    return 0;
}

/* Destroy chrome */
private void
chrome_destroy(chrome *cr)
{
    if (cr != 0) 
        p2mem_destroy(cr->mem);
}

/* Data structure recording current chrome position */
typedef struct chrome_pos_s {
    chrome_block *block;  /* Stored block pointer */
    byte *pos;            /* Stored position pointer */
    int len;              /* Length of chrome data (-1 if unknown) */
} chrome_pos;

/* Create a chrome_pos indentifying the next image in the chrome object. */
private chrome_pos *
chrome_getpos(chrome *cr)
{
    /* Allocation is controlled by the chrome object */
    chrome_pos *cp = p2mem_alloc(cr->mem, sizeof(chrome_pos));
    if (!cp) return 0;
    cp->block = cr->last;
    cp->pos = cr->wpos;
    cp->len = -1;
    return cp;
}

/* Internal forward: Close chrome_pos object */
private int chrome_put_end(chrome *cr);

/* Terminate the current image in the chrome object and 
   update the length in the corresponding chrome_pos object */
private void
chrome_closepos(chrome *cr, chrome_pos *cpos)
{
    /* Insert an end-of-object opcode */
    chrome_put_end(cr);
    /* Fix block and pos */
    if (!cpos->block) {
        cpos->block = cr->first;
        cpos->pos = cr->first->data;
    }
    /* Compute length */
    if (cpos->block == cr->last) {
        cpos->len = cr->wpos - cpos->pos;
    } else {
        chrome_block *b = cpos->block;
        cpos->len = chrome_block_len - (cpos->pos - b->data);
        while ((b = b->next) != cr->last)
            cpos->len += chrome_block_len;
        cpos->len += cr->wpos - b->data;
    }
}

/* Internal: opcodes */
enum chrome_opcodes {
    /* minor opcodes */
    op_end          = 0x01,  /* -- end of object       */
    op_pixmap       = 0x02,  /* -- followed by pixmap data */
    op_set_c0       = 0x03,  /* -- followed by r,g,b */
    op_set_c1       = 0x04,  /* -- followed by r,g,b */
    op_rest_x       = 0x05,  /* -- x = sx */
    op_rest_y       = 0x06,  /* -- y = sy */
    /* short opcodes */
    op_rect         = 0x10,  /* low nibble is repeat count, followed by r,g,b */
    op_add_x_short  = 0x30,  /* low nibble is positive x-offset */
    op_add_y_short  = 0x60,  /* low nibble is positive y-offset */
    op_set_w_short  = 0x90,  /* low nibble is small w */
    op_set_h_short  = 0xB0,  /* low nibble is small h */
    op_run_c0_short = 0xD0,  /* low nibble is run length */
    op_run_c1_short = 0xF0,  /* low nibble is run length */
    /* long opcodes */
    op_add_x        = 0x20,  /* argument is formed by using */
    op_add_y        = 0x50,  /* ... low nibble and the next */
    op_sub_x        = 0x40,  /* ... two bytes. */
    op_sub_y        = 0x70,  /* ... */
    op_set_w        = 0x80,  /* ... */
    op_set_h        = 0xA0,  /* ... */
    op_run_c0       = 0xC0,  /* ... */
    op_run_c1       = 0xE0   /* ... */
};

/* Internal: opcode transformation */
#define op_short(op)  (op + 0x10)  /* op_add_N, op_set_N, op_run_C */
#define op_minus(op)  (op + 0x20)  /* op_add_N */

/* --- Encoding subroutines */

/* Internal: append a new block to chrome data */
private chrome_block *
chrome_newblock(chrome *cr)
{
    chrome_block *cb;
    if (! (cb = p2mem_alloc(cr->mem, sizeof(chrome_block))))
        return 0;
    cb->next = 0;
    cr->last ? (cr->last->next = cb) : (cr->first = cb);
    cr->last = cb;
    cr->wpos = cb->data;
    cr->wend = cb->data + chrome_block_len;
    cr->wlastrect = 0;
    return cb;
}

/* Internal: append multiple bytes to chrome data */
private int
chrome_write(chrome *cr, const byte *p, int len)
{
    int nlen;
    while (len > 0) {
        if (cr->wpos >= cr->wend)
            if (! chrome_newblock(cr))
                return_VMerror;
        nlen = min(len, cr->wend - cr->wpos);
        memcpy(cr->wpos, p, nlen);
        cr->wpos += nlen;
        p += nlen;
        len -= nlen;
    }
    return 0;
}

/* Internal: append a byte to chrome data */
#define chrome_putc(cr, b, onerror) { \
  if ((cr)->wpos >= (cr)->wend && !chrome_newblock(cr)) { onerror;} \
  *(cr)->wpos++ = (byte)(b); }

/* Internal: append short opcode */
#define chrome_put_short(cr, op, arg, onerror) { \
   ASSERT((arg)>=1 && (arg)<=0x10); \
   chrome_putc((cr), (byte)(op) + (byte)(arg) - 1, onerror); }

/* Internal: append long opcode */
#define chrome_put_long(cr, op, arg, onerror) { \
   ASSERT((arg)>=0 && (arg)<=0xfffff); \
   chrome_putc((cr), (byte)((arg)>>16) + (byte)(op), onerror); \
   chrome_putc((cr), (byte)((arg)>>8), onerror); \
   chrome_putc((cr), (byte)(arg), onerror); }

/* Internal: append long or short opcode */
#define chrome_put_cmd(cr, op, arg, onerror) { \
   if ((arg)<=0x10 && (arg)>=1) chrome_put_short(cr, op_short(op), arg, onerror) \
   else chrome_put_long(cr, op, arg, onerror) }

/* Internal: append rgb color */
#define chrome_put_rgb(cr, c, onerror) { \
   chrome_putc((cr), (byte)((c)>>16), onerror); \
   chrome_putc((cr), (byte)((c)>>8), onerror); \
   chrome_putc((cr), (byte)(c), onerror); }

/* Internal: put colors */
private int
chrome_put_c0c1(chrome *cr, gx_color_index c0, gx_color_index c1)
{
    if (c0 != gx_no_color_index && cr->c0 != c0) {
        chrome_putc(cr, op_set_c0, return_VMerror);
        chrome_put_rgb(cr, c0, return_VMerror);
        cr->c0 = c0;
    }
    if (c1 != gx_no_color_index && cr->c1 != c1) {
        chrome_putc(cr, op_set_c1, return_VMerror);
        chrome_put_rgb(cr, c1, return_VMerror); 
        cr->c1 = c1;
    }
    return 0;
}

/* Internal: put x y */
private int
chrome_put_xy(chrome *cr, int x, int y)
{
    /* Set X */
    if (x != cr->x) {
        if (x == cr->sx) {
            chrome_putc(cr, op_rest_x, return_VMerror);
        } else if (x > cr->x && x <= cr->x + 0x10) {
            chrome_put_short(cr, op_add_x_short, (x - cr->x), return_VMerror);
        } else if (x > cr->sx && x <= cr->sx + 0x10) {
            chrome_putc(cr, op_rest_x, return_VMerror);
            chrome_put_short(cr, op_add_x_short, (x - cr->sx), return_VMerror);
        } else if (x > cr->x) {
            chrome_put_long(cr, op_add_x, (x - cr->x), return_VMerror);
        } else {
            chrome_put_long(cr, op_sub_x, (cr->x - x), return_VMerror);
            cr->sx = x;
        }
        cr->x = x;
    }
    /* Set Y */
    if (y != cr->y) {
        if (y == cr->sy) {
            chrome_putc(cr, op_rest_y, return_VMerror);
        } else if (y > cr->y && y <= cr->y + 0x10) {
            chrome_put_short(cr, op_add_y_short, (y - cr->y), return_VMerror);
        } else if (y > cr->sy && y <= cr->sy + 0x10) {
            chrome_putc(cr, op_rest_y, return_VMerror);
            chrome_put_short(cr, op_add_y_short, (y - cr->sy), return_VMerror);
        } else if (y > cr->y) {
            chrome_put_long(cr, op_add_y, (y - cr->y), return_VMerror);
        } else {
            chrome_put_long(cr, op_sub_y, (cr->y - y), return_VMerror);
            cr->sy = y;
        }
        cr->y = y;
    }
    return 0;
}

/* Internal: put x y w h */
private int
chrome_put_xywh(chrome *cr, int x, int y, int w, int h)
{
    if (cr->w != w) {
        chrome_put_cmd(cr, op_set_w, w, return_VMerror);
        cr->w = w; 
    }
    if (cr->h != h) {
        chrome_put_cmd(cr, op_set_h, h, return_VMerror);
        cr->h = h;
    }
    return chrome_put_xy(cr, x, y);
}

/* Put end-of-object */
private int
chrome_put_end(chrome *cr)
{
    chrome_putc(cr, op_end, return_VMerror);
    chrome_reset_ctx(cr);
    return 0;
}

/* Appends a filled rectangle to the current image */
private int
chrome_put_rect(chrome *cr, int x, int y, int w, int h, gx_color_index c)
{
    int code;
    if (cr->c0 == c && cr->x == x && cr->wpos 
        && cr->wpos == cr->wlastrect + 4
        && cr->y == y && cr->w == w && cr->h == h
        && cr->wlastrect[0] >= op_rect 
        && cr->wlastrect[0] < op_rect + 0xf ) {
        /* Shortcut for adjacent rectangles */
        cr->wlastrect[0] += 1;
        cr->x += w;
        return 0;
    }
    /* Regular code */
    if ((code = chrome_put_xywh(cr, x, y, w, h)) < 0)
        return code;
    cr->wlastrect = cr->wpos;
    chrome_putc(cr, op_rect, return_VMerror);
    chrome_put_rgb(cr, c, return_VMerror);
    /* Setup state */
    cr->x = x + w;
    cr->c0 = c;
    return 0;
}

/* Appends a monochrome runmap to the current image */
private int
chrome_put_runmap(chrome *cr, runmap *r, gx_color_index c)
{
    int y;
    int code;
    /* Prepare */
    if (r == 0) 
        return 0;
    if ((code = chrome_put_c0c1(cr, gx_no_color_index, c)) < 0)
        return code; 
    /* Loop on y */
    for (y=r->ymin; y<=r->ymax; y++) {
        byte *data = r->data + r->rows[y - r->ymin];
        byte *dataend = r->data + r->rows[y - r->ymin + 1];
        if (data < dataend) {
            int x = 0;
            /* Set x, y */
            if ((code = chrome_put_xy(cr, 0, y)) < 0)
                return code;
            /* Store runs */
            while (data < dataend) {
                int x1 = x;
                rl_decode(data, x);
                if (data < dataend) {
                    int x2 = x;
                    rl_decode(data, x);
                    if (x2 > x1) 
                        chrome_put_cmd(cr, op_add_x, x2-x1, return_VMerror);
                    if (x > x2) 
                        chrome_put_cmd(cr, op_run_c1, x-x2, return_VMerror);
                }
            }
            cr->x = x;
        }
    }
    return 0;
}

/* Appends a monochrome bitmap to the current image.  The bitmap is described
   by <base>, <sourcex>, <sraster>, <xx>, <yy>, <w>, <h> as in ghostscript.
   Its colors are specified by arguments <c0> and <c1>. Colors can be
   <gx_no_color_index> to indicate transparency. */
private int
chrome_put_bitmap(chrome *cr, const byte * base, int sourcex, int sraster,
                  int xx, int yy, int w, int h, 
                  gx_color_index c0, gx_color_index c1)
{
    int y;
    int code;
    int basemask = 0x80 >> (sourcex & 0x7);
    base = base + (sourcex >> 3);
    /* Prepare */
    if ((code = chrome_put_c0c1(cr, c0, c1)) < 0)
        return code;
    /* Loop */
    for (y=yy; y<yy+h; y++, base+=sraster) {
        const byte *row = base;
        int z = (*row++);
        int m = basemask;
        int x1 = 0;
        int x2 = 0;
        int x = 0;
        /* Set x, y */
        if ((code = chrome_put_xy(cr, xx, y)) < 0)
            return code;
        /* Search runs */
        while (x < w) {
            x1 = x;
            while (x < w) {
                if (m) {
                    if (z & m)
                        break;
                    x += 1;
                    m >>= 1;
                } else {
                    z = (*row++);
                    while (z==0x00 && x+7<w) {
                        z = (*row++);
                        x += 8;
                    }
                    m = 0x80;
                }
            }
            x2 = x;
            while (x < w) {
                if (m) {
                    if (! (z & m))
                        break;
                    x += 1;
                    m >>= 1;
                } else {
                    z = (*row++);
                    while (z==0xff && x+7<w) {
                        z = (*row++);
                        x += 8;
                    }
                    m = 0x80;
                }
            }
            /* Write runs */
            if (x2 > x1) {
                if (c0 == gx_no_color_index) {
                    chrome_put_cmd(cr, op_add_x, (x2-x1), return_VMerror);
                } else {
                    chrome_put_cmd(cr, op_run_c0, (x2-x1), return_VMerror);
                }
            }
            if (x > x2) {
                if (c1 == gx_no_color_index) {
                    chrome_put_cmd(cr, op_add_x, (x-x2), return_VMerror);
                } else {
                    chrome_put_cmd(cr, op_run_c1, (x-x2), return_VMerror);
                }
            }
        }
        cr->x = xx + x;
    }
    return 0;
}

/* Appends a 24 bit color pixmap to the current image.  The pixmap is
   described by <base>, <sourcex>, <sraster>, <xx>, <yy>, <w>, <h> as in
   ghostscript */
private int
chrome_put_pixmap(chrome *cr, const byte * base, int sourcex, int sraster,
                  int xx, int yy, int w, int h)
{
    int y;
    int code;
    base = base + 3*sourcex;
    /* Prepare */
    if ((code = chrome_put_xywh(cr, xx, yy, w, h)) < 0) 
        return code;
    /* Store opcode */
    chrome_putc(cr, op_pixmap, return_VMerror);
    for (y=0; y<h; y++, base+=sraster) {
        if ((code = chrome_write(cr, base, 3*w)) < 0)
            return code;
    }
    /* Set state after pixmap opcode. */
    cr->x = xx + w;
    cr->sx = xx;
    cr->y = yy;
    return 0;
}

/* --- Decoding subroutines */

/* Internal: Data structure for reading data from chrome */
typedef struct chrome_reader_s {
    /* Reading point */
    chrome_block *block;
    byte *pos;
    byte *end;
} chrome_reader;

/* Internal: Read one byte */
#define chrome_getc(crd) \
    ((crd)->pos < (crd)->end ? *(crd)->pos++ : chrome_getc_sub(crd))

/* Internal: Read one byte, helper subroutine */
private inline byte 
chrome_getc_sub(chrome_reader *crd)
{
    if (!crd->block->next) return 0;
    crd->block = crd->block->next;
    crd->pos = crd->block->data;
    crd->end = crd->pos + chrome_block_len;
    return * crd->pos ++;
}

/* Internal: Read long argument */
#define chrome_getarg(crd, op, argvar) { \
    argvar = ((op & 0xf) << 16) | (chrome_getc(crd) << 8); \
    argvar |= chrome_getc(crd); }

/* Internal: Read bytes */
private inline void
chrome_read(byte *dest, chrome_reader *crd, int len)
{
    while (len > 0) {
        int nlen;
        if (crd->pos >= crd->end) {
            if (!crd->block->next) return;
            crd->block = crd->block->next;
            crd->pos = crd->block->data;
            crd->end = crd->pos + chrome_block_len;
        }
        nlen = min(len, crd->end - crd->pos);
        if (dest) {
            memcpy(dest, crd->pos, nlen);
            dest += nlen;
        }
        crd->pos += nlen;
        len -= nlen;
    }
}

/* Replay chrome data into a raster pixmap.
   See runmap_play() for an explanation of arguments 
   <base>, <sraster>, <cx>, <cy>, <cw>, and <ch>. */
private void
chrome_play(chrome_pos *cpos, 
            byte *base, int sraster, 
            int cx, int cy, int cw, int ch)
{
    chrome_reader crd;
    int cx2 = cx + cw;
    int cy2 = cy + ch;
    /* Opcode context */
    int x, sx, y, sy, w, h;
    byte c0[3], c1[3];
    byte *row;
    byte *c;
    int arg;
    /* Prepare chrome reader */
    crd.block = cpos->block;
    crd.pos = cpos->pos;
    crd.end = crd.block->data + chrome_block_len;
    /* Initialize opcode context */
    base = base - cy * sraster - cx * 3;
    c0[0] = c0[1] = c0[2] = 0;
    c1[0] = c1[1] = c1[2] = 0;
    x = sx = y = sy = 0;
    w = h = 1;
    row = 0;
    if (y>=cy && y<cy+ch) 
        row = base;
    /* Execute opcodes */
    for(;;) {
        byte op = chrome_getc(&crd);
        switch(op & 0xf0) {
            /* ---- Change X */
        case op_add_x:
            chrome_getarg(&crd, op, arg);
            x += arg;
            break;
        case op_sub_x:
            chrome_getarg(&crd, op, arg);
            x -= arg;
            sx = x;
            break;
        case op_add_x_short:
            x += 1 + (op & 0xf);
            break;
            /* ---- Change Y */
        case op_add_y:
            chrome_getarg(&crd, op, arg);
            y += arg;
            row = ((y>=cy && y <cy2) ? base + y * sraster : 0);
            break;
        case op_sub_y:
            chrome_getarg(&crd, op, arg);
            y -= arg;
            sy = y;
            row = ((y>=cy && y <cy2) ? base + y * sraster : 0);
            break;
        case op_add_y_short:
            y += 1 + (op & 0xf);
            row = ((y>=cy && y <cy2) ? base + y * sraster : 0);
            break;
            /* ---- Change W,H */
        case op_set_w:
            chrome_getarg(&crd, op, w);
            break;
        case op_set_w_short:
            w = 1 + (op & 0xf);
            break;
        case op_set_h:
            chrome_getarg(&crd, op, h);
            break;
        case op_set_h_short:
            h = 1 + (op & 0xf);
            break;
            /* ----- Runs */
        case op_run_c0:
            chrome_getarg(&crd, op, arg);
            c = c0;
            goto case_run;
        case op_run_c0_short:
            arg = 1 + (op & 0xf);
            c = c0;
            goto case_run;
        case op_run_c1:
            chrome_getarg(&crd, op, arg);
            c = c1;
            goto case_run;
        case op_run_c1_short:
            arg = 1 + (op & 0xf);
            c = c1;
        case_run:
            /* Draw run of length arg with color c */
            if (row) {
                byte *p1 = row + 3 * max(x,cx);
                byte *p2 = row + 3 * min(x+arg,cx2);
                while (p1 < p2) {
                    *p1++ = c[0];
                    *p1++ = c[1];
                    *p1++ = c[2];
                }
            }
            x += arg;
            break;
            /* ---- FillRect */
        case op_rect:
            arg = w * (1 + (op & 0xf));
            c0[0] = chrome_getc(&crd);
            c0[1] = chrome_getc(&crd);
            c0[2] = chrome_getc(&crd);
            /* Fill rectangle color c0 */
            if (y<cy2 && y+h>cy)
            {
                int y1 = max(y,cy);
                int y2 = min(y+h,cy2);
                int s1 = 3 * max(x,cx);
                int s2 = 3 * min(x+arg,cx2);
                byte *r = (row ? row : base + y1*sraster);
                while (y1 < y2) {
                    byte *p1 = r + s1;
                    byte *p2 = r + s2;
                    while (p1 < p2) {
                        *p1++ = c0[0];
                        *p1++ = c0[1];
                        *p1++ = c0[2];
                    }
                    y1 += 1;
                    r += sraster;
                }
            }
            x += arg;
            break;
        default:
            /* ---- Minor opcodes */
            switch(op) {
                /* ---- Pixmaps */
            case op_pixmap: 
                /* Fill pixmap with inline data */
                if (w>0) {
                    int ny = y;
                    int n1, n2, n3;
                    byte *r = 0;
                    while (ny < y+h && ny < cy) {
                        chrome_read(0, &crd, 3*w);
                        ny += 1;
                    }
                    n1 = (cx>x ? 3*(cx-x) : 0);
                    n3 = (x+w>cx2 ? 3*(x+w-cx2) : 0);
                    n2 = (3 * w) - n1 - n3;
                    r = base + ny * sraster + 3*max(x,cx);
                    while (ny < y+h && ny < cy2) {
                        chrome_read(0, &crd, n1);
                        chrome_read(r, &crd, n2);
                        chrome_read(0, &crd, n3);
                        r += sraster;
                        ny += 1;
                    }
                    while (ny < y+h) {
                        chrome_read(0, &crd, 3*w);
                        ny += 1;
                    }
                }
                sx = x;
                x += w;
                break;
                /* ---- Miscellaneous */
            case op_set_c0:
                c0[0] = chrome_getc(&crd);
                c0[1] = chrome_getc(&crd);
                c0[2] = chrome_getc(&crd);
                break;
            case op_set_c1:
                c1[0] = chrome_getc(&crd);
                c1[1] = chrome_getc(&crd);
                c1[2] = chrome_getc(&crd);
                break;
            case op_rest_x:
                x = sx;
                break;
            case op_rest_y:
                y = sy;
                row = ((y>=cy && y <cy2) ? base + y * sraster : 0);
                break;
            case op_end:
                return;
            default:
                ASSERT(0);
                break;
            }
            break;
        }
    }
}

/* Replay chrome data intersecting area <cw*ch+cx+cy>,
   (1) accumulating a colormap into <colorhash>,
   (2) storing colordata pointers into a memory buffer <base>.
   Aborts if the total number of colors exceeds <maxcolors>. */
private int
chrome_play_indexed(chrome_pos *cpos, 
                    htable *colorhash, int maxcolors,
                    colordata **base, int sraster, 
                    int cx, int cy, int cw, int ch)
{
    chrome_reader crd;
    int cx2 = cx + cw;
    int cy2 = cy + ch;
    /* Opcode context */
    int x, sx, y, sy, w, h;
    gx_color_index c0, c1;
    colordata *cd;
    colordata **row;
    byte rgb[3];
    int arg;
    /* Prepare chrome reader */
    crd.block = cpos->block;
    crd.pos = cpos->pos;
    crd.end = crd.block->data + chrome_block_len;
    /* Initialize opcode context */
    c0 = c1 = 0x000000;
    base = base - cy * sraster - cx;
    x = sx = y = sy = 0;
    w = h = 1;
    row = 0;
    if (y>=cy && y<cy+ch) 
        row = base;
    /* Execute opcodes */
    for(;;) {
        byte op = chrome_getc(&crd);
        switch(op & 0xf0) {
            /* ---- Change X */
        case op_add_x:
            chrome_getarg(&crd, op, arg);
            x += arg;
            break;
        case op_sub_x:
            chrome_getarg(&crd, op, arg);
            x -= arg;
            sx = x;
            break;
        case op_add_x_short:
            x += 1 + (op & 0xf);
            break;
            /* ---- Change Y */
        case op_add_y:
            chrome_getarg(&crd, op, arg);
            y += arg;
            row = ((y>=cy && y <cy2) ? base + y * sraster : 0);
            break;
        case op_sub_y:
            chrome_getarg(&crd, op, arg);
            y -= arg;
            sy = y;
            row = ((y>=cy && y <cy2) ? base + y * sraster : 0);
            break;
        case op_add_y_short:
            y += 1 + (op & 0xf);
            row = ((y>=cy && y <cy2) ? base + y * sraster : 0);
            break;
            /* ---- Change W,H */
        case op_set_w:
            chrome_getarg(&crd, op, w);
            break;
        case op_set_w_short:
            w = 1 + (op & 0xf);
            break;
        case op_set_h:
            chrome_getarg(&crd, op, h);
            break;
        case op_set_h_short:
            h = 1 + (op & 0xf);
            break;
            /* ----- Runs */
        case op_run_c0:
            chrome_getarg(&crd, op, arg);
            goto case_run0;
        case op_run_c0_short:
            arg = 1 + (op & 0xf);
        case_run0:
            if (row) {
                colordata **p1 = row + max(x,cx);
                colordata **p2 = row + min(x+arg,cx2);
                if (p2 > p1) {
                    if (! (cd = htable_lookup(colorhash, &c0, hash_rgb(c0), true)))
                        return_VMerror;
                    if (colorhash->nelems > maxcolors)
                        return 1;
                    while (p1 < p2)
                        *p1++ = cd;
                }
            }
            x += arg;
            break;
        case op_run_c1:
            chrome_getarg(&crd, op, arg);
            goto case_run1;
        case op_run_c1_short:
            arg = 1 + (op & 0xf);
        case_run1:
            if (row) {
                colordata **p1 = row + max(x,cx);
                colordata **p2 = row + min(x+arg,cx2);
                if (p2 > p1) {
                    if (! (cd = htable_lookup(colorhash, &c1, hash_rgb(c1), true)))
                        return_VMerror;
                    if (colorhash->nelems > maxcolors)
                        return 1;
                    while (p1 < p2)
                        *p1++ = cd;
                }
            }
            x += arg;
            break;
            /* ---- FillRect */
        case op_rect:
            arg = w * (1 + (op & 0xf));
            rgb[0] = chrome_getc(&crd);
            rgb[1] = chrome_getc(&crd);
            rgb[2] = chrome_getc(&crd);
            c0 = (rgb[0]<<16)|(rgb[1]<<8)|(rgb[2]);
            /* Fill rectangle */
            if (y<cy2 && y+h>cy)
            {
                int y1 = max(y,cy);
                int y2 = min(y+h,cy2);
                int s1 = max(x,cx);
                int s2 = min(x+arg,cx2);
                colordata **r = (row ? row : base + y1*sraster);
                if (! (cd = htable_lookup(colorhash, &c0, hash_rgb(c0), true)))
                    return_VMerror;
                if (colorhash->nelems > maxcolors)
                    return 1;
                while (y1 < y2) {
                    colordata **p1 = r + s1;
                    colordata **p2 = r + s2;
                    while (p1 < p2)
                        *p1++ = cd;
                    y1 += 1;
                    r += sraster;
                }
            }
            x += arg;
            break;
        default:
            /* ---- Minor opcodes */
            switch(op) {
                /* ---- Pixmaps */
            case op_pixmap:
                for (arg=y; arg<y+h; arg++) {
                    colordata **r = 0;
                    if (arg >= cy && arg < cy2)
                        r = base + arg * sraster;
                    for(sx=x; sx<x+w; sx++) {
                        gx_color_index c;
                        rgb[0] = chrome_getc(&crd);
                        rgb[1] = chrome_getc(&crd);
                        rgb[2] = chrome_getc(&crd);
                        c = (rgb[0]<<16)|(rgb[1]<<8)|(rgb[2]);
                        cd = htable_lookup(colorhash, &c, hash_rgb(c), true);
                        if (!cd) return_VMerror;
                        if (colorhash->nelems > maxcolors) 
                            return 1;
                        if (r && sx>=cx && sx<cx2)
                            r[sx] = cd;
                    }
                }
                sx = x;
                x += w;
                break;
                /* ---- Misc */
            case op_set_c0:
                rgb[0] = chrome_getc(&crd);
                rgb[1] = chrome_getc(&crd);
                rgb[2] = chrome_getc(&crd);
                c0 = (rgb[0]<<16)|(rgb[1]<<8)|(rgb[2]);
                break;
            case op_set_c1:
                rgb[0] = chrome_getc(&crd);
                rgb[1] = chrome_getc(&crd);
                rgb[2] = chrome_getc(&crd);
                c1 = (rgb[0]<<16)|(rgb[1]<<8)|(rgb[2]);
                break;
            case op_rest_x:
                x = sx;
                break;
            case op_rest_y:
                y = sy;
                row = ((y>=cy && y<cy2) ? base + y * sraster : 0);
                break;
            case op_end:
                return 0;
            default:
                ASSERT(0);
                break;
            }
            break;
        }
    }
}



/* ======================================
       C O L O R   R L E   F I L E S
   ====================================== */

/* Color RLE images start with the magic text "R6" followed by three integers
   representing the image width, the image height, and the palette size.
   These integers may be separated by blanks or by comments.  Comments start
   with character '#' and extend until the end of the line.  The last integer
   is followed by exactly one blank character (typically a new line).  Then
   come <palettesize> groups of three bytes (r,g,b) representing palette
   entries.  This is followed by a sequence of 32 bit integers (MSB first)
   whose low 20 bits represent a run length and high 12 bits represent the run
   color index.  Color indices are in range <0...palettesize-1>.  Color
   indices greater than 0xff0 are reserved.  Indice 0xfff represents the
   transparent color.  Runs are arranged in lexicographical order (top to
   bottom).  A line ends when the sum of run lengths is a multiple of the
   image width. */


/* Internal: Transparent color */
#define crle_transparent (0xfff)

/* Internal: Size of buffer for run output */
#define crle_buffer_size (512*4)

/* Internal: Data structure to keep track of run output */
typedef struct crle_output_s {
    FILE *f;
    int x, y;
    int w, h;
    byte *bufpos;
    byte  buffer[crle_buffer_size];
} crle_output;

/* Internal: Initialize crle_output data structure */
private int
crle_open(crle_output *crle, FILE *f, int w, int h, int ncolors, 
          gx_color_index *palette, const char *comment)
{
    int i;
    crle->f = f;
    crle->x = crle->y = 0;
    crle->w = w;
    crle->h = h;
    crle->bufpos = crle->buffer;
    fprintf(f, "R6\n");
    if (comment) fprintf(f, "# %s\n", comment);
    fprintf(f, "%d %d\n%d\n", w, h, ncolors);
    for(i=0; i<ncolors; i++) {
        byte rgb[3];
        rgb[0] = palette[i]>>16;
        rgb[1] = palette[i]>> 8;
        rgb[2] = palette[i];
        if (fwrite(rgb, 1, 3, f) != 3)
            return_error(gs_error_ioerror);            
    }
    return 0;
}

/* Internal: Outputs data buffered in crle_output data structure */
private int
crle_flush(crle_output *crle)
{
    int len = crle->bufpos - crle->buffer;
    if (len > 0) {
        if (fwrite(crle->buffer, 1, len, crle->f) != len)
            return_error(gs_error_ioerror);
        crle->bufpos = crle->buffer;
    }
    return 0;
}

/* Internal: Appends a run to crle output */
private inline int 
crle_write(crle_output *crle, int index, int len) 
{
    int code;
    int datum;
    /* Handle overflows */
    while (len > 0xfffff) {
        len -= 0xfffff;
        if ((code = crle_write(crle, index, 0xfffff)) < 0)
            return code;
    }
    /* Append run description */
    crle->x += len;
    datum = (index<<20) | len;
    *crle->bufpos++ = datum >> 24;
    *crle->bufpos++ = datum >> 16;
    *crle->bufpos++ = datum >> 8;
    *crle->bufpos++ = datum >> 0;
    if (crle->bufpos >= crle->buffer + crle_buffer_size - 4)
        if ((code = crle_flush(crle)) < 0)
            return code;
    /* Check line change */
    if (crle->x >= crle->w) {
        ASSERT(crle->x == crle->w);
        ASSERT(crle->y < crle->h); 
        crle->y += 1;
        crle->x = 0;
    }
    return 0;
}

/* Internal: Sorted list of foreground components */
typedef struct crle_comp_s {
    runmap *rmap;
    gx_color_index color;
    struct crle_comp_s *next;
} crle_comp;

/* Internal: Y-sorting subroutine for foreground components */
private int
crle_comp_sub(const void *pa, const void *pb)
{
    const crle_comp *a = (const crle_comp*)pa;
    const crle_comp *b = (const crle_comp*)pb;
    return a->rmap->ymin - b->rmap->ymin;
}

/* Internal: Description of a run */
typedef struct crle_run_s {
    uint c; /* palette color index */
    uint x; /* start of run */
    uint l; /* length of run */
} crle_run;

/* Internal: X-sorting subroutine for runs */
private int
crle_run_sub(const void *pa, const void *pb)
{
    const crle_run *a = (const crle_run*)pa;
    const crle_run *b = (const crle_run*)pb;
    return a->x - b->x;
}

/* Save a CRLE file */
private int 
crle_save(p2mem *mem, 
          int width, int height,
          gx_color_index *palette, uint palettesize,
          int num, runmap **rmaps, gx_color_index *colors,
          const char *comment, 
          FILE *outfile)
{
    int i, x, y;
    int code;
    crle_output out;
    runmap *r;
    gx_color_index c;
    crle_comp *fgarray = 0;
    crle_comp *fglist;
    crle_comp **fgp;
    int maxruns = 0;
    int nruns = 0;
    crle_run *runs = 0;
    /* Allocate and fill fgarray */
    fgarray = p2mem_alloc(mem, sizeof(crle_comp) * num);
    if (num && !fgarray) return_VMerror;
    for (i=0; i<num; i++) {
        fgarray[i].rmap = rmaps[i];
        fgarray[i].color = colors[i];
    }
    /* Sort and link fgarray */
    if (num > 1)
        qsort(fgarray, num, sizeof(crle_comp), crle_comp_sub);
    fglist = fgarray;
    for (i=1; i<num; i++)
        fgarray[i-1].next = &fgarray[i];
    if (num>0)
        fgarray[num-1].next = 0;
    /* Prepare output */
    code = crle_open(&out, outfile, width, height, 
                     palettesize, palette, comment);
    if (code < 0) goto xit;
    /* Iterate on rows */
    for (y=0; y<height; y++) {
        nruns = 0;
        fgp = &fglist;
        while ( (*fgp) && (*fgp)->rmap->ymin <= y) {
            c = (*fgp)->color;
            r = (*fgp)->rmap;
            if (r->ymax < y) {
                /* This drawlist has been completely processed */
                (*fgp) = (*fgp)->next;
            } else {
                /* This drawlist intersects the current row */
                byte *data = r->data + r->rows[y - r->ymin];
                byte *dataend = r->data + r->rows[y - r->ymin + 1];
                int reqruns = nruns + (dataend - data);
                if (reqruns > maxruns) {
                    int nmaxruns = max(reqruns, maxruns+maxruns);
                    runs = p2mem_resize(mem, runs, sizeof(crle_run)*nmaxruns);
                    if (!runs) goto xitmem;
                    maxruns = nmaxruns;
                }
                x = 0;
                while (data < dataend) {
                    rl_decode(data, x);
                    runs[nruns].x = x;
                    if (data < dataend) {
                        rl_decode(data, x);
                        runs[nruns].l = x - runs[nruns].x;
                        runs[nruns].c = c;
                        nruns += 1;
                    }
                }
                fgp = &(*fgp)->next;
            }
        }
        /* Sort all runs on the row */
        if (nruns > 1)
            qsort(runs, nruns, sizeof(crle_run), crle_run_sub);
        /* Output row */
        x = 0;
        for (i=0; i<nruns; i++) {
            ASSERT(runs[i].x >= x);
            if (runs[i].x > x) {
                crle_write(&out, crle_transparent, runs[i].x - x);
                x = runs[i].x;
            }
            if (runs[i].l > 0) {
                crle_write(&out, runs[i].c, runs[i].l);
                x += runs[i].l;
            }
        }
        ASSERT(x <= width);
        if (x < width)
            crle_write(&out, crle_transparent, width - x);
    }
    /* Success */
    code = crle_flush(&out);
    goto xit;
 xitmem:
    /* Memory error */
    code = gs_note_error(gs_error_VMerror);
 xit:
    /* Cleanup and terminate */
    p2mem_free(mem, runs);
    p2mem_free(mem, fgarray);
    return code;
}



/* ======================================
       D R A W I N G     L I S T
   ====================================== */

/* The first job of the gdevdjvu devices is to build a list of drawlist records
   representing all drawing operations.  Each drawlist record essentially
   contains a flag word (see below), a mask runmap describing its extent, 
   and a color index representing the single color of the drawing element
   (or gx_no_color_index when the drawing element contains several colors). */

/* Drawlist head */
typedef struct drawlist_head_s {
    struct drawlist_s *first;
    struct drawlist_s *last;
} drawlist_head;

/* Drawlist element */
typedef struct drawlist_s {
    struct drawlist_s *next;   /* linked list pointer */
    struct drawlist_s *prev;   /* linked list pointer */
    runmap *mask;              /* shape of the object */
    uint flags;                /* type for the object (see below) */
    gx_color_index color;      /* pure color or gx_no_color_index */
    chrome_pos *chromepos;     /* non zero when color is gx_no_color_index */
    struct txtmark_s *text;    /* textual information when available */
    struct drawlist_s *bglink; /* cf. section drawlist classification */
    int perimeter;             /* cf. section drawlist classification */
    int area;                  /* cf. section drawlist classufication */
    int palcolor;              /* cf. section foreground colors */
} drawlist;

/* Drawlist flags. */
#define DLIST_FOREGROUND     (0x1)     /* Assigned to foreground */
#define DLIST_BACKGROUND     (0x2)     /* Assigned to background */
#define DLIST_INPURE         (0x10)    /* Has more than one color */
#define DLIST_TEXT           (0x20)    /* Drawing arises from text operators */
#define DLIST_PATH           (0x40)    /* Drawing arises from path operators */
#define DLIST_PATH_FILL      (0x80)    /* Drawing arises from path fill operator */
#define DLIST_HASGRAY        (0x1000)  /* Color(s) may include non b&w colors */
#define DLIST_HASCOLOR       (0x2000)  /* Color(s) may include non gray colors */
#define DLIST_RECURSIVE      (0x10000) /* Already processing high level routine */

/* Return flags indicating whether color represents
   black or white , gray shades, or actual colors. */
#define color_flags(color) \
  ((color^(color>>8))&0xffff) ? (DLIST_HASCOLOR|DLIST_HASGRAY) \
    : (color && color!=0xffffff) ? (DLIST_HASGRAY) : 0 

/* Appends a new drawlist element. */
private drawlist *
drawlist_append(p2mem *mem, drawlist_head *head)
{
    drawlist *d = p2mem_alloc(mem, sizeof(drawlist));
    if (!d) return 0;
    memset(d, 0, sizeof(drawlist));
    d->color = gx_no_color_index;
    d->next = 0;
    d->prev = head->last;
    head->last = d;
    (d->prev) ? (d->prev->next = d) : (head->first = d);
    return d;
}

/* Inserts a new drawlist element before a given element. */
private drawlist *
drawlist_insert(p2mem *mem, drawlist_head *head, drawlist *pos)
{
    drawlist *d = p2mem_alloc(mem, sizeof(drawlist));
    if (!d) return 0;
    memset(d, 0, sizeof(drawlist));
    d->color = gx_no_color_index;
    d->next = pos;
    d->prev = (pos ? pos->prev : head->last);
    (d->next) ? (d->next->prev = d) : (head->last  = d);
    (d->prev) ? (d->prev->next = d) : (head->first = d);
    return d;
}

/* Destroy a drawlist element. */
private void
drawlist_remove(p2mem *mem, drawlist_head *head, drawlist *d)
{
    (d->next) ? (d->next->prev = d->prev) : (head->last = d->prev);
    (d->prev) ? (d->prev->next = d->next) : (head->first = d->next);
    if (mem) {
        if (d->mask) 
            runmap_free(d->mask);
        p2mem_free(mem, d);
    }
}

/* Compact display for debugging purposes. */
#ifdef DEBUG
private void
drawlist_print(drawlist *dl, gx_device_djvu *cdev)
{
    fprintf(STDOUT,"  --");
    if (dl->flags & DLIST_FOREGROUND)
	fprintf(STDOUT," fg");
    else if (dl->flags & DLIST_BACKGROUND)
	fprintf(STDOUT," bg");
    if (dl->mask)
        fprintf(STDOUT," %dx%d+%d+%d", 
                dl->mask->xmax - dl->mask->xmin +1, 
                dl->mask->ymax - dl->mask->ymin +1, 
                dl->mask->xmin, 
                dl->mask->ymin );
    if (dl->color != gx_no_color_index)
        fprintf(STDOUT," purecolor=%06x", (uint)(dl->color));
    if (dl->flags & DLIST_HASCOLOR)
        fprintf(STDOUT," hascolor");
    else if (dl->flags & DLIST_HASGRAY)
        fprintf(STDOUT," hasgray");
    if (dl->flags & DLIST_INPURE)
        fprintf(STDOUT," inpure");
    if (dl->flags & DLIST_PATH)
        fprintf(STDOUT," path");
    if (dl->flags & DLIST_PATH_FILL)
        fprintf(STDOUT," fill");
    if (dl->flags & DLIST_TEXT)
        fprintf(STDOUT," text");
    if (dl->perimeter)
	fprintf(STDOUT," peri=%d", dl->perimeter);
    if (dl->area)
	fprintf(STDOUT," area=%d", dl->area);
    if (dl->palcolor)
	fprintf(STDOUT," palcolor=%d", dl->palcolor);
    if (dl->chromepos)
	fprintf(STDOUT," haschrome");
    fprintf(STDOUT,"\n");
}
#endif

/* Render drawlist entries that have specific flags set.
   See runmap_play() for an explanation of arguments 
   base, sraster, cx, cy, cw, and ch. */
private void
drawlist_play(drawlist_head *head, uint flags,
              byte *base, int sraster, int cx, int cy, int cw, int ch)
{
    /* Make up for negative coordinates in rectangle */
    if (cx < 0) {
        cw += cx;
        base -= cx;
        cx -= cx;
    }
    if (cy < 0) {
        ch += cy;
        base -= cy * sraster;
        cy -= cy;
    }
    /* Render */
    if (cw > 0 && ch > 0) {
        drawlist *dl;
        for (dl = head->first; dl; dl = dl->next) {
            if (!(dl->mask) || 
                !(dl->flags & flags) ||
                dl->mask->ymax < (uint)(cy) ||
                dl->mask->xmax < (uint)(cx) || 
                dl->mask->ymin >= (uint)(cy + ch) ||
                dl->mask->xmin >= (uint)(cx + cw)  )
                continue;
            /* Render drawlist element */
            if (dl->color != gx_no_color_index) {
                runmap_play(dl->mask, dl->color, base, sraster, cx, cy, cw, ch);
            } else {
                ASSERT(dl->chromepos);
                chrome_play(dl->chromepos, base, sraster, cx, cy, cw, ch);
            }
        }
    }
}




/* ======================================
       P D F M A R K S 
   ====================================== */

/* This code generates comments in the sep files
   on the basis of pdfmark sequences in the postscript.

   -- LNK pdfmark sequences are generated by the pdf interpreter. 
   More specifically:
        [ /Rect [ x1 y1 x2 y2 ] /URI (uri) /LNK pdfmark
     or [ /Rect [ x1 y1 x2 y2 ] /Page (pageno) /LNK pdfmark
     or [ /Rect [ x1 y1 x2 y2 ] /A << ... /URI (uri) .. >> /LNK pdfmark
     or [ /Rect [ x1 y1 x2 y2 ] ... /Action /Launch /Subtype /Link /ANN pdfmark
   produces the following comment in the sep file:
        # L <w>x<h>+<x>+<y> (uri)
   It also converts uri starting with page:// into uri starting with #.

   -- OUT pdfmarks sequences are generated by the pdf interpreter
   More specifically
        [ /Count <c> /Page <p> /Title (<ttl>) /OUT pdfmark
   produces the following comment in the sep file:
        # B <c> (<ttl>) (#<p>)

   -- PAGELABEL pdfmarks (pdf interpreter does not do this well)
        [ /Label (<ttl>) /PAGELABEL pdfmark
   produces the following comment in the sep file:
        # T (<ttl>)
*/


/* The txtmark structure.
   -- These are used for textmarks (see djvutext.ps)
   and are pointed to by the field "text" in the 
   drawlist structure. */
typedef struct txtmark_s txtmark;
struct txtmark_s {
    int x,y,w,h;
    char *arg;
};

/* The pdfmark structure.  
   -- These are stored in a linked list pointed to by 
      the field "marks" in the djvu device */
typedef struct pdfmark_s pdfmark;
struct pdfmark_s {
    pdfmark *next;
    int type;
    union {
        struct pdfmark_l_s {
            int x,y,w,h;
            char *uri;
        } l;
        struct pdfmark_b_s {
            int count;
            char *title;
            int page;
        } b;
        struct pdfmark_p_s {
	    char *title;
	} p;
    } u;
};

/* Helper: Convert unicode to utf8 */
private int
unicode_to_utf8(gs_char unicode, byte *utf)
{
    int i = 0;
    if (unicode <= 0) {
        utf[0] = 0;
    } else if (unicode <= 0x7f) {
        utf[i++] = unicode;
    } else if (unicode <= 0x7ff) {
        utf[i++] = 0xc0 + ((unicode >> 6) & 0x1f);
        utf[i++] = 0x80 + (unicode & 0x3f);
    } else if (unicode <= 0xffff) {
        utf[i++] = 0xe0 + ((unicode >> 12) & 0x0f);
        utf[i++] = 0x80 + ((unicode >> 6) & 0x3f);
        utf[i++] = 0x80 + (unicode & 0x3f);
    } else if (unicode <= 0x1fffff) {
        utf[i++] = 0xf0 + ((unicode >> 18) & 0x07);
        utf[i++] = 0x80 + ((unicode >> 12) & 0x3f);
        utf[i++] = 0x80 + ((unicode >> 6) & 0x3f);
        utf[i++] = 0x80 + (unicode & 0x3f);
    } else if (unicode <= 0x3ffffff) {
        utf[i++] = 0xf8 + ((unicode >> 24) & 0x03);
        utf[i++] = 0x80 + ((unicode >> 18) & 0x3f);
        utf[i++] = 0x80 + ((unicode >> 12) & 0x3f);
        utf[i++] = 0x80 + ((unicode >> 6) & 0x3f);
        utf[i++] = 0x80 + (unicode & 0x3f);
    } else {
        utf[i++] = 0xfc + ((unicode >> 30) & 0x01);
        utf[i++] = 0x80 + ((unicode >> 24) & 0x3f);
        utf[i++] = 0x80 + ((unicode >> 18) & 0x3f);
        utf[i++] = 0x80 + ((unicode >> 12) & 0x3f);
        utf[i++] = 0x80 + ((unicode >> 6) & 0x3f);
        utf[i++] = 0x80 + (unicode & 0x3f);
    }
    utf[i] = 0;
    return i;
}

/* Helper-- convert pdfdocencoding to utf8 */
private int
pdfdocencoding_to_utf8(byte b, byte *utf)
{
    static const gs_char at0x18[] = { 
        0x02d8, 0x02c7, 0x02c6, 0x02d9, 0x02dd, 0x02db, 0x02da, 0x02dc };
    static const gs_char at0x80[] = {
        0x2022, 0x2020, 0x2021, 0x2026, 0x2014, 0x2013, 0x0192, 0x2044, 
        0x2039, 0x203a, 0x2212, 0x2030, 0x201e, 0x201c, 0x201d, 0x2018,
        0x2019, 0x201a, 0x2122, 0xfb01, 0xfb02, 0x0141, 0x0152, 0x0160,
        0x0178, 0x017d, 0x0131, 0x0142, 0x0153, 0x0161, 0x017e, 0xfffd,
        0x20ac };
    gs_char unicode = (gs_char)b;
    if (b >= 0x18 && b < 0x20)
        unicode = at0x18[b-0x18];
    else if (b >= 0x80 && b <= 0xa0)
        unicode = at0x80[b-0x80];
    return unicode_to_utf8(unicode, utf);
}

/* Helper: return postscript string encoded on 7 bit ascii. */
private char *
make_ps_string(p2mem *mem, const byte *p, uint l)
{
    int i;
    int n = l+3;
    char *str, *s;
    static const char odigit[] = "01234567";
    for (i=0; i<l; i++)
        if (p[i]=='(' || p[i]==')' || p[i]=='\\')
            n += 1;
        else if (! (p[i]>=0x20 && p[i]<0x7f))
            n += 3;
    if (! (str = s = p2mem_alloc(mem, n)))
        return 0;
    *s++ = '(';
    for (i=0; i<l; i++) {
        byte b = p[i];
        if (b=='(' || b==')' || b=='\\') {
            *s++ = '\\';
            *s++ = b;
        } else if ( b>=0x20 && b<0x7f ) {
            *s++ = b;
        } else {
            *s++ = '\\';
            *s++ = odigit[(b>>6)&0x3];            
            *s++ = odigit[(b>>3)&0x7];
            *s++ = odigit[(b>>0)&0x7];
        }
    }
    *s++ = ')';
    *s++ = 0;
    return str;
}

/* Internal */
private int
pdfmark_eq(const gs_param_string *p, const char *q)
{
    return ( p->size == strlen(q) && !memcmp(p->data, q, p->size) );
}

/* Internal */
private int
pdfmark_dup(p2mem *mem, const gs_param_string *p, char **out)
{
    if (! (*out = p2mem_alloc(mem, p->size+1))) return_VMerror;
    memcpy(*out, p->data, p->size);
    (*out)[p->size] = 0;
    return 0;
}

/* Internal -- recode psstring as utf8 */
private int
pdfmark_recode(p2mem *mem, const gs_param_string *p, char **out)
{
    int len = p->size;
    const byte *data = p->data;
    byte *tmp = p2mem_alloc(mem, len+1);
    byte *ptr = tmp;
    byte *utf = 0;
    int ret = 0;
    int i;
    if (!tmp) goto vmerror;
    if (len < 2) goto error;
    if (data[0]=='<' && data[len-1]=='>' && !(len&1)) {
        for (i=1; i<len-1; i++) {
            byte b;
            if (data[i]>='0' && data[i]<='9') 
                b = data[i] - '0';
            else if (data[i]>='A' && data[i]<='F') 
                b = data[i] - 'A' + 10;
            else if (data[i]>='a' && data[i]<='f') 
                b = data[i] - 'a' + 10;
            else 
                goto error;
            if (i&1) 
                *ptr = (b<<4); 
            else 
                *ptr++ |= b;
        }
    } else if (data[0] == '(' && data[len-1] == ')') {
        for (i=1; i<len-1; i++) {
            byte b = data[i];
            if (b=='\\') {
                char *s;
                static const char *s1 = "nrtbf";
                static const char *s2 = "\n\r\t\b\f";
                if (++i >= len-1) goto error;
                b = data[i];
                if ((s = strchr(s1, b))) 
                    b = s2[s-s1];
                else if (b >= '0' || b <= '7') {
                    b = b - '0';
                    if (i+1 < len-1 && data[i+1] >= '0' && data[i+1] <= '7')
                        b = (b << 3) + (data[++i]-'0');
                    if (i+1 < len-1 && data[i+1] >= '0' && data[i+1] <= '7')
                        b = (b << 3) + (data[++i]-'0');
                } 
            }
            *ptr++ = b;
        }
    } else
        goto error;
    len = ptr - tmp;
    if (! (utf = p2mem_alloc(mem, 3*len))) goto vmerror;
    ptr = utf;
    if (len>=2 && tmp[0]==0xfe && tmp[1]==0xff) // unicode
        for (i=0; i<len-1; i+=2)
            ptr += unicode_to_utf8((gs_char)(tmp[i+1]+((int)tmp[i]<<8)), ptr);
    else
        for (i=0; i<len; i++)
            ptr += pdfdocencoding_to_utf8(tmp[i], ptr);
    len = ptr - utf;
    if ((*out = make_ps_string(mem, utf, len)))
        goto xit;
 error:
    ret = 1;
    goto xit;
 vmerror:
    ret = -1;
 xit:
    p2mem_free(mem, tmp); 
    p2mem_free(mem, utf); 
    if (ret != 0) { *out = 0; }
    if (ret < 0) return_VMerror;
    return ret;
}

/* Internal */
private int
pdfmark_find(p2mem *mem, gs_param_string_array *pma, 
             const char *key, char **out)
{
    int i;
    for (i=0; i<pma->size-2; i+=2)
        if (pdfmark_eq(&pma->data[i], key)) 
            return pdfmark_dup(mem, &pma->data[i+1], out);
    *out = 0;
    return 1;
}

/* Internal */
private int
pdfmark_find_recode(p2mem *mem, gs_param_string_array *pma, 
                    const char *key, char **out)
{
    int i;
    for (i=0; i<pma->size-2; i+=2)
        if (pdfmark_eq(&pma->data[i], key)) 
            return pdfmark_recode(mem, &pma->data[i+1], out);
    *out = 0;
    return 1;
}

/* Internal -- obtain current matrix */
private int
pdfmark_get_ctm(p2mem *mem, gs_param_string_array *pma, gs_matrix *pctm)
{
    char c;
    int code;
    char *temp = 0;
    if ((code = pdfmark_dup(mem, &pma->data[pma->size-2], &temp)) < 0)
        return code;
    code = gs_error_rangecheck;
    if (temp && sscanf(temp, "[ %g %g %g %g %g %g %c",
                       &pctm->xx, &pctm->xy, &pctm->yx, &pctm->yy, 
                       &pctm->tx, &pctm->ty, &c) == 7)
        code = 0;
    p2mem_free(mem, temp);
    return code;
}

/* Internal -- obtain rect argument */
private int
pdfmark_get_rect(p2mem *mem, gs_param_string_array *pma, 
                 const char *argname, gs_rect *prect)
{
    char *temp = 0;
    int code;
    char c;
    if ((code = pdfmark_find(mem, pma, argname, &temp)) < 0)
        return code;
    code = gs_error_rangecheck;
    if (temp && sscanf(temp, "[ %lg %lg %lg %lg %c",
                       &prect->p.x, &prect->p.y,
                       &prect->q.x, &prect->q.y, &c ) == 5)
        code = 0;
    p2mem_free(mem, temp);
    return code;
} 

/* Internal -- process LNK mark */
private int
pdfmark_do_lnk(p2mem *mem, gs_param_string_array *pma, pdfmark **pmark)
{
    gs_matrix ctm;
    gs_rect rect, urect;
    pdfmark *gmark = 0;
    struct pdfmark_l_s *mark = 0;
    char *temp0 = 0;
    char *temp1 = 0;
    char *temp2 = 0;
    int code = 0;
    /* Allocate mark */
    code = gs_error_VMerror; 
    if (! (gmark = p2mem_alloc(mem, sizeof(pdfmark)))) 
        goto xit;
    memset(gmark, 0, sizeof(pdfmark));
    mark = &gmark->u.l;
    /* Search Rect */
    if (((code = pdfmark_get_ctm(mem, pma, &ctm)) < 0) ||
        ((code = pdfmark_get_rect(mem, pma, "/Rect", &urect)) < 0) )
        goto xit;
    gs_bbox_transform(&urect, &ctm, &rect);
    mark->x = (int)floor(rect.p.x + 0.5);
    mark->y = (int)floor(rect.p.y + 0.5);
    mark->w = (int)floor(rect.q.x + 0.5) - mark->x;
    mark->h = (int)floor(rect.q.y + 0.5) - mark->y;
    /* Search URI */
    if (((code = pdfmark_find(mem, pma, "/A", &temp0)) < 0) ||
        ((code = pdfmark_find(mem, pma, "/URI", &temp1)) < 0) ||
        ((code = pdfmark_find(mem, pma, "/Page", &temp2)) < 0 ) )
        goto xit;
    if (temp0) {  /* Action: "<< ... /URI (...) >>" */
        int state = 0;
        const char *c, *s;
        for (c=s=temp0; *s; s++) {
            if (state<4 && *s == ("/URI")[state]) {
                state += 1;
            } else if (state==4 && *s=='(') {
                c = s;
                state = 5;
            } else if (state==4 && strchr(" \n\r\t", *s)) {
                state = 4;
            } else if (state==5 && *s=='\\') {
                state = 6;
            } else if (state==5 && *s==')') {
                code = gs_error_VMerror; 
                if (! (mark->uri = p2mem_alloc(mem, s-c+2)))
                    goto xit;
                memcpy(mark->uri, c, s-c+1);
                mark->uri[s-c+1] = 0;
                break;
            } else if (state>=5) {
                state = 5;
            } else {
                state = ((*s == '/') ? 1 : 0);
            }
        }
    } else if (temp1) { /* Format: "(...)" */
        mark->uri = temp1;
        temp1 = 0;
    } else if (temp2) { /* Format: "123" */
        code = gs_error_VMerror; 
        if (! (mark->uri = p2mem_alloc(mem, strlen(temp2) + 4)))
            goto xit;
        strcpy(mark->uri, "(#");
        strcat(mark->uri, temp2);
        strcat(mark->uri, ")");
    }
    /* Check if we were able to understand this LNK */
    code = 2;
    if (! mark->uri) goto xit;
    /* Postprocess page:// urls */
    if (! strncmp(mark->uri, "(page://", 8)) {
        memcpy(mark->uri, "(#", 2);
        memmove(mark->uri+2, mark->uri+8, 1+strlen(mark->uri+8));
    }
    /* The end */
    gmark->type = 'L';
    *pmark = gmark;
    gmark = 0;
    code = 0;
    /* PDF defines numerous other fields definining borders and hiliting.  Few
       can be used because Ghostscript only passes /View and /Border, and
       because DjVu borders and hiliting modes are different. */
 xit:
    p2mem_free(mem, gmark);
    p2mem_free(mem, temp0);
    p2mem_free(mem, temp1);
    p2mem_free(mem, temp2);
    return code;
}


/* Internal -- process OUT mark */
private int
pdfmark_do_out(p2mem *mem, gs_param_string_array *pma, pdfmark **pmark)
{
    int code;
    pdfmark *gmark = 0;
    struct pdfmark_b_s *mark = 0;
    char *temp0 = 0;
    char *temp1 = 0;
    char *temp2 = 0;
    /* Allocate mark */
    code = gs_error_VMerror; 
    if (! (gmark = p2mem_alloc(mem, sizeof(pdfmark)))) 
        goto xit;
    memset(gmark, 0, sizeof(pdfmark));
    mark = &gmark->u.b;
    /* Search count */
    if (((code = pdfmark_find(mem, pma, "/Count", &temp0)) < 0) ||
        ((code = pdfmark_find(mem, pma, "/Page", &temp1)) < 0) ||
        ((code = pdfmark_find_recode(mem, pma, "/Title", &temp2)) < 0) )
        goto xit;
    code = 2;
    if (temp0 && sscanf(temp0, "%d", &mark->count) != 1) goto xit;
    mark->count = abs(mark->count);
    if (! temp1 || sscanf(temp1, "%d", &mark->page) != 1) goto xit;
    if (mark->page < 0 || ! temp2) goto xit;
    mark->title = temp2;
    temp2 = 0;
    /* The end */
    gmark->type = 'B';
    *pmark = gmark;
    gmark = 0;
    code = 0;
 xit:
    p2mem_free(mem, gmark);
    p2mem_free(mem, temp0);
    p2mem_free(mem, temp1);
    p2mem_free(mem, temp2);
    return code;
}



/* ======================================
       D J V U    D E V I C E
   ====================================== */

/* Declaring the process procedure. */
#define djvu_proc_process(process) int process(gx_device_djvu *dev)

/* The device structure. */
struct gx_device_djvu_s {
    gx_device_common;
    /* Procedure vector */
    djvu_proc_process((*process));
    /* Parameters */
    int    threshold;
    int    fgcolors;
    int    fgimgcolors;
    int    bgsubsample;
    long   maxbitmap;
    char   outputfilename[prn_fname_sizeof];
    bool   reopenperpage;
    bool   autohires;
    bool   extracttext;
    bool   quiet;
    /* Device data */
    p2mem          *pmem;          /* stable memory */
    FILE           *outputfile;    /* output filename */
    drawlist_head   head;          /* drawlist */
    chrome         *gchrome;       /* chrome data */
    gx_color_index *fgpalette;     /* fgcolors */
    uint            fgpalettesize; /* number of fgcolors */
    pdfmark        *marks;         /* hyperlinks etc. */
    /* Current drawlist component */
    runmapbld      *curmask;       /* accumulated mask */
    gx_color_index  curcolor;      /* color or gx_no_color_index */
    uint            curflags;      /* flags */
    chrome_pos     *curchrome;     /* optional chrome position */
    bool            curbreak;      /* see comment in begin_typed_image */
};

gs_private_st_suffix_add0_final(st_device_djvu,gx_device_djvu,"gx_device_djvu", 
                                device_djvu_enum_ptrs, device_djvu_reloc_ptrs,
                                gx_device_finalize, st_device );

/* Device procedures. */
private dev_proc_open_device(djvu_open);
private dev_proc_output_page(djvu_output_page);
private dev_proc_close_device(djvu_close);
private dev_proc_get_params(djvu_get_params);
private dev_proc_put_params(djvu_put_params);
private dev_proc_fill_rectangle(djvu_fill_rectangle);
private dev_proc_copy_mono(djvu_copy_mono);
private dev_proc_copy_color(djvu_copy_color);
private dev_proc_stroke_path(djvu_stroke_path);
private dev_proc_fill_path(djvu_fill_path);
private dev_proc_begin_typed_image(djvu_begin_typed_image);
private dev_proc_text_begin(djvu_text_begin);

/* Device procs */
#define djvu_device_procs { \
    djvu_open, 0, 0, djvu_output_page, djvu_close, \
    gx_default_rgb_map_rgb_color, gx_default_rgb_map_color_rgb, \
    djvu_fill_rectangle, 0, djvu_copy_mono, djvu_copy_color, \
    0, 0, djvu_get_params, djvu_put_params, 0, 0, 0, 0, \
    gx_page_device_get_page_device, 0, 0, 0, 0, \
    djvu_fill_path, djvu_stroke_path, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    djvu_begin_typed_image, 0, 0, 0, 0, djvu_text_begin }

/* Setting up a mask device */
#define djvu_device_body(name,procs,proc) \
   std_device_color_stype_body(gx_device_djvu, 0, name, &st_device_djvu, \
                               (int)((long)DEFAULT_WIDTH_10THS*300/10),  \
	   		       (int)((long)DEFAULT_HEIGHT_10THS*300/10), \
				300, 300, 24, 255, 256 ), \
   procs, proc, 80, 256, 256, 3, 10000000 

/* Open the device. */
private int
djvu_open(gx_device * dev)
{
    gx_device_djvu *cdev = (gx_device_djvu *)dev;
    cdev->head.first = 0;
    cdev->head.last = 0;
    cdev->gchrome = 0;
    cdev->curmask = 0;
    cdev->curflags = 0;
    cdev->curcolor = 0;
    cdev->curchrome = 0;
    cdev->curbreak = false;
    cdev->marks = 0;
    cdev->pmem = p2mem_create();
    if (!cdev->pmem)
        return_VMerror;
    return 0;
}

/* Flush current runmapbld into drawlist. */
private int
djvu_flush_current(gx_device_djvu *cdev)
{
    int code = 0;
    /* Proceed */
    if (cdev->curmask) {
        runmap *mask;
        drawlist *last;
        if ((code = runmapbld_close(cdev->curmask, &mask)) < 0)
            return code;
        /* Terminate chrome data if present */
        if (cdev->curchrome) 
            chrome_closepos(cdev->gchrome, cdev->curchrome);
        /* Attempt to merge path operations */
        last = cdev->head.last;
        if (mask && last && !cdev->curchrome 
            && cdev->curcolor == last->color && ! last->text
            && ((cdev->curflags ^ last->flags) & ~DLIST_PATH_FILL) == 0
            && runmap_hittest(mask, last->mask) ) {
            /* Merge with previous */
            last->flags |= cdev->curflags;
            code = runmap_or(cdev->pmem, mask, last->mask, 3, &last->mask);
        } else if (mask) {
            /* Create new element */
            drawlist *dl = drawlist_append(cdev->pmem, &cdev->head);
            if (! dl) return_VMerror;
            dl->flags = cdev->curflags;
            dl->color = cdev->curcolor;
            dl->chromepos = cdev->curchrome;
            dl->mask = mask;
#ifdef DEBUG
            if (cdev->curchrome && gs_debug_c(DEBUG_CHAR_CHROME))
                drawlist_print(dl, cdev);
#endif
        }
    }
    cdev->curflags = 0;
    cdev->curmask = 0;
    cdev->curcolor = 0;
    cdev->curchrome = 0;
    cdev->curbreak = false;
    return code;
}

/* Flushes only if curbreak is set and if 
   the current object bounding box does 
   not touch the specified rectangle WxH+X+Y. */
private int
djvu_flush_maybe(gx_device_djvu *cdev, int x, int y, int w, int h)
{
    int code;
    runmap *r;
    if (cdev->curbreak) {
        cdev->curbreak = false;
        if (cdev->curmask) {
            if ((code = runmapbld_flush(cdev->curmask)) < 0)
                return code;
            if ((r = cdev->curmask->rmap))
                if (x>r->xmax+1 || y>r->ymax+1 || x+w<r->xmin || y+h<r->ymin)
                    return djvu_flush_current(cdev);
        }
    }
    return 0;
}

/* Internal: Setup chrome data for this object */
private int
djvu_check_chrome(gx_device_djvu *cdev, gx_color_index color)
{
    int code;
    if (! cdev->curchrome) {
        /* Mark object has having inpure color */
        cdev->curflags |= DLIST_INPURE;
        /* Create chrome data if needed */
        if (!cdev->gchrome && !(cdev->gchrome = chrome_create()))
            return_VMerror;
        if (! (cdev->curchrome = chrome_getpos(cdev->gchrome)))
            return_VMerror;
        /* Store existing runmap into chrome */
        if (cdev->curmask) {
            runmapbld_flush(cdev->curmask);
            code = chrome_put_runmap(cdev->gchrome, 
                                     cdev->curmask->rmap, cdev->curcolor);
            if (code < 0) return code;
        }
    }
    return 0;
}

/* Make sure that there is a current runmapbld. */
private int
djvu_check_current(gx_device_djvu *cdev, gx_color_index color)
{
    int code = 0;
    cdev->curbreak = false;
    if (! cdev->curmask) {
        /* Start a new object */
        if (color == gx_no_color_index)
            if ((code = djvu_check_chrome(cdev, color)) < 0)
                return code;
        code = runmapbld_open(cdev->pmem, runmapbld_buffer_size, &cdev->curmask);
        cdev->curcolor = color;
        cdev->curflags = 0;
    } else if (color != cdev->curcolor) {
        /* Color is no longer a single color */
        code = djvu_check_chrome(cdev, color);
        cdev->curcolor = gx_no_color_index;
    }
    return code;
}

/* Clean up page data. */
private void
djvu_close_page(gx_device_djvu *cdev)
{
    p2mem_freeall(cdev->pmem);
    chrome_destroy(cdev->gchrome);
    cdev->curmask = 0;
    cdev->head.first = 0;
    cdev->head.last = 0;
    cdev->gchrome = 0;
    cdev->fgpalette = 0;
    cdev->fgpalettesize = 0;
    cdev->marks = 0;
}

/* GS callback: Output a page. */
private int
djvu_output_page(gx_device * dev, int num_copies, int flush)
{
    int code = 0;
    gx_device_djvu *cdev = (gx_device_djvu *)dev;
    code = djvu_flush_current(cdev);
    if (code < 0) return code;
    /* open output file if not already open */
    if (!cdev->outputfile && cdev->outputfilename[0]) {
        code = gx_device_open_output_file(dev, cdev->outputfilename, 
					  true, false, &cdev->outputfile);
	if (code < 0) return code;
    }
    /* process */
    if (!cdev->process) return_error(gs_error_unregistered);
    code = (*cdev->process)(cdev);
    if (code < 0) return code;
    /* cleanup and return */
    if (flush)
        djvu_close_page(cdev);
    if (cdev->reopenperpage) {
        if (cdev->outputfile)
            gx_device_close_output_file(dev, cdev->outputfilename, 
                                        cdev->outputfile);
        cdev->outputfile = 0;
    }
    return gx_finish_output_page(dev, num_copies, flush);
}

/* GS callback: Close the device. */
private int
djvu_close(gx_device * dev)
{
    gx_device_djvu *cdev = (gx_device_djvu *) dev;
#ifndef NO_IMPLICIT_FINAL_SHOWPAGE
    /* Force an implicit final showpage */
    if (cdev->head.first || cdev->curmask) {
        int code = djvu_output_page(dev, 1, true);
        if (code < 0) return code;
    }
#endif
    /* Cleanup */
    djvu_close_page(cdev);
    if (cdev->outputfile)
        gx_device_close_output_file(dev, cdev->outputfilename, 
                                    cdev->outputfile);
    cdev->outputfile = 0;
    p2mem_destroy(cdev->pmem);
    cdev->pmem = 0;
    return 0;
}

/* GS callback: Get parameters. */
private int
djvu_get_params(gx_device * dev, gs_param_list * plist)
{
    int code;
    gx_device_djvu *cdev = (gx_device_djvu *) dev;
    gs_param_string ofns;
    ofns.data = (const byte *)cdev->outputfilename;
    ofns.size = strlen(cdev->outputfilename);
    ofns.persistent = false;
    code = gx_default_get_params(dev, plist);
#define GET(param_write,param_name, place) \
    if (code >= 0) code=param_write(plist, param_name, place)
    GET(param_write_string, "OutputFile", &ofns);
    GET(param_write_int, "Threshold", &cdev->threshold);
    GET(param_write_int, "FgColors", &cdev->fgcolors);
    GET(param_write_int, "FgImgColors", &cdev->fgimgcolors);
    GET(param_write_int, "BgSubsample", &cdev->bgsubsample);
    GET(param_write_long, "MaxBitmap", &cdev->maxbitmap);
    GET(param_write_bool, "ReopenPerPage", &cdev->reopenperpage);
    GET(param_write_bool, "AutoHires", &cdev->autohires);
    GET(param_write_bool, "ExtractText", &cdev->extracttext);
    GET(param_write_bool, "QUIET", &cdev->quiet);
    if (code >= 0) code=param_write_null(plist, "pdfmark");
#undef GET
    return code;
}

/* Internal: Process pdfmarks */
private int
djvu_pdfmark(gx_device_djvu *cdev, gs_param_string_array *pma)
{
    int code = 1;
    int size = pma->size;
    if (size>=1 && pdfmark_eq(&pma->data[size-1], "LNK")) {
        pdfmark *mark = 0;
        code = pdfmark_do_lnk(cdev->pmem, pma, &mark);
        if (code >= 0 && mark) {
            mark->next = cdev->marks;
            cdev->marks = mark;
        }
    } else if (size>=1 && pdfmark_eq(&pma->data[size-1], "ANN")) {
        char *action = 0;
        char *subtype = 0;
        if (pdfmark_find(cdev->pmem, pma, "/Subtype", &subtype) >= 0 &&
            pdfmark_find(cdev->pmem, pma, "/Action", &action) >= 0 &&
            subtype && !strcmp(subtype,"/Link") &&
            action && !strcmp(action,"/Launch") ) {
                pdfmark *mark = 0;
                code = pdfmark_do_lnk(cdev->pmem, pma, &mark);
                if (code >= 0 && mark) {
                        mark->next = cdev->marks;
                        cdev->marks = mark;
                }
        }
        p2mem_free(cdev->pmem, subtype);
        p2mem_free(cdev->pmem, action);
    } else if (size>=1 && pdfmark_eq(&pma->data[size-1], "PAGELABEL")) {
	char *label = 0;
	if (pdfmark_find_recode(cdev->pmem, pma, "/Label", &label) >= 0 && label) {
	    pdfmark *mark = p2mem_alloc(cdev->pmem, sizeof(pdfmark));
	    if (mark) {
		memset(mark, 0, sizeof(pdfmark));
		mark->type = 'P';
		mark->u.p.title = label;
		mark->next = cdev->marks;
		cdev->marks = mark;
		label = 0;
	    }
	}
	p2mem_free(cdev->pmem, label);
    } else if (size>=1 && pdfmark_eq(&pma->data[size-1], "OUT")) {
        pdfmark *mark = 0;
        code = pdfmark_do_out(cdev->pmem, pma, &mark);
        if (code >= 0 && mark) {
            mark->next = cdev->marks;
            cdev->marks = mark;
        }
    }
    return code;
}

/* Internal: Validate outputfilename */
private bool
validate_outputfilename(const char *data, int size, bool *rpp, gs_memory_t *mem)
{
    const char *hasformat = 0;
    gs_parsed_file_name_t parsed;
#if GS_VERSION >= 900
    int code = gx_parse_output_file_name(&parsed, &hasformat, data, size, mem);
#else
    int code = gx_parse_output_file_name(&parsed, &hasformat, data, size);
#endif
    if (code<0 || size>=prn_fname_sizeof)
        return false;
    *rpp = (hasformat ? true: false);
    return true;
}

/* GS callback: Put parameters. */
private int
djvu_put_params(gx_device * dev, gs_param_list * plist)
{
    gx_device_djvu *cdev = (gx_device_djvu *) dev;
    bool was_open = dev->is_open;
    int ncode = 0;
    int code = 0;
    /* Parameters and defaults */
    gs_param_string ofs;
    int thr = cdev->threshold;
    int fgc = cdev->fgcolors;
    int fgi = cdev->fgimgcolors;
    int bgs = cdev->bgsubsample;
    long mbm = cdev->maxbitmap;
    bool rpp = cdev->reopenperpage;
    bool ahi = cdev->autohires;
    bool etx = cdev->extracttext;
    bool quiet = cdev->quiet;
    /* Handle pseudo-parameter pdfmark */
    gs_param_string_array pma;
    if (param_read_string_array(plist, "pdfmark", &pma) == 0) {
        code = djvu_pdfmark(cdev, &pma);
	return (code < 0) ? code : 0;
    }
    /* Validation macro */
#define PUT(param_read, param_name, place, predicate) \
    ncode = param_read(plist, param_name, place); \
    if (ncode == 0 && !(predicate)) { \
        ncode = gs_error_limitcheck; } \
    if (ncode < 0) { code=ncode; \
        param_signal_error(plist, param_name, code); } 
    /* Validate OutputFile */
    ofs.data = 0;
    PUT(param_read_string, "OutputFile", &ofs, 
        validate_outputfilename((const char*)ofs.data, ofs.size, &rpp, cdev->memory));
    /* Validate remaining parms */
    PUT(param_read_int, "Threshold", &thr, (thr>=0 && thr<=100));
    PUT(param_read_int, "FgColors", &fgc, (fgc>=1 && fgc<=4000));
    PUT(param_read_int, "FgImgColors", &fgi, (fgi>=0 && fgi<=4000));
    PUT(param_read_int, "BgSubsample", &bgs, (bgs>=1 && bgs<=12));
    PUT(param_read_long, "MaxBitmap", &mbm, (mbm<=1000000000));
    PUT(param_read_bool, "ReopenPerPage", &rpp, true);
    PUT(param_read_bool, "AutoHires", &ahi, true);
    PUT(param_read_bool, "ExtractText", &etx, true);
    PUT(param_read_bool, "QUIET", &quiet, true);
#undef PUT
    /* Terminate validation */
    if (code >= 0) {
	cdev->is_open = false;
	code = gx_default_put_params(dev, plist);
	cdev->is_open = was_open;
    }
    if (code < 0)
        return code;
    /* Install parameters */
    cdev->threshold = thr;
    cdev->fgcolors = fgc;
    cdev->fgimgcolors = fgi;
    cdev->bgsubsample = bgs;
    cdev->maxbitmap = max(1000000, mbm);
    cdev->reopenperpage = rpp;
    cdev->autohires = ahi;
    cdev->extracttext = etx;
    cdev->quiet = quiet;
    /* Install ``OutputFile'' */
    if (ofs.data) {
        if (ofs.size != strlen(cdev->outputfilename) ||
            memcmp(ofs.data, cdev->outputfilename, ofs.size) ) {
            /* Modify outputfilename */
            if (was_open)  
                discard(gs_closedevice(dev));
            memcpy(cdev->outputfilename, ofs.data, ofs.size);
            cdev->outputfilename[ofs.size] = 0;
            if (was_open)  
                code = gs_opendevice(dev);
        }
    }
    return code;
}

/* ------ low-level rendering procs ------ */

/* GS callback: fill a rectangle with solid color */
private int
djvu_fill_rectangle(gx_device *dev, int x, int y, int w, int h,
		    gx_color_index color)
{
    int code;
    gx_device_djvu * const cdev = (gx_device_djvu *)dev;
    /* Ignore initial white fills */
    if (!cdev->head.last && !cdev->curmask)
        if (color == gx_device_white(dev))
            return 0;
    /* Check for break */
    fit_fill(dev, x, y, w, h);
    if (cdev->curbreak) {
        code = djvu_flush_maybe(cdev, x, y, w, h);
        if (code < 0) return code;
    }
    /* Accumulate into current mask */
    cdev->curflags |= color_flags(color);
    if ((code = djvu_check_current(cdev, color)) < 0)
        return code;
    if ((code = runmapbld_put_rect(cdev->curmask, x, y, w, h)) < 0)
        return code;
    /* Save chrome info */
    if (cdev->curchrome)
        if ((code = chrome_put_rect(cdev->gchrome, x, y, w, h, color)) < 0)
            return code;
    return 0;
}

/* GS callback: copy bitmap with specified fg/bg colors */
private int
djvu_copy_mono(gx_device *dev,
	       const byte * base, int sourcex, int sraster, gx_bitmap_id id,
	       int x, int y, int w, int h, 
	       gx_color_index color0, gx_color_index color1)
{
    int code;
    gx_device_djvu * const cdev = (gx_device_djvu *)dev;
    /* Check for break */
    fit_copy(dev, base, sourcex, sraster, id, x, y, w, h);
    if (cdev->curbreak) {
        code = djvu_flush_maybe(cdev, x, y, w, h);
        if (code < 0) return code;
    }
    /* Attempt to do something smart about scanned bilevel images */
    if (!cdev->head.last) {
        /* We test that (a) this is the first page component, (b) this
           component has been drawn with pure black color, and (c) this image
           is drawn below all other drawings (banding). */
        bool maybescan = true;
        if (cdev->curcolor != 0x000000)
            maybescan = false;
        if (cdev->curmask) {
            runmapbld *bld = cdev->curmask;
            if (bld->rmap && y<=bld->rmap->ymax)
                maybescan = false;
            else if (bld->nruns>0 && y<=bld->runs[bld->nruns-1][0])
                maybescan = false;
        }                
        if (maybescan) {
            /* Convert b&w image to imagemask */
            if (color0 == 0xffffff && color1 == 0x000000)
                color0 = gx_no_color_index;
            if (color1 == 0xffffff && color0 == 0x000000)
                color1 = gx_no_color_index;
        }
    }
    /* Three cases really matter */
    if (color0 == gx_no_color_index) {
        /* One color */
        cdev->curflags |= color_flags(color1);
        if ((code = djvu_check_current(cdev, color1)) >= 0)
            code = runmapbld_put_bitmap(cdev->curmask, base, sourcex, sraster, 
					0x00, x, y, w, h);
    } else if (color1 == gx_no_color_index) {
        /* One color, reverse video */
        cdev->curflags |= color_flags(color0);
        if ((code = djvu_check_current(cdev, color0)) >= 0)
            code = runmapbld_put_bitmap(cdev->curmask, base, sourcex, sraster, 
					0xff, x, y, w, h);
    } else {
        /* Two colors */
        cdev->curflags |= color_flags(color0);
        cdev->curflags |= color_flags(color1);
        if ((code = djvu_check_current(cdev, gx_no_color_index)) >= 0)
            code = runmapbld_put_rect(cdev->curmask, x, y, w, h);
    }
    if (code < 0) 
        return code;
    /* Save chrome info */
    if (cdev->curchrome)
        if ((code = chrome_put_bitmap(cdev->gchrome, base, sourcex, sraster, 
                                      x, y, w, h, color0, color1)) < 0)
            return code;
    return 0;
}

/* GS callback: copy color pixmap */
private int
djvu_copy_color(gx_device *dev,
                   const byte * base, int sourcex, int sraster, gx_bitmap_id id,
                   int x, int y, int w, int h)
{
    gx_device_djvu * const cdev = (gx_device_djvu *)dev;
    byte hascolor = 0;
    int code; 
    /* Check for break */
    fit_copy(dev, base, sourcex, sraster, id, x, y, w, h);
    if (cdev->curbreak) {
        code = djvu_flush_maybe(cdev, x, y, w, h);
        if (code < 0) return code;
    }
    /* Accumulate mask */
    code = djvu_check_current(cdev, gx_no_color_index);
    if (code < 0) return code;
    code = runmapbld_put_rect(cdev->curmask, x, y, w, h);
    if (code < 0) return code;
    /* Save chrome info */
    if (cdev->curchrome)
        if ((code = chrome_put_pixmap(cdev->gchrome, base, sourcex, sraster,
                                      x, y, w, h )) < 0)
            return code;
    /* Color check */
    base = base + 3*sourcex;
    while (!hascolor && h-- > 0) {
        int x = 0;
        const byte *row = base;
        while (x++ < w) {
            byte r = *row++;
            byte g = *row++;
            byte b = *row++;
            hascolor |= ((r^g)|(g^b));
        }
        base += sraster;
    }
    cdev->curflags |= (DLIST_HASGRAY|DLIST_INPURE);
    if (hascolor) cdev->curflags |= DLIST_HASCOLOR;
    return 0;
}

/* ------ path procs ------ */

/* GS callback: postscript operator ``fill'' */
private int
djvu_fill_path(gx_device *dev, const gs_imager_state * pis,
	       gx_path * ppath, const gx_fill_params * params,
	       const gx_device_color *pdcolor, const gx_clip_path * pcpath)
{
    int code;
    gx_device_djvu * const cdev = (gx_device_djvu *)dev;
    gx_color_index color = gx_no_color_index;
    /* Sometimes high level subroutines call each other */
    if (cdev->curflags & DLIST_RECURSIVE)
        return gx_default_fill_path(dev, pis, ppath, params, pdcolor, pcpath);
    /* Segment driver stream */
    code = djvu_flush_current(cdev);
    if (code < 0) return code;
    if (  gx_dc_writes_pure(pdcolor, pis->log_op) 
	  && pis->alpha==gx_max_color_value)
        color = gx_dc_pure_color(pdcolor);
    code = djvu_check_current(cdev, color);
    if (code < 0) return code;
    /* Call default routine */
    cdev->curflags |= DLIST_PATH | DLIST_PATH_FILL | DLIST_RECURSIVE;
    code = gx_default_fill_path(dev, pis, ppath, params, pdcolor, pcpath);
    cdev->curflags &= ~DLIST_RECURSIVE;
    if (code < 0) return code;
    return djvu_flush_current(cdev);
}

/* GS callback:  postscript operator ``stroke'' */
private int
djvu_stroke_path(gx_device *dev, const gs_imager_state * pis,
		 gx_path * ppath, const gx_stroke_params * params,
		 const gx_drawing_color *pdcolor, const gx_clip_path * pcpath)
{
    int code;
    gx_device_djvu * const cdev = (gx_device_djvu *)dev;
    gx_color_index color = gx_no_color_index;
    /* Sometimes high level subroutines call each other */
    if (cdev->curflags & DLIST_RECURSIVE)
        return gx_default_stroke_path(dev, pis, ppath, params, pdcolor, pcpath);
    /* Segment driver stream */
    code = djvu_flush_current(cdev);
    if (code < 0) return code;
    if (  gx_dc_writes_pure(pdcolor, pis->log_op) 
	  && pis->alpha==gx_max_color_value)
        color = gx_dc_pure_color(pdcolor);
    code = djvu_check_current(cdev, color);
    if (code < 0) return code;
    /* Call default routine */
    cdev->curflags |= DLIST_PATH | DLIST_RECURSIVE;
    code = gx_default_stroke_path(dev, pis, ppath, params, pdcolor, pcpath);
    cdev->curflags &= ~DLIST_RECURSIVE;
    return djvu_flush_current(cdev);
}

/* ------ image procs ------ */

/* GS callback:  begin processing an image */
private int 
djvu_begin_typed_image(gx_device *dev,
                       const gs_imager_state *pis, const gs_matrix *pmat,
                       const gs_image_common_t *pim, const gs_int_rect *prect,
                       const gx_drawing_color *pdcolor,
                       const gx_clip_path *pcpath,
                       gs_memory_t *memory, gx_image_enum_common_t **pinfo)
{
    gx_device_djvu * const cdev = (gx_device_djvu *)dev;
    /* Segment driver stream (only when the call is non recursive) */
    if (! (cdev->curflags & DLIST_RECURSIVE)) {
        /* We should here call djvu_flush_current() in order to keep images
           separate.  There are many crazy pdf files however that encode
           images using a number of small images at different resolutions.
           In order to avoid breaking these, we just set a flag that will
           cause the low level routines to call djvu_flush_maybe(). */
        cdev->curbreak = true;
    }
    /* Note: we cannot easily set DLIST_RECURSIVE around image drawing 
       callsbecause we should then override the image procedure vector. */
    return gx_default_begin_typed_image(dev, pis, pmat, pim, prect,
                                        pdcolor, pcpath, memory, pinfo);
}

/* ------ text procs ------ */

/* We want to delimit the text operations by executing custom code
   before and after running function gs_text_enum_procs_t::process().

   PROBLEM --- The ``normal'' solution (cf. gdevbbox.c) consists in having
   djvu_text_begin() create a subclass of gs_text_enum_t that delegates its
   work to another gs_text_enum_t created with the default function.  This is
   illustrated in gdevbbox and gdevpdft.  Although this solution basically
   works, it seems that something in the ghostscript font renderer works
   differently when the text enumerator is not a gs_show_enum_t.

   SOLUTION --- The default text wrapper is based on ugly hack.  Function
   djvu_text_begin() creates a default text enumerator and replaces its
   statically allocated procedure vector by a dynamically allocated procedure
   vector containing the usual procedure, plus a pointer to the original
   procedures, plus whatever additional fields are needed by the wrapper.  It
   works as long as the default text enumerator does not extend the procedure
   vector itself.  The fat procedure vector and anything inside must be
   allocated with an allocator without garbage collection.  Otherwise the
   garbage collector would not be able to properly mark these objects and
   would reclaim the memory at once.  
*/

typedef struct fat_text_enum_procs_s {
    gs_text_enum_procs_t procs;
    const gs_text_enum_procs_t *origprocs;
    void (*origfree)(gs_memory_t*, void*, client_name_t);
    struct fat_text_enum_procs_s *myself;
    p2mem *pmem;
    gs_state *pgs;
    gs_int_point lastpoint;
    bool lastpointvalid;
} fat_text_enum_procs_t;


/* Utility: finds the current point in image coordinates */
private bool
djvu_currentpoint(gs_state *pgs, gs_int_point *ipoint)
{
    gs_matrix ctm;
    gs_point point, tpoint;
    if (gs_currentpoint(pgs, &point) >= 0 &&
        gs_currentmatrix(pgs, &ctm) >= 0) {
        gs_point_transform(point.x, point.y, &ctm, &tpoint);
        ipoint->x = (int)floor(tpoint.x + 0.5);
        ipoint->y = (int)floor(tpoint.y + 0.5);
        return true;
    }
    return false;
}


/* Utility: finds what character we are displaying */
private void
djvu_text_extract(gs_text_enum_t *pte)
{
    fat_text_enum_procs_t *fat = (fat_text_enum_procs_t *)pte->procs;
    gx_device_djvu * const cdev = (gx_device_djvu *)pte->dev;
    gs_font *font = gs_text_current_font(pte);
    gs_glyph glyph = gs_text_current_glyph(pte);
    gs_char unicode = GS_NO_CHAR;
    gs_int_point lastpoint = fat->lastpoint;
    bool lastpointvalid = fat->lastpointvalid;
    gs_int_point point;
    bool pointvalid = djvu_currentpoint(fat->pgs, &point);
    if (font && glyph) {
#if GS_VERSION >= 903
        unicode = font->procs.decode_glyph(font, glyph, -1);
#else
        unicode = font->procs.decode_glyph(font, glyph);
#endif
    }
    if (unicode != GS_NO_CHAR && lastpointvalid && pointvalid) {
        txtmark *mark = 0;
        char *arg = 0;
        byte buffer[8];
        int fcode;
        int len;
        /* flush current drawing */
        if (cdev->curmask) {
            gx_color_index color;
            cdev->curflags &= ~DLIST_RECURSIVE;
            fcode = djvu_flush_current(cdev);
            len = unicode_to_utf8(unicode, buffer);
            if (len > 0 && fcode >= 0) {
                drawlist *dl = cdev->head.last;
                mark = p2mem_alloc(cdev->pmem, sizeof(txtmark));
                arg = make_ps_string(cdev->pmem, buffer, len);
                if (dl && dl->mask && !dl->text && mark && arg) {
                    mark->x = lastpoint.x;
                    mark->y = lastpoint.y;
                    mark->w = point.x - mark->x;
                    mark->h = point.y - mark->y;
                    mark->arg = arg;
                    dl->text = mark;
                    mark = 0;
                    arg = 0;
                }
                if (mark)
                    p2mem_free(cdev->pmem, mark);
                if (arg)
                    p2mem_free(cdev->pmem, arg);            
            }
            /* prepare new context */
            color = gx_no_color_index;
            if (  pte->pdcolor
                  && gx_dc_writes_pure(pte->pdcolor, pte->pis->log_op)
                  && pte->pis->alpha==gx_max_color_value )
                color = gx_dc_pure_color(pte->pdcolor);
            djvu_check_current(cdev, color);
            cdev->curflags |= DLIST_TEXT | DLIST_RECURSIVE;
        }
    }
    if (pointvalid) {
        fat->lastpoint = point;
        fat->lastpointvalid = pointvalid;
    }
}


/* GS finalization callback */
private void
djvu_text_free(gs_memory_t * memory, void *vpte, client_name_t cname)
{
    gs_text_enum_t *pte = (gs_text_enum_t *)vpte;
    fat_text_enum_procs_t *fat = (fat_text_enum_procs_t *)pte->procs;
    gx_device_djvu * const cdev = (gx_device_djvu *)pte->dev;
    if (cdev->curmask)
        if (cdev->curflags & DLIST_TEXT) {
            cdev->curflags &= ~DLIST_RECURSIVE;
            djvu_flush_current(cdev);
        }
    ASSERT(fat == fat->myself);
    pte->procs = fat->origprocs;
    pte->rc.free = fat->origfree;
    fat->origfree(memory, vpte, cname);
    p2mem_free(fat->pmem, fat->myself);
}

/* GS callback: process text */
private int
djvu_text_process(gs_text_enum_t * pte)
{
    fat_text_enum_procs_t *fat = (fat_text_enum_procs_t *)pte->procs;
    gx_device_djvu * const cdev = (gx_device_djvu *)pte->dev;
    int code;
    int oper;
    /* delegate */
    ASSERT(fat == fat->myself);
    oper = pte->text.operation;
    if (cdev->extracttext)
        pte->text.operation |= TEXT_INTERVENE;
    do {
        code = fat->origprocs->process(pte);
        if (cdev->extracttext)
            if (code == TEXT_PROCESS_INTERVENE || code == 0)
                djvu_text_extract(pte);
    } while ((code==TEXT_PROCESS_INTERVENE) && !(oper & TEXT_INTERVENE));
    if (cdev->extracttext && !(oper & TEXT_INTERVENE))
        pte->text.operation &= ~TEXT_INTERVENE;
    /* is text operation finished */
    if (code <= 0) {
        int fcode;
        cdev->curflags &= ~DLIST_RECURSIVE;
        fcode = djvu_flush_current(cdev);
        if (fcode<0 && code>=0) 
            code = fcode;
    }
    return code;
}

/* GS callback: start processing text */
private int
djvu_text_begin(gx_device *dev, gs_imager_state * pis,
		const gs_text_params_t * text, gs_font * font,
		gx_path * path, const gx_device_color * pdcolor,
		const gx_clip_path * pcpath,
		gs_memory_t * mem, gs_text_enum_t ** ppenum)
{
    gx_device_djvu * const cdev = (gx_device_djvu *)dev;
    fat_text_enum_procs_t *fat;
    gs_text_enum_t *pte;
    gx_color_index color;
    int code;
    /* Create default text enumerator */
    code = gx_default_text_begin(dev, pis, text, font, path, 
                                 pdcolor, pcpath, mem, ppenum );
    pte = *ppenum;
    /* Return immediately when */
    if (code < 0) 
        return code; /* failed */
    if (cdev->curflags & DLIST_RECURSIVE)
        return code; /* called within another text operation */
    if (! (pte->text.operation & TEXT_DO_DRAW))
        return code; /* not a drawing operation */
    /* Flush current object */
    code = djvu_flush_current(cdev);
    if (code < 0) return code;
    /* Create fat procedure vector */
    if (! (fat = p2mem_alloc(cdev->pmem, sizeof(fat_text_enum_procs_t))))
	return_VMerror;
    fat->procs = *pte->procs;
    fat->origprocs = pte->procs;
    fat->origfree = pte->rc.free;
    fat->procs.process = djvu_text_process;
    fat->myself = fat;
    fat->pmem = cdev->pmem;
    fat->pgs = (gs_state*)pis;
    fat->lastpointvalid = false;
    if (cdev->extracttext)
        if (djvu_currentpoint((gs_state*)pis, &fat->lastpoint))
            fat->lastpointvalid = true;
    /* Prepare text object */
    color = gx_no_color_index;
    if (  pte->pdcolor
          && gx_dc_writes_pure(pte->pdcolor, pte->pis->log_op)
	  && pte->pis->alpha==gx_max_color_value )
        /* Pure color so far */
        color = gx_dc_pure_color(pte->pdcolor);
    code = djvu_check_current(cdev, color);
    cdev->curflags |= DLIST_TEXT | DLIST_RECURSIVE;
    /* Mutate returned pte */
    pte->rc.free = djvu_text_free;
    pte->procs = (gs_text_enum_procs_t *)fat;
    return code;
}



/* ======================================
      L O W   C O L O R   I M A G E S
   ====================================== */

/* Drawlist components with multiple colors usually represent images.  Some
   images have a small number of colors (especially if they were encoded using
   an indexed color model).  The following code examines these images and
   sometimes breaks them into pure color components. */

/* Internal: size of per-color runmapbld */
#define lowcolor_bldsz 170

/* Internal: subclassed payload for color hash table */
typedef struct lowcolordata_s {
    colordata cdata;
    runmapbld *bld;
    runmap *map;
} lowcolordata;

/* Internal: release function for lowcolordata */
private void
lowcolordata_release(p2mem *mem, void *payload)
{
    lowcolordata *cd = payload;
    runmapbld_close(cd->bld, 0);
    runmap_free(cd->map);
}

/* Internal: sorting function for lowcolordata pointers */
private int
lowcolordata_sortsub(const void *a, const void *b)
{
    const lowcolordata *cda = *(lowcolordata * const *)a;
    const lowcolordata *cdb = *(lowcolordata * const *)b;
    return cdb->cdata.w - cda->cdata.w;
}

/* Internal: color distance test */
private inline int
lowcolor_distance(lowcolordata *cd1, lowcolordata *cd2)
{
    gx_color_index c1 = cd1->cdata.color;
    gx_color_index c2 = cd2->cdata.color;
    if ((c1 ^ c2) & 0xF0F8E0) { /* quick filter */
        int dr = (byte)(c1>>16) - (byte)(c2>>16);
        int dg = (byte)(c1>>8) - (byte)(c2>>8);
        int db = (byte)(c1) - (byte)(c2);
        /* Following difference range from 0 to 0xf */
        return ( 5*dr*dr + 9*dg*dg + 2*db*db ) >> 16;
    }
    return 0;
}

/* Internal: perimeter cost (threshold adjustment) */
#define lowcolor_perimeter_cost 8

/* Internal: play the chrome data, build a colormap,
   compute various perimeters, build runmap for each 
   color drawlist component.  Returns a positive
   status if the image must remain in the background. */
private int
lowcolor_separate(gx_device_djvu *cdev, drawlist *dl, 
                  htable *colorhash, int maxcolors)
{
    int code;
    p2mem *mem = cdev->pmem;
    lowcolordata **band = 0;
    int size;
    int w, h, xb, yb;
    lowcolordata *cd;
    int raster;
    int nbands;
    int bandh;
    uint perim;
    uint hardperim;
    /* Determine band size and allocate buffer */
    w = dl->mask->xmax - dl->mask->xmin + 1;
    h = dl->mask->ymax - dl->mask->ymin + 1;
    raster = w;
    nbands = 1;
    for (;;) {
        bandh = (h + nbands - 1) / nbands;
        size = raster * bandh * sizeof(lowcolordata*);
        if (size <= cdev->maxbitmap)
            if ((band = p2mem_alloc(mem, size)))
                break;
        if (bandh < 16)
            return 1;
        nbands += 1;
    }
    /* Compute perimeters */
    if (dl->mask && !dl->area)
        runmap_info(dl->mask, &dl->area, &dl->perimeter);
    perim = dl->perimeter * lowcolor_perimeter_cost;
    hardperim = dl->perimeter * 0xf;
    xb = dl->mask->xmin;
    for (yb = dl->mask->ymin; yb <= dl->mask->ymax; yb += bandh) {
        int y;
        int hb = min(bandh, dl->mask->ymax - yb + 1);
        lowcolordata **row;
        /* render band using indexed model */
        memset(band, 0, size);
        code = chrome_play_indexed(dl->chromepos, colorhash, maxcolors,
                                   (colordata**)band, raster, 
                                   xb, yb, w, hb);
        if (code != 0) goto xit;
        /* compute perimeter data */
        row = band;
        for (y = yb; y < yb+hb; y++, row+=raster) {
            int x;
            lowcolordata **r = row;
            lowcolordata *lcd = 0;
            for (x=0; x<w; x++,r++) {
                lowcolordata *cd = r[0];
                if (cd!=lcd && cd && lcd) {
                    perim += lowcolor_perimeter_cost;
                    hardperim += lowcolor_distance(cd,lcd);
                }
                if (y > yb) {
                    lcd = r[-raster];
                    if (cd!=lcd && cd && lcd) {
                        perim += lowcolor_perimeter_cost;
                        hardperim += lowcolor_distance(cd,lcd);
                    }
                }
                lcd = cd;
            }
        }
    }
    /* Perform test on perimeter ratio */
    code = 1;
#ifdef DEBUG
    if (gs_debug_c(DEBUG_CHAR_LOCOLOR))
        fprintf(STDOUT,"LC: perim=%d hardperim=%d ncolors=%d\n", 
                perim, hardperim, colorhash->nelems);
#endif
    while (perim >= 0x1000000) { perim >>= 1; hardperim >>= 1; }
    if ( hardperim * 100 <= perim * (100 - cdev->threshold) )
        goto xit;
    /* Separate color runmaps */
    for (yb = dl->mask->ymin; yb <= dl->mask->ymax; yb += bandh) {
        int y;
        int hb = min(bandh, dl->mask->ymax - yb + 1);
        lowcolordata **row;
        /* render band (already there if there is only one band) */
        if (nbands > 1) {
            memset(band, 0, size);
            code = chrome_play_indexed(dl->chromepos, colorhash, maxcolors,
                                       (colordata**)band, raster, 
                                       xb, yb, w, hb);
            if (code != 0) goto xit;
        }
        /* extract runs */
        row = band;
        for (y = yb; y < yb+hb; y++, row+=raster) {
            int x = 0;
            lowcolordata **r = row;
            while (x < w) {
                int rx = x;
                lowcolordata *lcd = r[0];
                do { x++; r++; } while (x < w && r[0] == lcd);
                if (lcd) {
                    if (!lcd->bld) {
                        code = runmapbld_open(mem, lowcolor_bldsz, &lcd->bld);
                        if (code < 0) goto xit;
                    }
                    code = runmapbld_write(lcd->bld, y, xb+rx, xb+x-1);
                    if (code < 0) goto xit;
                }
            }
        }
    }
    /* Close runmap builders */
    htable_begin_loop(cd, colorhash) {
        code = runmapbld_close(cd->bld, &cd->map);
        if (code < 0) goto xit;
        runmap_info(cd->map, (int*)(&cd->cdata.w), 0); 
        cd->bld = 0;
    } htable_end_loop(cd, colorhash);
    /* Success */
    code = 0;
 xit:
    p2mem_free(mem, band);
    return code;
}

/* Internal: process one lowcolor image */
private int
lowcolor_process_one(gx_device_djvu *cdev, drawlist *dl)
{
    int code;
    p2mem *mem = cdev->pmem;
    htable *colorhash = 0;
    int cdmapsize;
    lowcolordata **cdmap = 0;
    lowcolordata *cd;
    drawlist *nextdl;
    int ncolors;
    int i;
    /* Check */
    ASSERT(dl->mask);
    ASSERT(dl->chromepos);
#ifdef DEBUG
    /* Number of colors is small enough */
    if (gs_debug_c(DEBUG_CHAR_LOCOLOR)) {
        fprintf(STDOUT, "LC: examining");
        drawlist_print(dl, cdev);
    }
#endif
    /* Create color histogram */
    colorhash = htable_alloc(mem, sizeof(lowcolordata));
    if (!colorhash) goto xitmem;
    colorhash->payload_release = lowcolordata_release;
    /* Separate colors */
    code = lowcolor_separate(cdev, dl, colorhash, cdev->fgimgcolors); 
    if (code != 0) goto xit;
    ncolors = colorhash->nelems;
    if (ncolors < 1) goto xit;
    /* Select colors to be moved into foreground layer */
    cdmapsize = 0;
    if (! (cdmap = p2mem_alloc(mem, ncolors*sizeof(lowcolordata*))))
        goto xitmem;
    htable_begin_loop(cd, colorhash) {
        if (cd->map) cdmap[cdmapsize++] = cd;
    } htable_end_loop(cd, colorhash);
    if (cdmapsize < 1) goto xit;
    /* Sort separated color masks by area */
    qsort(cdmap, cdmapsize, sizeof(lowcolordata*), lowcolordata_sortsub);
    /* Transform the original component with the heaviest color */
    dl->flags = dl->flags & ~DLIST_INPURE;
    dl->color = cdmap[0]->cdata.color;
    dl->chromepos = 0;
#ifdef DEBUG
    if (gs_debug_c(DEBUG_CHAR_LOCOLOR)) {
        fprintf(STDOUT,"LC: becomes");
        drawlist_print(dl, cdev);
    }
#endif      
    /* Insert all remaining components */
    nextdl = dl->next;
    for (i=1; i<cdmapsize; i++) {
        drawlist *newdl = drawlist_insert(mem, &cdev->head, nextdl);
        if (!newdl) goto xitmem;
        newdl->flags = (dl->flags | DLIST_FOREGROUND) & ~DLIST_INPURE;
        newdl->color = cdmap[i]->cdata.color;
        newdl->mask = cdmap[i]->map;
        cdmap[i]->map = 0;
#ifdef DEBUG
        if (gs_debug_c(DEBUG_CHAR_LOCOLOR)) {
            fprintf(STDOUT,"LC: plus");
            drawlist_print(newdl, cdev);
        }
#endif
    }
    /* Success */
    code = 0;
    goto xit;
 xitmem:
    /* Memory error */
    code = gs_note_error(gs_error_VMerror);
 xit:
    /* Cleanup and return */
    p2mem_free(mem, cdmap);
    htable_free(colorhash);
    return code;
}

/* Main processing function for lowcolor images */
private int
lowcolor_process(gx_device_djvu *cdev)
{
    if (cdev->fgimgcolors > 0) {
        drawlist *dl;
        drawlist *dlnext = cdev->head.first;
        while ((dl = dlnext)) {
            int code;
            dlnext = dl->next;
            if (dl->chromepos && dl->mask) 
                if ((code = lowcolor_process_one(cdev, dl)) < 0)
                    return code;
        }
    }
    return 0;
}





/* ======================================
               D R A W L I S T
         C L A S S I F I C A T I O N
   ====================================== */

/* Each drawlist record must be classified as either a foreground or a
   background object.  This can become a highly non trivial task when a lot of
   objects overlap.  Imagine for instance a road on a map.  The road might be
   drawn by first drawing a fat black segment, and then drawing a thinner red
   segment over the center of the black segment.  The resulting image shows a
   red road surrounded by two black edges.  The edges result from the
   occlusion of the fat black segment by the thinner red segment.  They
   logically belong to the black segment and are drawn before the red segment.
   But they should be encoded as foreground objects...

   The classification algorithm is based on the minimum description length
   principle.  Each decision affects the number of bits necessary to encode
   the image.  Each decision also affects the quality of the resulting image,
   represented by the number of bits necessary to encode the difference
   between the encoded image and the original image.

   We first assume that all objects are classified as background.  We then
   proceed from the topmost object towards the bottommost object and estimate
   whether classifying the object as foreground would reduce the bit costs of
   the overall image.  Although greedy optimization is not perfect, we can
   obtain satisfactory results (more on this later).

   - An object classified as foreground will be encoded very accurately.  The
     number of bits required by this encoding is roughly proportional to the
     perimeter of the object after removing its occluded parts.  The number of
     bits to describe the residual differences is zero.

   - An object classified as background will be encoded using a lossy scheme
     and less resolution.  Both the encoding and the difference bit costs are
     roughly proportional to the length of the perimeter, *except* for those
     parts of the perimeter that result from occlusions by foreground objects
     located above the current object.  Furthermore, the proportionality
     coefficient depends on the color differences along the object boundary.

   We are then left with several practical problems:

   - Computing the perimeter of a runmap.  This is achieved very efficiently
     by function runmap_info().  The basic idea consists first in adding the
     perimeters of each run considered in isolation, and secondly to remove
     twice the length of the contact lines between runs located on adjacent
     rows.

   - Computing the length of the part of the perimeter which does not result
     from an occlusion by objects located above.  First we compute the
     occlusion as the intersection of the initial mask with a runmap
     representing all objects located above.  Then we compute the final mask
     by removing the occlusion using runmap_andnot().  The following equality
     allows us to efficiently compute the relevant perimeter lengths:
          perimeter(initial) + perimeter(occlusion) = 
              perimeter(final) + 2 * masked_part_of_inital_perimeter

   - Evaluating the color differences along the object boundaries that are not
     arising from an occlusion.  An exact computation would be highly
     time-consuming.  Therefore we make the following simplifications: (a) we
     compute the largest color difference only (instead of averaging the color
     differences along the boundary), and (b) we only compute the color
     differences with respect to a small number of background candidates
     preselected during a preliminary pass on the drawlist.

     This approximation works because it undoes some nasty side effects of the
     greedy optimization procedure.  At each time indeed, we assume that the
     objects located below the current object are all going to be background
     objects.  When we compute the bit costs for encoding an object as
     background, we should also consider that it touches objects located
     below, which eventually might be classified as foreground objects (c.f.
     the red road with black edges).  We should not consider these boundaries
     when evaluating the cost of encoding the current object as background,
     but we cannot do it in a single pass because we do not know yet how these
     objects will turn out.  The above approximation somehow works like an
     oracle that guesses which objects will turn out being background objects. */


/* Must be called before classify_drawlist. 
   Spot good background candidates and link them via pointer bglink.
   This function searches large path fills that define a local background
   color. Every drawlist element can access potential background elements
   located below using the bglink linked list. */
private int
identify_bg_candidates(gx_device_djvu *cdev)
{
    /* Educated guesses */
    int min_fill = 3;
    int min_area = cdev->HWResolution[0] * cdev->HWResolution[1] / 4;
    drawlist *dl = cdev->head.last;
    /* Search candidate background fills */
    drawlist *kdl = dl;
    while (dl) {
        /* Search */
        while (dl) {
            if (dl->mask && !dl->area)
                runmap_info(dl->mask, &dl->area, &dl->perimeter);
            if (dl->mask && (dl->flags & DLIST_PATH_FILL)) {
                uint bboxw = dl->mask->xmax - dl->mask->xmin + 1;
                uint bboxh = dl->mask->ymax - dl->mask->ymin + 1;
                if (dl->area > min_area)
                    if (dl->area * min_fill >= bboxw * bboxh)
                        break; /* Got a bg candidate */
            }
            dl = dl->prev;
        }
#ifdef DEBUG
        if (gs_debug_c(DEBUG_CHAR_DLIST) && dl) {
            fprintf(STDOUT,"bgcan ");
            drawlist_print(dl, cdev);
        }
#endif
        /* Link */
        while (kdl != dl) {
            kdl->bglink = dl;
            kdl = kdl->prev;
        }
        /* Proceed */
        dl = (dl ? dl->prev : dl);
    }
    return 0;
}

/* Internal: Color distance */
private inline int
color_distance(gx_color_index c1, gx_color_index c2)
{
    int dr = (byte)(c1>>16) - (byte)(c2>>16);
    int dg = (byte)(c1>>8) - (byte)(c2>>8);
    int db = (byte)(c1) - (byte)(c2);
    /* Following difference range from 0 to 0xfffff */
    int d = 5*dr*dr + 9*dg*dg + 2*db*db; 
    /* Compute final distance. Very small distances are hardly
       visible and actually must be considered as small.  Beyond that,
       the mere fact that items were drawn using different operators
       indicates that their border probably is perceptually more
       meaningful than their color indicate.  The magic formula below
       attempts to quantify these facts. */
    return min(255, (d >> 8));
}

/* Internal: Return the maximal color distance between the specified object
   and the potential background candidates. */
private int
scan_bg_candidates(p2mem *mem, drawlist *dl, int *pdist)
{
    drawlist *kl = dl->bglink;
    runmap *mask = dl->mask;
    int code = 0;
    int dist = 0;
    int del = 0;
    /* Scan bg candidates */
    while (kl && mask) {
        /* This code sucks because it considers all candidate background
           objects that are overlapped by the visible part of the current
           object.  It should only consider those candidates intersecting the
           part of the boundary that does not result from an occlusion by a
           foreground object.  Possible but tough.  */
        if (runmap_hittest(mask, kl->mask)) {
            int ndist = 160;
            if (kl->color != gx_no_color_index) 
                ndist = color_distance(dl->color, kl->color);
            dist = max(dist, ndist);
            code = runmap_andnot(mem, mask, kl->mask, del, &mask);
            if (code < 0) return code;
            del = 1;
        }
        kl = kl->bglink;
    }
    /* Compare with default white background */
    if (mask) {
        int ndist = color_distance(dl->color, 0xffffff);
        dist = max(dist, ndist);
        if (del) runmap_free(mask);
    }
    /* Return final distance.  */
    *pdist = dist;
    return 0;
}

/* Scan the drawing list and set flag DLIST_FOREGROUND or DLIST_BACKGROUND in
   each element.  The mask of each foreground element is clipped in order to
   only keep the visible part of the element in either the background or
   foreground layer.  Returns a runmap representing the pixels classified as
   foreground. Also returns the accrued flags describing the colors of the
   foreground and background pixels. */
private int
classify_drawlist(p2mem *mem, gx_device_djvu *cdev,
                  runmap **fgmapout, uint *fgflagsout, uint *bgflagsout )
{
    int code = 0;
    runmap *fgmap = 0;
    runmap *fgmapi = 0;
    uint fgflags = 0;
    runmap *bgmap = 0;
    runmap *bgmapi = 0;
    uint bgflags = 0;
    drawlist *dl = cdev->head.last;
    int nthreshold = 100 - cdev->threshold;
#ifdef DEBUG
    bool debug = gs_debug_c(DEBUG_CHAR_DLIST);
#endif
    /* Iterate */
    while (dl) {
        drawlist *dlprev = dl->prev;
        runmap *occlusion = 0;
        runmap *clipped = 0;
	runmap *tmp = 0;
        int dist = -1;
        int initial_perim;
        int occlusion_perim;
        int clipped_perim;
        int bg_perim;
	/* Objects with inpure color or already marked as background */
	if (dl->flags & (DLIST_INPURE|DLIST_BACKGROUND)) {
            dl->flags &= ~DLIST_FOREGROUND;
	    dl->flags |= DLIST_BACKGROUND;
	    goto cleanup;
	}
        /* Determine occlusions by background objects */
        if ( (code = runmap_and(mem, dl->mask, bgmapi, 0, &occlusion)) < 0 ||
             (code = runmap_and(mem, dl->mask, bgmap, 0, &tmp)) < 0        ||
             (code = runmap_or(mem, tmp, occlusion, 3, &occlusion)) < 0     )
            return code;
        /* Applies background occlusion */
        if (occlusion) {
            code = runmap_andnot(mem, dl->mask, occlusion, 3, &dl->mask);
            if (code < 0) return code;
            /* Although we should logically recompute the perimeter here,
               experience indicates that it helps to simply keep the old
               perimeter because occlusion by background objects often hint
               that this object should go into the background as well. */
        }
        /* Determine occlusions by foreground objects */
        if ( (code = runmap_and(mem, dl->mask, fgmapi, 0, &occlusion)) < 0 ||
             (code = runmap_and(mem, dl->mask, fgmap, 0, &tmp)) < 0        ||
             (code = runmap_or(mem, tmp, occlusion, 3, &occlusion)) < 0     )
            return code;
        /* Applies foreground occlusion */
        if (occlusion) {
            /* Compute clipped mask */
            initial_perim = dl->perimeter;
            runmap_info(occlusion, 0, &occlusion_perim);
            code = runmap_andnot(mem, dl->mask, occlusion, 2, &clipped);
            if (code < 0) return code;
            /* Compute perimeter after occlusion */
            runmap_info(clipped, &dl->area, &dl->perimeter);
            clipped_perim = dl->perimeter;
            /* Compute length of background boundary */
            bg_perim = (initial_perim - occlusion_perim + clipped_perim ) / 2;
        } else {
            /* Compute perimeters without occlusions */
            occlusion_perim = 0;
            initial_perim = clipped_perim = bg_perim = dl->perimeter;
        }
        /* Object marked as text or already marked as foreground */
        if (dl->flags & (DLIST_TEXT|DLIST_FOREGROUND)) {
            dl->flags |= DLIST_FOREGROUND;
            goto cleanup;               
        }
        /* It is true that internal boundaries caused by occlusions do not
           contribute to the cost of coding the object in the background
           (e.g. text on a solid background).  The same cannot be said from
           external boundaries caused by occlusions.  Consider for instance an
           object delineated by a thin foreground border: all its boundaries
           are caused by the border occlusion, and yet it is untrue that such
           an object can be coded in the background at zero cost.  This cost
           critically depends on the border width, but this is too difficult
           to measure.  The following adjustement attempts to count such
           external boundaries for a fraction of the cost of a regular
           background boundary.  This is not perfect.  Yet counting a non zero
           value gives the user a chance to correct possible problems by
           adjusting the gobal threshold. */
        bg_perim += min(initial_perim, clipped_perim) / 2;
        /* Determine color distance with related bg candidates */
        code = scan_bg_candidates(mem, dl, &dist);
        if (code < 0) return code;
#ifdef DEBUG
        if (debug)
            fprintf(STDOUT,"---- %3d * %6d ? %6d\n", 
                    dist, bg_perim, clipped_perim);
#endif        
        /* Renormalize in order to avoid integer overflows */
        while (bg_perim >= 0x10000 || clipped_perim >= 0x10000) {
            bg_perim >>= 1;
            clipped_perim >>= 1;
        }
        /* Choose option with minimum estimated description length */
        if (bg_perim * dist * 100 < clipped_perim * nthreshold * 256) {
            dl->flags |= DLIST_BACKGROUND;
        } else { 
            dl->flags |= DLIST_FOREGROUND;
        }
        /* --- Foreground/background decision has now been made --- */
    cleanup:
#ifdef DEBUG
	if (debug) drawlist_print(dl, cdev);
#endif
	/* Add element to the proper mask/flags */
	if (dl->flags & DLIST_FOREGROUND) {
            /* Replace foreground masks by clipped masks */
            if (occlusion) {
                runmap_free(dl->mask);
                dl->mask = clipped;
            }
            /* Check for invisible objects */
            if (! dl->mask)
                dl->flags = 0;
            /* Update global foreground flags and mask */
	    fgflags |= dl->flags;
	    code = runmap_or(mem, fgmapi, dl->mask, 1, &fgmapi);
            if (code >= 0 && runmap_becoming_big(fgmapi)) {
                code = runmap_or(mem, fgmap, fgmapi, 3, &fgmap);
                fgmapi = 0;
            }
	} else {
            /* Get rid of clipped mask (area and perimeter still refer to it) */
            runmap_free(clipped);
            /* Check for invisible objects */
            if (! dl->mask)
                dl->flags = 0; 
            /* Update global background flags and mask */
            if (dl->color != 0xffffff)
                bgflags |= dl->flags;
	    code = runmap_or(mem, bgmapi, dl->mask, 1, &bgmapi);
            if (code >= 0 && runmap_becoming_big(bgmapi)) {
                code = runmap_or(mem, bgmap, bgmapi, 3, &bgmap);
                bgmapi = 0;
            }
        }
        if (code < 0) return code;
	/* Remove invisible elements */
	if (! dl->flags)
            drawlist_remove(mem, &cdev->head, dl);
        /* Continue */
        dl = dlprev;
    }
    /* Gather global bitmaps */
    code = runmap_or(mem, fgmap, fgmapi, 3, &fgmap);
    if (code < 0) return code;
    runmap_free(bgmapi);
    runmap_free(bgmap);
    /* Remove white foreground objects over empty background.
       These are invisible and yet might prevent bitonal encoding. */
    if (! bgflags) {
        drawlist *dlprev = cdev->head.last;
        while ((dl = dlprev)) {
            dlprev = dl->prev;
	    if ((dl->color == 0xffffff) && (dl->flags & DLIST_FOREGROUND)) {
		code = runmap_andnot(mem, fgmap, dl->mask, 1, &fgmap);
                drawlist_remove(mem, &cdev->head, dl);
		if (code < 0) return code;
	    } 
        }
    }
    /* Return */
    *fgmapout = fgmap;
    *fgflagsout = fgflags;
    *bgflagsout = bgflags;
    return 0; 
}

/* Prepare string describing fg and bg content */
private char *
describe_fg_bg(uint fgflags, uint bgflags)
{
    static char out[40];
    /* Foreground */
    if (!fgflags)
        strcpy(out," fg-empty");
    else if (fgflags & DLIST_HASCOLOR)
        strcpy(out," fg-color");
    else if (fgflags & DLIST_HASGRAY)
        strcpy(out," fg-gray");	    
    else
        strcpy(out," fg-bw");	    	    
    /* Background */
    if (!bgflags)
        strcat(out," bg-empty");
    else if (bgflags & DLIST_HASCOLOR)
        strcat(out," bg-color");
    else
        strcat(out," bg-gray");	  
    /* Background composed of flat colors only */
    if (bgflags && !(bgflags & DLIST_INPURE))
        strcat(out," bg-lineart");
    /* Return in static buffer */
    return out;
}

/* Calls everything in the correct order */
private int
process_drawlist(gx_device_djvu *cdev,
                 runmap **fgmapout, uint *fgflagsout, uint *bgflagsout, 
                 const char **commentout )
{
    int code;
    runmap *fgmap = 0;
    uint fgflags;
    uint bgflags;
#ifdef DEBUG
    if (gs_debug_c(DEBUG_CHAR_P2MEM)) 
        p2mem_diag(cdev->pmem, cdev);
#endif
    /* Break potential lowcolor images */
    code = lowcolor_process(cdev);
    if (code < 0) goto xit;
    /* Indentify potential background candidates */
    code = identify_bg_candidates(cdev);
    if (code < 0) goto xit;
    /* Classify drawlist components */
    code = classify_drawlist(cdev->pmem, cdev, &fgmap, &fgflags, &bgflags);
    if (code < 0) goto xit;
#ifdef DEBUG
    if (gs_debug_c(DEBUG_CHAR_P2MEM)) 
        p2mem_diag(cdev->pmem, cdev);
#endif
    /* Terminate */
    if (fgflagsout)
        *fgflagsout = fgflags;
    if (bgflagsout)
        *bgflagsout = bgflags;
    if (commentout)
        *commentout = describe_fg_bg(fgflags, bgflags);
    if (fgmapout) {
        *fgmapout = fgmap;
        fgmap = 0;
    }
 xit:
    runmap_free(fgmap);
    return code;
}





/* ======================================
       D J V U M A S K    D E V I C E
   --------------------------------------
    Output foreback mask as a PBM image.
   ====================================== */

/* Device template. */
private djvu_proc_process(djvumask_process);
gx_device_djvu gs_djvumask_device = {
    djvu_device_body("djvumask", djvu_device_procs, djvumask_process)
};

/* Device processing function. */
private int
djvumask_process(gx_device_djvu *cdev)
{
    int code;
    runmap *fgmap;
    uint fgflags;
    uint bgflags;
    const char *comment;
    /* Classify drawlist components */
    code = process_drawlist(cdev, &fgmap, &fgflags, &bgflags, &comment);
    if (code < 0) return code;
    /* Save mask as PBM */
    if (cdev->outputfile) {
        code = runmap_save(fgmap, cdev->outputfile, 
                           0, cdev->width-1, 0, cdev->height-1, comment);
        if (code < 0) return code;
    }
    /* Print message */
    if (! cdev->quiet) {
        fprintf(STDOUT, "Page %dx%d (%s )\n", 
                cdev->width, cdev->height, comment);
        fflush(STDOUT);
    }
    /* Terminate */
    runmap_free(fgmap);
#ifdef DEBUG
    if (gs_debug_c(DEBUG_CHAR_P2MEM))
        p2mem_diag(cdev->pmem, cdev);
#endif
    return 0;
}



/* ======================================
    F O R E G R O U N D   P A L E T T E
   ====================================== */


/* Compute foreground palette */
private int
quantize_fg_colors(gx_device_djvu *cdev)
{
    /* Variables */
    int code;
    p2mem *mem = cdev->pmem;
    htable *colorhash = 0;
    colordata *cd;
    drawlist *dl;
    /* Compute color histogram */
    if (! (colorhash = colorhash_alloc(mem))) {
        code = gs_note_error(gs_error_VMerror);
        goto xit;
    }
    for (dl=cdev->head.last; dl; dl=dl->prev) 
	if (dl->flags & DLIST_FOREGROUND) 
            if ((code = colorhash_add(colorhash, dl->color, dl->area)) < 0)
                goto xit;
    /* Make a fake palette if foreground is empty */
    if (colorhash->nelems == 0) {
        if (! (cdev->fgpalette = p2mem_alloc(mem, sizeof(gx_color_index)))) {
            code = gs_note_error(gs_error_VMerror);
            goto xit;
        }
        cdev->fgpalette[0] = 0x000000;
        cdev->fgpalettesize = 1;
        code = 0;
        goto xit;
    }
    /* Call color quantizer */
    code = color_quantize(mem, colorhash, cdev->fgcolors,
                          &cdev->fgpalette, &cdev->fgpalettesize);
    if (code < 0)
        goto xit;
#ifdef DEBUG
    if (gs_debug_c(DEBUG_CHAR_COLOR))
	fprintf(STDOUT,"Found %d colors, reduced to %d colors\n",
		colorhash->nelems, cdev->fgpalettesize);
#endif
    /* Fill palcolor in foreground drawlist elements */
    for (dl=cdev->head.last; dl; dl=dl->prev) {
	if (dl->flags & DLIST_FOREGROUND) {
	    gx_color_index color = dl->color;
	    cd = htable_lookup(colorhash, &color, hash_rgb(color), false);
	    dl->palcolor = cd->w;
	}
    }
    /* Success */
    code = 0;
 xit:
    /* Free everything and return */
    colorhash_free(colorhash);
    return code;
}



/* ======================================
     F O R E G R O U N D   O U T P U T
   ====================================== */

/* Output foreground runs as a CRLE image */
private int
save_foreground(gx_device_djvu *cdev, const char *comment)
{
    int code;
    int fgnum = 0;
    runmap **rmaps = 0;
    gx_color_index *colors = 0;
    drawlist *dl;
    /* Count foreground components */
    for (dl=cdev->head.last; dl; dl=dl->prev)
	if ((dl->flags & DLIST_FOREGROUND) && dl->mask) 
            fgnum += 1;
    /* Allocate rmap and color arrays */
    rmaps = p2mem_alloc(cdev->pmem, sizeof(runmap*) * fgnum);
    colors = p2mem_alloc(cdev->pmem, sizeof(gx_color_index) * fgnum);
    if (fgnum && !(rmaps && colors)) {
        code = gs_note_error(gs_error_VMerror);
        goto xit;
    }        
    /* Fill rmap and color arrays */    
    fgnum = 0;
    for (dl=cdev->head.last; dl; dl=dl->prev)
        if ((dl->flags & DLIST_FOREGROUND) && dl->mask) {
            rmaps[fgnum] = dl->mask;
            colors[fgnum] = dl->palcolor;
            fgnum += 1;
        }
    /* Save */
    code = crle_save(cdev->pmem,
                     cdev->width, cdev->height,
                     cdev->fgpalette, cdev->fgpalettesize,
                     fgnum, rmaps, colors, comment, 
                     cdev->outputfile);
    /* Done */
 xit:
    p2mem_free(cdev->pmem, rmaps);
    p2mem_free(cdev->pmem, colors);
    return code;
}





/* ======================================
     B A C K G R O U N D    O U T P U T
   ====================================== */

/* Save background without subsampling */
private int
save_background_no_subsampling(gx_device_djvu *cdev)
{
    FILE *f = cdev->outputfile;
    int sraster = cdev->width * 3;
    int bandh = cdev->height;
    byte *band = 0;
    int bandsize;
    int nbands = 1;
    int bandy;
    /* Compute bandsize and allocate band buffer */
    for(;;) {
        bandh = (cdev->height + nbands - 1) / nbands;
        bandsize = sraster * bandh;
        if (bandh == 1 || bandsize < cdev->maxbitmap) 
            if ((band = p2mem_alloc(cdev->pmem, bandsize)))
                break;
        if (bandh == 1)
            return_VMerror;
        nbands += 1;
    }
    /* Write PPM header */
    fprintf(f, "P6\n%d %d\n255\n", cdev->width, cdev->height);
    /* Write data */
    for(bandy=0; bandy < cdev->height; bandy += bandh) {
#ifdef DEBUG
        if (gs_debug_c(DEBUG_CHAR_BG))
            fprintf(STDOUT,"band %d .. %d\n",bandy, bandy+bandh);
#endif
        bandh = min(bandh, cdev->height - bandy);
        memset(band, 0xff, bandsize);
        drawlist_play(&cdev->head, DLIST_BACKGROUND,
                      band, sraster, 0, bandy, cdev->width, bandh);
        if (fwrite(band, sraster, bandh, f) < bandh) {
            p2mem_free(cdev->pmem, band);
            return_error(gs_error_ioerror);
        }
    }
    /* End */
    p2mem_free(cdev->pmem, band);
    return 0;
}

/* Internal: in-place masked subsampling on a single row */
private void 
masksub(byte *rgb, int rgbraster, const byte *mask, int maskraster, 
        int subsample, int subw)
{
    byte *out = rgb;
    static int inv[12*12+1];
    int i,j;
    /* Create inversion table */
    if (! inv[1])
        for (i=1; i<=12*12; i++)
            inv[i] = (0x10000 + (i-1)/2) / i;
    ASSERT(subsample>1 && subsample<=12);
    /* Subsample */
    while (subw > 0) {
        int r, g, b, c;
        byte *rgbrow = rgb;
        const byte *maskrow = mask;
        r = g = b = c = 0;
        for (i=0; i<subsample; i++) {
            for (j=0; j<subsample; j++)
                if (! maskrow[j]) {
                    byte *pix = rgbrow + j + j + j;
                    r += pix[0];
                    g += pix[1];
                    b += pix[2]; 
                    c += 1;
                }
            rgbrow += rgbraster;
            maskrow += maskraster;
        }
        if (c == 0) {
            /* Pixel is entirely masked */
            rgbrow = rgb;
            for (i=0; i<subsample; i++) {
                for (j=0; j<subsample; j++) {
                    byte *pix = rgbrow + j + j + j;
                    r += pix[0];
                    g += pix[1];
                    b += pix[2];
                    c += 1;
                }
                rgbrow += rgbraster;
            }
        }
        /* Store result */
        out[0] = (r * inv[c]) >> 16;
        out[1] = (g * inv[c]) >> 16;
        out[2] = (b * inv[c]) >> 16;
        /* Next pixel */
        out += 3;
        rgb += 3 * subsample;
        mask += subsample;
        subw -= 1;
    }
}

/* Save background with masked subsampling */
private int
save_background(gx_device_djvu *cdev, runmap *mask, int subsample)
{
    FILE *f = cdev->outputfile;
    int subw = (cdev->width + subsample - 1) / subsample;
    int subh = (cdev->height + subsample - 1) / subsample;
    int deadh = subh * subsample - cdev->height;
    int bmapraster = subw * subsample;
    int bmapsize = bmapraster * subsample;
    byte *bmap = 0;
    int bandraster = bmapraster * 3;
    int bandsize;
    int bandh;
    byte *band = 0;
    int nbands = 1;
    int bandy;
    /* Handle full resolution case */
    if (subsample == 1)
        return save_background_no_subsampling(cdev);
    /* Allocate bmap */
    if (! (bmap = p2mem_alloc(cdev->pmem, bmapsize)))
        return_VMerror;
    /* Compute bandsize and allocate band buffer */
    for(;;) {
        bandh = (int)((subh + nbands - 1) / nbands) * subsample;
        bandsize = bandraster * bandh;
        if (bandh <= subsample || bandsize < cdev->maxbitmap) 
            if ((band = p2mem_alloc(cdev->pmem, bandsize)))
                break;
        if (bandh <= subsample) {
            p2mem_free(cdev->pmem, bmap);
            return_VMerror;
        }
        nbands += 1;
    }
    /* Write PPM header */
    fprintf(f, "P6\n%d %d\n255\n", subw, subh);
    /* Iterate on bands */
    for(bandy = -deadh; bandy < cdev->height; bandy += bandh) {
        int ry;
#ifdef DEBUG
        if (gs_debug_c(DEBUG_CHAR_BG))
            fprintf(STDOUT,"band %d .. %d\n", bandy, bandy+bandh);
#endif
        /* Render band */
        bandh = min(bandh, cdev->height - bandy);
        memset(band, 0xff, bandsize);
        drawlist_play(&cdev->head, DLIST_BACKGROUND,
                      band, bandraster, 0, bandy, cdev->width, bandh);
        /* Iterate on subsampled rows */
        for (ry = bandy; ry < bandy + bandh; ry += subsample) {
            byte *bandrow = band + (ry - bandy) * bandraster;
            byte *bmaprow = bmap; 
            uint y;
            /* Render bmap */
            memset(bmap, 1, bmapsize);
            for (y = max(0,ry); 
                 y < (uint)(ry+subsample); 
                 y++, bmaprow += bmapraster) {
                memset(bmaprow, 0, cdev->width);
                if (mask && y>=mask->ymin && y<=mask->ymax) {
                    byte *data = mask->data + mask->rows[y - mask->ymin];
                    byte *dataend = mask->data + mask->rows[y - mask->ymin + 1];
                    uint x = 0;
                    while (data < dataend) {
                        rl_decode(data, x);
                        if (data < dataend) {
                            uint x1 = x;
                            rl_decode(data, x);
                            if (x > x1 && x <= bmapraster)
                                memset(bmaprow + x1, 1, x - x1);
                        }
                    }
                }
            }
            /* Subsample row */
            masksub(bandrow, bandraster, bmap, bmapraster, subsample, subw);
            /* Save subsampled row */
            if (fwrite(bandrow, 3, subw, f) < subw) {
                p2mem_free(cdev->pmem, band);
                p2mem_free(cdev->pmem, bmap);
                return_error(gs_error_ioerror);
            }
        }
    }
    /* End */
    p2mem_free(cdev->pmem, band);
    p2mem_free(cdev->pmem, bmap);
    return 0;
}




/* ======================================
       D J V U S E P    D E V I C E
   --------------------------------------
   Output separation file:
   - foreground as a CRLE image file
   - background (possibly) as a PPM image file
   - optional comment lines
   ====================================== */

/* Device template. */
private djvu_proc_process(djvusep_process);
gx_device_djvu gs_djvusep_device = {
    djvu_device_body("djvusep", djvu_device_procs, djvusep_process)
};

/* Device processing function. */
private int
djvusep_process(gx_device_djvu *cdev)
{
    int code;
    runmap *fgmap;
    uint fgflags;
    uint bgflags;
    const char *comment;
    /* Classify drawlist components */
    code = process_drawlist(cdev, &fgmap, &fgflags, &bgflags, &comment);
    if (code < 0) return code;
    /* Perform color quantization on foreground */
    code = quantize_fg_colors(cdev);
    if (code < 0) return code;
#ifdef DEBUG
    if (gs_debug_c(DEBUG_CHAR_P2MEM))
        p2mem_diag(cdev->pmem, cdev);
#endif
    /* Save separation */
    if (cdev->outputfile) {
        /* Generate RLE encoded foreground */
        code = save_foreground(cdev, comment);
        /* Generate PPM background */
        if (code >= 0 && bgflags) {
            int bgsubsample = cdev->bgsubsample;
            if (cdev->autohires && !fgflags)
                bgsubsample = 1;
            code = save_background(cdev, fgmap, bgsubsample);
        }
        /* Generate comments with mark information */
        if (code >=0 ) {
            drawlist *dl;
            pdfmark *mark = cdev->marks;
            cdev->marks = 0;
            while (mark) { /* reverse list */
                pdfmark *next = mark->next;
                mark->next = cdev->marks;
                cdev->marks = mark;
                mark = next;
            }
            for (mark=cdev->marks; mark; mark=mark->next) {
                if (mark->type == 'L') {
                    struct pdfmark_l_s *m = &mark->u.l;
                    fprintf( cdev->outputfile, 
                             "# L %dx%d%+d%+d %s\n", 
                             m->w, m->h, m->x, m->y, m->uri );
                } else if (mark->type == 'B') {
                    struct pdfmark_b_s *m = &mark->u.b;
                    fprintf( cdev->outputfile,
                             "# B %d %s (#%d)\n",
                             m->count, m->title, m->page );
                } else if (mark->type == 'P') {
		    fprintf( cdev->outputfile,
                             "# P %s\n", mark->u.p.title);
                }
            }
            for (dl=cdev->head.first; dl; dl=dl->next) {
		runmap *mask = 0;
		txtmark *mark = 0;
                if ((dl->flags & DLIST_TEXT) && 
                    (mask = dl->mask) &&  (mark = dl->text) )
                    fprintf( cdev->outputfile, 
                             "# T %d:%d %d:%d %dx%d%+d%+d %s\n", 
                             mark->x, mark->y, mark->w, mark->h,
                             mask->xmax - mask->xmin +1, 
                             mask->ymax - mask->ymin +1, 
                             mask->xmin, mask->ymin, mark->arg );
	    }
        }
        /* Check for errors */
        if (ferror(cdev->outputfile)) 
            code = gs_error_ioerror;
        if (code < 0) 
            return code; 
    }
    /* Print message */
    if (! cdev->quiet) {
        fprintf(STDOUT,"Page %dx%d (%s )\n",
                cdev->width, cdev->height, comment);
        fflush(STDOUT);
    }
    /* Terminate */
    runmap_free(fgmap);
#ifdef DEBUG
    if (gs_debug_c(DEBUG_CHAR_P2MEM))
        p2mem_diag(cdev->pmem, cdev);
#endif
    return 0;
}



/* ======================================
         T H E   E N D 
   ====================================== */
