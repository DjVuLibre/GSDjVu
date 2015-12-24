# initializations
function _esc_init() {
    _ctrl = ""
    for (i=1; i<32; i++) { _ctrl = _ctrl sprintf("%c",i) }
    _esc = "[\"\\&" _ctrl "]"
    _ctrl = "[" _ctrl "]"
}
function _ord_init() {
    for (i=0; i<256; i++) {
	_ord[sprintf("%c",i)]=i
    }
}
function _amp_init() {
    _amp["&amp;"] = "&"
    _amp["&lt;"] = "<"
    _amp["&gt;"] = ">"
    _amp["&apos;"] = "'"
    _amp["&quot;"] = "\\\""
}
BEGIN {
    # inits
    _esc_init()
    _ord_init()
    _amp_init()
    # xml parsing trick
    RS=">"
    # variables
    delete meta
    pheight=0
    pwidth=0
    location=""
    context=""
    content=""
    pageno=0
    dpi = 300 # use awk -f xml2dsed.awk dpi=xxx to override
    dometa=1  # use awk -f xml2dsed.awk dometa=0 to override
    dotext=1  # use awk -f xml2dsed.awk dotext=0 to override
}
# return character code
function ord(str) {
    return _ord[substr(str,1,1)]
}
# print c string
function pstr(str,tmp) {
    printf("\"")
    while (str) {
	tmp = match(str,_esc) # char classes do not always work
	if (tmp == 0) {
	    printf("%s",str)
	    str = ""
	} else if (tmp > 1) {
	    printf("%s",substr(str,1,tmp-1))
	    str = substr(str,tmp)
	} else {
	    tmp = match(str,"^&[a-z]*;")
	    if (tmp != 1) { tmp = "" }
	    if (tmp) { tmp = _amp[substr(str,RSTART,RLENGTH)] }
	    if (tmp) {
		printf("%s",tmp)
		str = substr(str,RLENGTH+1)
	    } else {
		printf("\\%03o", ord(str))
		str = substr(str,2)
	    }
	}
    }
    printf("\"")
}
# sax-like callbacks
function charData(str) {
    if (context == "title") {
	meta["Title"] = str
    } else if (context == "word") {
	gsub(_ctrl," ",str)      # kill control characters
	gsub("\302\240"," ",str) # nbsp in utf-8
	gsub(/  */," ",str)      # simplify spaces
	gsub(/^ /,"",str)        # simplify spaces
	gsub(/ $/,"",str)        # simplify spaces
	if (match(str,/[^ ]/)) { content = str } else { content = "" }
    }
}
function startElement(tag,attrs) {
    if (tag=="head" && ! context && dometa) {
	context = "head"
    } else if (tag=="title" && context=="head") {
	context = "title"
    } else if (tag=="meta" && attrs["name"] && attrs["content"] ) {
	meta[attrs["name"]] = attrs["content"]
    } else if (tag=="body" && ! context && dotext) {
	context = "body"
    } else if (tag=="page" && context == "body") {
	pageno = pageno + 1
	pwidth = attrs["width"] * dpi / 72
	pheight = attrs["height"] * dpi / 72
	context = "page"
	printf("select %d\n", pageno)
	printf("set-txt\n")
	printf("(page %d %d %d %d\n", 0, 0, pwidth, pheight)
    } else if (tag=="word" && context == "page") {
	if (attrs["xMin"] && attrs["xMax"] && attrs["yMin"] && attrs["yMax"]) {
	    context = "word"
	    location = sprintf("%d %d %d %d",
			       attrs["xMin"] * dpi / 72,
			       pheight - attrs["yMin"] * dpi / 72,
			       attrs["xMax"] * dpi / 72,
			       pheight - attrs["yMax"] * dpi / 72 )
	}
    }
}
function endElement(tag) {
    if (tag == "title" && context == "title") {
	context = "head"
    } else if (tag == "head" && context == "head" && length(meta)>0) {
	context = ""
	printf("create-shared-ant\n")
	printf("set-meta\n")
	for (i in meta) {
	    printf("%s\t",i)
	    pstr(meta[i])
	    printf("\n")
	}
	printf(".\n")
    } else if (tag == "page" && context == "page") {
	context = "body"
	printf(")\n.\n")
    } else if (tag == "word" && context == "word" && content) {
	context = "page"
	printf("  (word %s ", location)
	pstr(content)
	printf(")\n")
	content = ""
    }
}
# xml parser (for pdftotext!)
{
    str = $0
    match(str,/^[ \n\r\t\f]*/)
    if (RSTART == 1 && RLENGTH > 0) {
	str = substr(str,RLENGTH+1)
    }
    if (! match(str,/</)) {
	charData(str)
    } else {
	if (RSTART > 1) {
	    arg = substr(str,1,RSTART-1)
	    str = substr(str,RSTART)
	    charData(arg)
	}
	match(str, "^</?[a-zA-Z][a-zA-Z0-9_]*")
	if (RSTART == 1) {
	    tag = substr(str,2,RLENGTH-1)
	    arg = substr(str,RLENGTH+1)
	    # parse attrs
	    delete attrs
	    if (match(arg,"/$")) { arg = substr(arg,1,RSTART-1) }
	    while (arg && match(arg,/[a-zA-Z][a-zA-Z0-9_]*=/)) {
		attr = substr(arg,RSTART,RLENGTH-1)
		arg = substr(arg,RSTART+RLENGTH)
		if (match(arg,/^[ \n\r\f\t]*"/) && match(arg,/"[^"]*"/)) {
		    attrs[attr] = substr(arg,RSTART+1,RLENGTH-2)
		    arg = substr(arg,RSTART+RLENGTH)
		} else if (match(arg,/^[^ \n\r\g\t]*/)) {
		    attrs[attr] = substr(arg,RSTART,RLENGTH)
		    arg = substr(arg,RSTART+RLENGTH)
		}
	    }
	    # callbacks
	    if (match(tag,"^/")) {
		endElement(substr(tag,2))
	    } else {
		startElement(tag,attrs)
		if (match(str,"/$")) { endElement(tag) }
	    }
	} 
    }	    
}
