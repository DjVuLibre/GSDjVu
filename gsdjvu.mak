### ---------------- DjVu files ---------------- ###

djvu_=$(GLOBJ)gdevdjvu.$(OBJ)

$(DD)djvumask.dev : $(djvu_)
	$(SETDEV) $(DD)djvumask $(djvu_)

$(DD)djvusep.dev : $(djvu_)
	$(SETDEV) $(DD)djvusep $(djvu_)

$(GLOBJ)gdevdjvu.$(OBJ) : $(GLSRC)gdevdjvu.c
	$(GLCC) $(GLO_)gdevdjvu.$(OBJ) $(C_) $(GLSRC)gdevdjvu.c

