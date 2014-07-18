### ---------------- DjVu files ---------------- ###

djvu_=$(@@djvu@@OBJ)gdevdjvu.$(OBJ)

$(DD)djvumask.dev : $(djvu_)
	$(SETDEV) $(DD)djvumask $(djvu_)

$(DD)djvusep.dev : $(djvu_)
	$(SETDEV) $(DD)djvusep $(djvu_)

$(@@djvu@@OBJ)gdevdjvu.$(OBJ) : $(@@djvu@@SRC)gdevdjvu.c $(GLGEN)arch.h
	$(@@djvu@@CC) $(@@djvu@@O_)gdevdjvu.$(OBJ) $(C_) $(@@djvu@@SRC)gdevdjvu.c \
	    -DGS_VERSION=$(GS_VERSION)

