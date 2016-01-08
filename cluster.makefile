CC     = cc
OFILES = cluster.o
LFILES = bioplib/fgetsany.o bioplib/array2.o bioplib/upstrncmp.o bioplib/GetWord.o bioplib/OpenStdFiles.o

all : cluster


cluster :  $(OFILES) $(LFILES)
	$(CC) -o cluster $(OFILES) $(LFILES) -lm

.c.o  :
	$(CC) -o $@ -c $<

clean :
	\rm -f $(OFILES) $(LFILES)
