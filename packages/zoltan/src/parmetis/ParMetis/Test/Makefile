include ../Makefile.in

BINDIR = ../Graphs

INCLUDES = -I../ParMETISLib $(INCDIR) 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)


LIBSDIR = -L.. $(LIBDIR) 
LIBS = -lparmetis -lmetis $(XTRALIBS) -lm
PARMETISLIB = ../libparmetis.a ../libmetis.a


PTESTOBJS = main.o ptest.o tio.o tstadapt.o
          
          
.c.o:
	$(CC) $(CFLAGS) -c $*.c


default: $(BINDIR)/ptest

$(BINDIR)/ptest: $(PTESTOBJS) $(PARMETISLIB)
	$(LD)  -o $@ $(PTESTOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@


clean:
	rm -f *.o ;\
	rm -f $(BINDIR)/ptest$(VERNUM) 
