include ../make.inc

LUSRC   = ../SRC
HEADER	= -I$(LUSRC) -I$(MATLAB)/extern/include 
OFLAGS  = -O

SOLVEOBJS 	= mexlusolve.o
FACTOROBJS	= mexsuperlu.o

all:            mexlusolve mexsuperlu

mexlusolve:	$(SOLVEOBJS) ../$(SUPERLULIB)
		${MATLAB}/bin/cmex -lm ${OFLAGS} ${SOLVEOBJS} \
		../$(SUPERLULIB) $(BLASLIB)

mexsuperlu:	$(FACTOROBJS) ../$(SUPERLULIB)
		${MATLAB}/bin/cmex -lm ${OFLAGS} ${FACTOROBJS} \
		../$(SUPERLULIB) $(BLASLIB)

.c.o:
	$(CC) $(CFLAGS) $(HEADER) -c $< $(VERBOSE)

clean:	
	rm -f *.o *.mex*
