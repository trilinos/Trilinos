include ../make.inc

#######################################################################
#  This makefile creates the example programs for the linear equation
#  routines in SuperLU.  The files are grouped as follows:
#
#       SLINEXM -- Single precision real example routines
#       DLINEXM -- Double precision real example routines
#       CLINEXM -- Double precision complex example routines
#       ZLINEXM -- Double precision complex example routines
#
#  Example programs can be generated for all or some of the four different
#  precisions.  Enter make followed by one or more of the data types
#  desired.  Some examples:
#       make single
#       make single double
#  Alternatively, the command
#       make
#  without any arguments creates all four example programs.
#  The executable files are called
#       slinsol		slinsolx
#       dlinsol		dlinsolx
#       clinsol		clinsolx
#       zlinsol		zlinsolx
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################

HEADER   = ../SRC

SLINEXM = slinsol.o
SLINXEXM = slinsolx.o

DLINEXM = dlinsol.o
DLINXEXM = dlinsolx.o
SUPERLUEXM = superlu.o sp_ienv.o

CLINEXM = clinsol.o
CLINXEXM = clinsolx.o

ZLINEXM = zlinsol.o
ZLINXEXM = zlinsolx.o


all:    single double complex complex16

single: slinsol slinsolx
double: dlinsol dlinsolx superlu
complex: clinsol clinsolx
complex16: zlinsol zlinsolx

slinsol: $(SLINEXM) ../$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SLINEXM) \
        ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

slinsolx: $(SLINXEXM) ../$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SLINXEXM) \
        ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

dlinsol: $(DLINEXM) ../$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DLINEXM) \
        ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

dlinsolx: $(DLINXEXM) ../$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DLINXEXM) \
        ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

superlu: $(SUPERLUEXM) ../$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SUPERLUEXM) \
        ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

clinsol: $(CLINEXM) ../$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CLINEXM) \
        ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

clinsolx: $(CLINXEXM) ../$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CLINXEXM) \
        ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

zlinsol: $(ZLINEXM) ../$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZLINEXM) \
        ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

zlinsolx: $(ZLINXEXM) ../$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZLINXEXM) \
        ../$(SUPERLULIB) $(BLASLIB) -lm -o $@


F77	= f77
FFLAGS	= -O
F77EXM	= f77_main.o hbcode1.o c_bridge_dgssv.o

f77exm: $(F77EXM) ../$(SUPERLULIB)
	$(F77) $(F77EXM) ../$(SUPERLULIB) $(BLASLIB) -o $@

c_bridge_dgssv.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER) -c $< $(VERBOSE)

.f.o:
	$(F77) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o *linsol *linsolx f77exm superlu
