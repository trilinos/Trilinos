#######################################################################
#
#  This makefile creates the example programs for the linear equation
#  routines in SuperLU_Dist.
#
#  The command
#       make
#  without any arguments creates all the example programs.
#  The command
# 	make double
#  creates double precision real example programs.
#  The command
#       make complex16
#  creates double precision complex example programs.
#
#  The executable files are called
#       double real:	pddrive pddrive1 pddrive2 pddrive3 pddrive4
#	double complex: pzdrive pzdrive1 pzdrive2 pzdrive3 pzdrive4 
#
#  Alternatively, you can create example programs individually by
#  typing the command (for example)
#	make pddrive
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#
#######################################################################
include ../make.inc
INCLUDEDIR = -I../SRC

DEXM	= pddrive.o
DEXM1	= pddrive1.o
DEXM2	= pddrive2.o
DEXM3	= pddrive3.o
DEXM4	= pddrive4.o
ZEXM	= pzdrive.o
ZEXM1	= pzdrive1.o
ZEXM2	= pzdrive2.o
ZEXM3	= pzdrive3.o
ZEXM4	= pzdrive4.o

all: pddrive pddrive1
#all: double complex16
#double: 	pddrive pddrive1 pddrive2 pddrive3 pddrive4
#complex16: 	pzdrive pzdrive1 pzdrive2 pzdrive3 pzdrive4

pddrive: $(DEXM) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DEXM) $(LIBS) -lm -o $@

pddrive1: $(DEXM1) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DEXM1) $(LIBS) -lm -o $@

pddrive2: $(DEXM2) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DEXM2) $(LIBS) -lm -o $@

pddrive3: $(DEXM3) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DEXM3) $(LIBS) -lm -o $@

pddrive4: $(DEXM4) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DEXM3) $(LIBS) -lm -o $@

pzdrive: $(ZEXM) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZEXM) $(LIBS) -lm -o $@

pzdrive1: $(ZEXM1) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZEXM1) $(LIBS) -lm -o $@

pzdrive2: $(ZEXM2) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZEXM2) $(LIBS) -lm -o $@

pzdrive3: $(ZEXM3) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZEXM3) $(LIBS) -lm -o $@

pzdrive4: $(ZEXM4) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZEXM4) $(LIBS) -lm -o $@

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o p[dz]drive p[dz]drive[1-9]

