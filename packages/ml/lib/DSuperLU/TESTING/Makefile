include ../make.inc

#######################################################################
#  This makefile creates the test programs for the linear equation
#  routines in SuperLU.  The test files are grouped as follows:
#
#       ALINTST -- Auxiliary test routines
#       SLINTST -- Single precision real test routines
#       DLINTST -- Double precision real test routines
#       CLINTST -- Double precision complex test routines
#       ZLINTST -- Double precision complex test routines
#
#  Test programs can be generated for all or some of the four different
#  precisions.  Enter make followed by one or more of the data types
#  desired.  Some examples:
#       make single
#       make single double
#  Alternatively, the command
#       make
#  without any arguments creates all four test programs.
#  The executable files are called
#       sp_lintsts
#       sp_lintstd
#       sp_lintstc
#       sp_lintstz
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################

HEADER  = ../SRC

ALINTST = sp_ienv.o

SLINTST = sdrive.o sp_sconvert.o sgssv.o sgssvx.o \
          sp_sget01.o sp_sget02.o sp_sget04.o sp_sget07.o 

DLINTST = ddrive.o sp_dconvert.o dgssv.o dgssvx.o \
          sp_dget01.o sp_dget02.o sp_dget04.o sp_dget07.o 

CLINTST = cdrive.o sp_cconvert.o cgssv.o cgssvx.o \
	  sp_cget01.o sp_cget02.o sp_cget04.o sp_cget07.o

ZLINTST = zdrive.o sp_zconvert.o zgssv.o zgssvx.o \
	  sp_zget01.o sp_zget02.o sp_zget04.o sp_zget07.o

all:    single double complex complex16

single: ./sp_lintsts stest.out

./sp_lintsts: $(SLINTST) $(ALINTST) ../$(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(SLINTST) $(ALINTST) \
        $(TMGLIB) ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

stest.out: sp_lintsts stest.csh
	@echo Testing SINGLE PRECISION linear equation routines 
	csh stest.csh

double: ./sp_lintstd dtest.out

./sp_lintstd: $(DLINTST) $(ALINTST) ../$(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(DLINTST) $(ALINTST) \
        $(TMGLIB) ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

dtest.out: sp_lintstd dtest.csh
	@echo Testing DOUBLE PRECISION linear equation routines 
	csh dtest.csh

complex: ./sp_lintstc ctest.out

./sp_lintstc: $(CLINTST) $(ALINTST) ../$(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(CLINTST) $(ALINTST) \
	$(TMGLIB) ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

ctest.out: sp_lintstc ctest.csh
	@echo Testing SINGLE COMPLEX linear equation routines 
	csh ctest.csh

complex16: ./sp_lintstz ztest.out

./sp_lintstz: $(ZLINTST) $(ALINTST) ../$(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(ZLINTST) $(ALINTST) \
	$(TMGLIB) ../$(SUPERLULIB) $(BLASLIB) -lm -o $@

ztest.out: sp_lintstz ztest.csh
	@echo Testing DOUBLE COMPLEX linear equation routines 
	csh ztest.csh

##################################
# Do not optimize this routine   #
##################################
dlamch.o: dlamch.c ; $(CC) -c $<

timer.o:  timer.c ; $(CC) -O -c $<

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

clean:	
	rm -f *.o sp_lintst* *.out

