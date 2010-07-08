# Simple Makefile for TKLU.
# Need to set TRILINOS_HOME correctly for amesos includes.

default: all

ccode: all

TRILINOS_HOME = /home/srajama/proj/Trilinos_Apr14
SS_HOME = $(TRILINOS_HOME)/packages/amesos/src/SuiteSparse

#include $(SS_HOME)/UFconfig/amesos_UFconfig.mk

# for testing only:
# TEST = -DTESTING

CC = cc
CPLUSPLUS = g++
CFLAGS = -O3
RANLIB = ranlib
AR = ar cr

# delete and rename a file
RM = rm -f
CLEAN = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d *.gcda *.gcno
MV = mv -f

# Fortran compiler (not normally required)
F77 = f77
F77FLAGS = -O
F77LIB =

# C and Fortran libraries
LIB = -lm

# For compiling MATLAB mexFunctions (MATLAB 7.5 or later)
MEX = mex -O -largeArrayDims -lmwlapack -lmwblas

C = $(CPLUSPLUS) $(CFLAGS)


INC = ../Include/tklu.h ../Include/tklu_internal.h ../Include/tklu_version.h \
    $(SS_HOME)/UFconfig/amesos_UFconfig.h Makefile

I = -I../../AMD/Include -I../../COLAMD/Include -I../../BTF/Include \
    -I../Include -I$(SS_HOME)/UFconfig -I$(SS_HOME)/AMD/Include \
    -I$(SS_HOME)/BTF/Include -I$(SS_HOME)/COLAMD/Include

all: library

library: libtklu.a

KLU_D = tklu_d.o tklu_d_kernel.o tklu_d_dump.o \
    tklu_d_factor.o tklu_d_free_numeric.o tklu_d_solve.o \
    tklu_d_scale.o tklu_d_refactor.o \
    tklu_d_tsolve.o tklu_d_diagnostics.o tklu_d_sort.o tklu_d_extract.o

KLU_Z = tklu_z.o tklu_z_kernel.o tklu_z_dump.o \
    tklu_z_factor.o tklu_z_free_numeric.o tklu_z_solve.o \
    tklu_z_scale.o tklu_z_refactor.o \
    tklu_z_tsolve.o tklu_z_diagnostics.o tklu_z_sort.o tklu_z_extract.o

COMMON = \
    tklu_free_symbolic.o tklu_defaults.o tklu_analyze_given.o \
    tklu_analyze.o tklu_memory.o 

OBJ = $(COMMON) $(KLU_D) $(KLU_Z)

libtklu.a: $(OBJ)
	$(AR) libtklu.a $(OBJ)
	$(RANLIB) libtklu.a

$(OBJ): $(INC)

#-------------------------------------------------------------------------------

tklu_d.o: ../Source/tklu.cpp
	$(C) -c $(I) $< -o $@

tklu_z.o: ../Source/tklu.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_kernel.o: ../Source/tklu_kernel.cpp
	$(C) -c $(I) $< -o $@

tklu_z_kernel.o: ../Source/tklu_kernel.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_sort.o: ../Source/tklu_sort.cpp
	$(C) -c $(I) $< -o $@

tklu_z_sort.o: ../Source/tklu_sort.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_diagnostics.o: ../Source/tklu_diagnostics.cpp
	$(C) -c $(I) $< -o $@

tklu_z_diagnostics.o: ../Source/tklu_diagnostics.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_dump.o: ../Source/tklu_dump.cpp
	$(C) -c $(I) $< -o $@

tklu_z_dump.o: ../Source/tklu_dump.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_factor.o: ../Source/tklu_factor.cpp
	$(C) -c $(I) $< -o $@

tklu_z_factor.o: ../Source/tklu_factor.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_free_numeric.o: ../Source/tklu_free_numeric.cpp
	$(C) -c $(I) $< -o $@

tklu_z_free_numeric.o: ../Source/tklu_free_numeric.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_extract.o: ../Source/tklu_extract.cpp
	$(C) -c $(I) $< -o $@

tklu_z_extract.o: ../Source/tklu_extract.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_refactor.o: ../Source/tklu_refactor.cpp
	$(C) -c $(I) $< -o $@

tklu_z_refactor.o: ../Source/tklu_refactor.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_scale.o: ../Source/tklu_scale.cpp
	$(C) -c $(I) $< -o $@

tklu_z_scale.o: ../Source/tklu_scale.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_solve.o: ../Source/tklu_solve.cpp
	$(C) -c $(I) $< -o $@

tklu_z_solve.o: ../Source/tklu_solve.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

tklu_d_tsolve.o: ../Source/tklu_tsolve.cpp
	$(C) -c $(I) $< -o $@

tklu_z_tsolve.o: ../Source/tklu_tsolve.cpp
	$(C) -c -DCOMPLEX $(I) $< -o $@

#-------------------------------------------------------------------------------

tklu_analyze.o: ../Source/tklu_analyze.cpp
	$(C) -c $(I) $< -o $@

tklu_analyze_given.o: ../Source/tklu_analyze_given.cpp
	$(C) -c $(I) $< -o $@

tklu_defaults.o: ../Source/tklu_defaults.cpp
	$(C) -c $(I) $< -o $@

tklu_free_symbolic.o: ../Source/tklu_free_symbolic.cpp
	$(C) -c $(I) $< -o $@

tklu_memory.o: ../Source/tklu_memory.cpp
	$(C) -c $(I) $< -o $@

#-------------------------------------------------------------------------------

purge: distclean

distclean: clean
	- $(RM) libtklu.a 

clean:
	- $(RM) $(CLEAN)

#-------------------------------------------------------------------------------

