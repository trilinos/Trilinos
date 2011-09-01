# Simple Makefile for TKLU.
# Need to set TRILINOS_HOME correctly for amesos includes.

default: all

ccode: all

# Need to set this !!! TODO
TRILINOS_HOME = /home/srajama/proj/Trilinos_Apr14
INCLUDE_HOME = /home/srajama/proj/Trilinos_Apr14/hash_build/include

SS_HOME = $(TRILINOS_HOME)/packages/amesos/src/SuiteSparse
TEUCHOS_HOME = $(TRILINOS_HOME)/packages/teuchos/src/

#include $(SS_HOME)/UFconfig/amesos_UFconfig.mk

# for testing only:
# TEST = -DTESTING

CC = g++
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


INC = ../Include/klu2.h ../Include/klu2_internal.h ../Include/klu2_version.h \
    $(SS_HOME)/UFconfig/amesos_UFconfig.h Makefile

I = -I../../AMD/Include -I../../COLAMD/Include -I../../BTF/Include \
    -I../Include -I$(SS_HOME)/UFconfig -I$(SS_HOME)/AMD/Include \
    -I$(SS_HOME)/BTF/Include -I$(SS_HOME)/COLAMD/Include -I$(INCLUDE_HOME)

all: example

$(OBJ): $(INC)

#-------------------------------------------------------------------------------

example: klu2_simple.cpp libamesos.a
	$(C) $(I) klu2_simple.cpp libamesos.a -o klu2_simple.exe

allcpp:
	$(C) $(I) *.cpp

purge: distclean

distclean: clean
	- $(RM) klu2_simple.exe

clean:
	- $(RM) $(CLEAN)

#-------------------------------------------------------------------------------

