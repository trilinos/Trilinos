############################################################################
#
#  Program:         SuperLU_DIST
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   February 4, 1999  version alpha
#
#  Modified:	    September 1, 1999  version 1.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT		= _linux

#
#  The name of the libraries to be created/linked to
#
#DSuperLUroot 	= $(HOME)/ml_dsuper/ml/lib/SuperLU_DIST
DSuperLUroot 	= ../../SuperLU_DIST
DSUPERLULIB   	= $(DSuperLUroot)/superlu$(PLAT).a
#
BLASDEF	     	= -DUSE_VENDOR_BLAS
BLASLIB      	= -lessl
IMPI            = -I/usr/local/mpich/include
#MPILIB		= -L/usr/lpp/ppe.poe/lib -lmpi
#PERFLIB     	= -L/vol1/VAMPIR/lib -lVT
LIBS	     	= $(DSUPERLULIB) $(BLASLIB) $(PERFLIB) $(MPILIB) -lm

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         	= ar
ARCHFLAGS    	= ruv
RANLIB       	= ranlib

CC           	= cc
CFLAGS          = -D_LINUX -O $(IMPI)
NOOPTS		=
FORTRAN         = f77
FFLAGS          = -DLINUX -O 
LOADER	        = cc
LOADOPTS	= 
#
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
CDEFS        = -DAdd_

