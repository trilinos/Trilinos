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
PLAT		= _sp

#
#  The name of the libraries to be created/linked to
#
DSuperLUroot 	= $(HOME)/SuperLU_DIST
DSUPERLULIB   	= $(DSuperLUroot)/superlu$(PLAT).a
#
BLASDEF	     	= -DUSE_VENDOR_BLAS
BLASLIB      	= -lessl
LIBS	     	= $(DSUPERLULIB) $(BLASLIB) $(PERFLIB) $(MPILIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         	= ar
ARCHFLAGS    	= cr
RANLIB       	= ranlib

CC           	= mpcc
CFLAGS          = -D_SP -O3 -qarch=PWR3 -qalias=allptrs -qmaxmem=-1 -DDEBUGlevel=0 -DPRNTlevel=0
NOOPTS		=
FORTRAN         = xlf90
FFLAGS          = -WF,-Dsp -O3 -Q -qfixed -qinit=f90ptr -qarch=pwr3 -qmaxmem=-1
LOADER	        = mpxlf
LOADOPTS	= -bmaxdata:0x80000000
#
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
CDEFS        = -DNoChange

