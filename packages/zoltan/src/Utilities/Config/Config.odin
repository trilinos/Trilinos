##############################################################################
# Zoltan Library for Parallel Applications                                   #
# Copyright (c) 2000,2001,2002, Sandia National Laboratories.                #
# For more info, see the README file in the top-level Zoltan directory.      # 
##############################################################################
##############################################################################
# CVS File Information
#    $RCSfile$
#    $Author$
#    $Date$
#    $Revision$
##############################################################################

##############################################################################
#  Environment variables for compiling the Zoltan and test drivers
#  on the linux cluster qed.sandia.gov.
#  For explanations of these variables, please see Config.generic.
#
#  Note: jobs on this machine should be run using
#           /opt/mpiexec/bin/mpiexec and
#           /opt/torque/bin/qsub
##############################################################################

DEFS 		= 

RANLIB		= ranlib
AR		= ar r

CC          	= mpicc
CPPC       	= mpiCC

INCLUDE_PATH	=
DBG_FLAGS		= -g -O0
OPT_FLAGS		= -O
CFLAGS 			= $(DBG_FLAGS)

F90			    = f90
LOCAL_F90		= f90
F90CFLAGS 		= -DFMANGLE=UNDERSCORE -DNO_MPI2
FFLAGS    		= 
SPPR_HEAD 		= spprinc.most
FARG      		= farg_typical
F90_MODULE_PREFIX 	= -M

#MPI_LIBS			= -lmpich -lpmpich -lnsl -lgcc
#MPI_LIBPATH		= -L/usr/local/mpich-1.2.6-eth/lib

PARMETIS_LIBPATH 	= -L/home/kddevin/code/ParMETIS3_1/
PARMETIS_INCPATH 	= -I/home/kddevin/code/ParMETIS3_1/
#
#JOSTLE_LIBPATH 		= -L/Net/local/proj/zoltan/arch/solaris/lib
#JOSTLE_INCPATH 		= -I/Net/local/proj/zoltan/arch/all/src
#JOSTLE_SRCPATH 		= /Net/local/proj/zoltan/arch/all/src
#
#PATOH_LIBPATH		= -L/Net/local/proj/zoltan/arch/solaris/lib
#PATOH_INCPATH		= -I/Net/local/proj/zoltan/arch/all/src
#
#NEMESIS_LIBPATH		= -L/Net/local/proj/zoltan/arch/solaris/lib
#NEMESIS_INCPATH		= -I/Net/local/proj/zoltan/arch/solaris/include

