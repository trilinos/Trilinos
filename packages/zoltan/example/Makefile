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
# Makefile for simple examples which use the Zoltan library
#
# To compile the examples, type
#    gmake ZOLTAN_ARCH=<targetmachine> 
# where
#    targetmachine = the suffix of one of the configuration files Config.*
#                    in Zoltan/Utilities/Config.
##############################################################################

#
# relative to this directory
#
ZOLTAN_TOP = ..
ZOLTAN_UTIL = ../Utilities

#
# relative to the subdirectories
#
ZTOP = ../$(ZOLTAN_TOP)
ZUTIL = ../$(ZOLTAN_UTIL)

##############################################################################
# Process architecture flag
##############################################################################

ifndef ZOLTAN_ARCH
ZOLTAN_ARCH	= generic
endif

ifndef ZOLTAN_CONFIG
ZOLTAN_CONFIG = $(ZOLTAN_TOP)/Utilities/Config/Config.$(ZOLTAN_ARCH)
endif

include $(ZOLTAN_CONFIG)

ifndef ZOLTAN_OBJ_DIR
ZOLTAN_OBJ_DIR = $(ZTOP)/Obj_$(ZOLTAN_ARCH)
endif

ifndef ZOLTAN_LIB_DIR
ZOLTAN_LIB_DIR = $(ZOLTAN_OBJ_DIR)
endif

############################################################################
#                 SUB DIRECTORIES
############################################################################
#

SUBDIRS = C CPP lib

############################################################################
#                 LIBRARIES
############################################################################
#
ifdef PARMETIS_LIBPATH
  PARMETIS_LIBS = -lparmetis -lmetis
endif

LIBDIR = -L../lib -L$(ZOLTAN_LIB_DIR) $(PARMETIS_LIBPATH) $(MPI_LIBPATH) 

LIBS = -lexzoltan -lzoltan $(PARMETIS_LIBS) $(MPI_LIB) -lm

ZOLTAN_INCLUDE= -I$(ZUTIL)/Communication \
                -I$(ZUTIL)/DDirectory \
                -I$(ZUTIL)/Memory \
                -I$(ZUTIL)/Timer \
                -I$(ZUTIL)/shared \
                -I$(ZTOP)/zz \
                -I$(ZTOP)/par \
                -I$(ZTOP)/lb \
                -I$(ZTOP)/include

INCLUDE = -I../lib $(ZOLTAN_INCLUDE)

##############################################################################
# 
##############################################################################

ifdef CPPC
TARGETS= exampleLibrary C_Examples CPP_Examples
else
TARGETS= exampleLibrary C_Examples
endif

all: $(TARGETS)

exampleLibrary:
	cd lib; $(MAKE) all "CC=$(CC)" "INCLUDE=$(INCLUDE)" 

C_Examples:
	cd C; $(MAKE) all "CC=$(CC)" "INCLUDE=$(INCLUDE)" \
                          "LIBDIR=$(LIBDIR)" "LIBS=$(LIBS)"

CPP_Examples:
	cd CPP; $(MAKE) all "CC=$(CPPC)" "INCLUDE=$(INCLUDE)" \
                          "LIBDIR=$(LIBDIR)" "LIBS=$(LIBS)"


clean:
	cd lib; $(MAKE) clean
	cd C; $(MAKE) clean
	cd CPP; $(MAKE) clean

