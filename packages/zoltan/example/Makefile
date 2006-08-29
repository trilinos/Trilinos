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
LIBDIR = -L../lib -L$(ZOLTAN_LIB_DIR) $(MPI_LIBPATH)

THIRD_PARTY_LIBS=

ifdef PARMETIS_LIBPATH
  THIRD_PARTY_LIBS += -lparmetis -lmetis
  LIBDIR += $(PARMETIS_LIBPATH)
endif

ifdef JOSTLE_LIBPATH
  THIRD_PARTY_LIBS += -ljostle
  LIBDIR += $(JOSTLE_LIBPATH)
endif

ifdef PATOH_LIBPATH
  THIRD_PARTY_LIBS += -lpatoh
  LIBDIR += $(PATOH_LIBPATH)
endif


LIBS = -lexzoltan -lzoltan $(THIRD_PARTY_LIBS) $(MPI_LIB) -lm

ZOLTAN_INCLUDE= -I$(ZUTIL)/Communication \
                -I$(ZUTIL)/DDirectory \
                -I$(ZUTIL)/Memory \
                -I$(ZUTIL)/Timer \
                -I$(ZUTIL)/shared \
                -I$(ZTOP)/zz \
                -I$(ZTOP)/ha \
                -I$(ZTOP)/par \
                -I$(ZTOP)/lb \
                -I$(ZTOP)/include

INCLUDE = -I../lib $(ZOLTAN_INCLUDE) $(INCLUDE_PATH)

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
	@cd lib; $(MAKE) all "CC=$(CC)" "AR=$(AR)" "CFLAGS=$(CFLAGS)"\
                            "INCLUDE=$(INCLUDE)" 

C_Examples: exampleLibrary
	@cd C; $(MAKE) all "ZOLTAN_LIB_DIR=$(ZOLTAN_LIB_DIR)" \
           "CC=$(CC)" "CFLAGS=$(CFLAGS)" "INCLUDE=$(INCLUDE)" \
                          "LIBDIR=$(LIBDIR)" "LIBS=$(LIBS)"

CPP_Examples: exampleLibrary
	@cd CPP; $(MAKE) all "ZOLTAN_LIB_DIR=$(ZOLTAN_LIB_DIR)" \
           "CC=$(CPPC)" "CFLAGS=$(CFLAGS)" "INCLUDE=$(INCLUDE)" \
                          "LIBDIR=$(LIBDIR)" "LIBS=$(LIBS)"


clean:
	@cd lib; $(MAKE) clean
	@cd C; $(MAKE) clean
	@cd CPP; $(MAKE) clean

