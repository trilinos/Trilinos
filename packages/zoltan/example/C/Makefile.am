# @HEADER
##############################################################################
# Zoltan Library for Parallel Applications                                   #
# Copyright (c) 2000,2001,2002, Sandia National Laboratories.                #
# This software is distributed under the GNU Lesser General Public License.  #
# For more info, see the README file in the top-level Zoltan directory.      #
##############################################################################
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:

include $(top_builddir)/Makefile.export.zoltan

EXEEXT = .exe

ZOLTAN_EX_DEPS = $(top_builddir)/src/libzoltan.a 

ZOLTAN_EX_INCLUDES = $(ZOLTAN_INCLUDES)

noinst_PROGRAMS = simpleBLOCK simpleRCB simpleGRAPH simplePHG 

simpleBLOCK_SOURCES      = simpleBLOCK.c
simpleBLOCK_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
simpleBLOCK_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
simpleBLOCK_LDADD        = $(ZOLTAN_LIBS)

simpleRCB_SOURCES      = simpleRCB.c
simpleRCB_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
simpleRCB_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
simpleRCB_LDADD        = $(ZOLTAN_EX_LIBS)

simpleGRAPH_SOURCES      = simpleGRAPH.c
simpleGRAPH_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
simpleGRAPH_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
simpleGRAPH_LDADD        = $(ZOLTAN_EX_LIBS)

simplePHG_SOURCES      = simplePHG.c
simplePHG_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
simplePHG_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
simplePHG_LDADD        = $(ZOLTAN_EX_LIBS)

EXTRA_DIST = classicMakefile objects.txt hypergraph.txt graph.txt mesh.txt
