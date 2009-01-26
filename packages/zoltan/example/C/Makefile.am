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

ZOLTAN_EX_LIBS = -L$(top_builddir)/example/lib -lexzoltan $(ZOLTAN_LIBS)

ZOLTAN_EX_DEPS = $(top_builddir)/src/libzoltan.a $(top_builddir)/example/lib/libexzoltan.a

ZOLTAN_EX_INCLUDES = -I$(top_srcdir)/example/lib $(ZOLTAN_INCLUDES)

noinst_PROGRAMS = simpleBLOCK simpleRCB simpleGRAPH simpleHG zoltanSimple zoltanExample1 phgExample

simpleBLOCK_SOURCES      = simpleBLOCK.c simpleGraph.h simpleQueries.h
simpleBLOCK_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
simpleBLOCK_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
simpleBLOCK_LDADD        = $(ZOLTAN_LIBS)

simpleRCB_SOURCES      = simpleRCB.c simpleGraph.h simpleQueries.h
simpleRCB_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
simpleRCB_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
simpleRCB_LDADD        = $(ZOLTAN_EX_LIBS)

simpleGRAPH_SOURCES      = simpleGRAPH.c simpleGraph.h simpleQueries.h
simpleGRAPH_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
simpleGRAPH_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
simpleGRAPH_LDADD        = $(ZOLTAN_EX_LIBS)

simpleHG_SOURCES      = simpleHG.c simpleGraph.h simpleQueries.h
simpleHG_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
simpleHG_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
simpleHG_LDADD        = $(ZOLTAN_EX_LIBS)

zoltanSimple_SOURCES      = zoltanSimple.c
zoltanSimple_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
zoltanSimple_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
zoltanSimple_LDADD        = $(ZOLTAN_EX_LIBS)


zoltanExample1_SOURCES      = zoltanExample1.c
zoltanExample1_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
zoltanExample1_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
zoltanExample1_LDADD        = $(ZOLTAN_EX_LIBS)

phgExample_SOURCES      = phgExample.c
phgExample_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
phgExample_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
phgExample_LDADD        = $(ZOLTAN_EX_LIBS)

SUBDIRS = sparse_matrix

EXTRA_DIST = makefile
