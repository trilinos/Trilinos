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

noinst_PROGRAMS = simpleBLOCK

simpleBLOCK_SOURCES      = simpleBLOCK.cpp
simpleBLOCK_DEPENDENCIES = $(ZOLTAN_DEPS)
simpleBLOCK_CXXFLAGS     = $(ZOLTAN_INCLUDES)
simpleBLOCK_LDADD        = $(ZOLTAN_LIBS)

EXTRA_DIST = classicMakefile objects.txt
