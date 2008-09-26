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


noinst_PROGRAMS = ex1 

ex1_SOURCES      = ex1.c queries1.c matrix.h 
ex1_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
ex1_CFLAGS     = $(ZOLTAN_EX_INCLUDES)
ex1_LDADD        = $(ZOLTAN_EX_LIBS)

EXTRA_DIST = makefile
