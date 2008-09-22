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

noinst_PROGRAMS = zCPPExample1 zCPPExample2 simpleRCB 

zCPPExample1_SOURCES      = zCPPExample1.cpp
zCPPExample1_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
zCPPExample1_CXXFLAGS     = $(ZOLTAN_EX_INCLUDES)
zCPPExample1_LDADD        = $(ZOLTAN_EX_LIBS)

zCPPExample2_SOURCES      = zCPPExample2.cpp
zCPPExample2_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
zCPPExample2_CXXFLAGS     = $(ZOLTAN_EX_INCLUDES)
zCPPExample2_LDADD        = $(ZOLTAN_EX_LIBS)

simpleRCB_SOURCES      = simpleRCB.cpp rectangularMesh.h
simpleRCB_DEPENDENCIES = $(ZOLTAN_EX_DEPS)
simpleRCB_CXXFLAGS     = $(ZOLTAN_EX_INCLUDES)
simpleRCB_LDADD        = $(ZOLTAN_EX_LIBS)

EXTRA_DIST = makefile
