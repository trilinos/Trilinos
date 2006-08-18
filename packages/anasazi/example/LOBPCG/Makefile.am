#@HEADER
# ************************************************************************
#
#                 Anasazi: Block Eigensolvers Package
#                 Copyright (2004) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.anasazi

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ANASAZI_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ANASAZI_TEST_INCLUDES))
else
EXPORT_LIBS = $(ANASAZI_TEST_LIBS)
EXPORT_INC_PATH = $(ANASAZI_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)  

EXEEXT = .exe

noinst_PROGRAMS = LOBPCGEpetraEx LOBPCGEpetraExSimple LOBPCGEpetraExGen

LOBPCGEpetraEx_SOURCES = $(srcdir)/LOBPCGEpetraEx.cpp

LOBPCGEpetraExSimple_SOURCES = $(srcdir)/LOBPCGEpetraExSimple.cpp

LOBPCGEpetraExGen_SOURCES = $(srcdir)/LOBPCGEpetraExGen.cpp

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
LOBPCGEpetraEx_DEPENDENCIES = \
        $(top_builddir)/src/libanasazi.a \
        $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a 

LOBPCGEpetraEx_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/util/ModeLaplace

#np# List libraries that np is dependent on.
LOBPCGEpetraEx_LDADD = $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
                       $(EXPORT_LIBS) 

LOBPCGEpetraExSimple_DEPENDENCIES = \
        $(top_builddir)/src/libanasazi.a \
        $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a 

LOBPCGEpetraExSimple_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/util/ModeLaplace

#np# List libraries that np is dependent on.
LOBPCGEpetraExSimple_LDADD = $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
                       $(EXPORT_LIBS) 

LOBPCGEpetraExGen_DEPENDENCIES = \
        $(top_builddir)/src/libanasazi.a \
        $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

LOBPCGEpetraExGen_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/util/ModeLaplace

#np# List libraries that np is dependent on.
LOBPCGEpetraExGen_LDADD = $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
                          $(EXPORT_LIBS)

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
SUBDIRS = 

# Dummy for new export makefiles
LIBS = 
