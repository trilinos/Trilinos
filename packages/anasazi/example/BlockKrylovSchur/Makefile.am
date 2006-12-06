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

if BUILD_BELOS
if BUILD_IFPACK
ANASAZI_BIRAM_EX_SI_BELOS = BlockKrylovSchurEpetraExGenBelos
endif
endif

if BUILD_AZTECOO
if BUILD_IFPACK
ANASAZI_BIRAM_EX_SI_AZTECOO = BlockKrylovSchurEpetraExGenAztecOO
endif
endif

if BUILD_AMESOS
ANASAZI_BIRAM_EX_SI_AMESOS = BlockKrylovSchurEpetraExGenAmesos
endif

if BUILD_EPETRAEXT
ANASAZI_BIRAM_EX_FILE = BlockKrylovSchurEpetraExFile
endif 

noinst_PROGRAMS = BlockKrylovSchurEpetraEx BlockKrylovSchurEpetraExSVD $(ANASAZI_BIRAM_EX_FILE) $(ANASAZI_BIRAM_EX_SI_AZTECOO) $(ANASAZI_BIRAM_EX_SI_BELOS) $(ANASAZI_BIRAM_EX_SI_AMESOS)

BlockKrylovSchurEpetraEx_SOURCES = $(srcdir)/BlockKrylovSchurEpetraEx.cpp

BlockKrylovSchurEpetraExSVD_SOURCES = $(srcdir)/BlockKrylovSchurEpetraExSVD.cpp

BlockKrylovSchurEpetraExFile_SOURCES = $(srcdir)/BlockKrylovSchurEpetraExFile.cpp

BlockKrylovSchurEpetraExGenAztecOO_SOURCES = $(srcdir)/BlockKrylovSchurEpetraExGenAztecOO.cpp

BlockKrylovSchurEpetraExGenBelos_SOURCES = $(srcdir)/BlockKrylovSchurEpetraExGenBelos.cpp

BlockKrylovSchurEpetraExGenAmesos_SOURCES = $(srcdir)/BlockKrylovSchurEpetraExGenAmesos.cpp

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
BlockKrylovSchurEpetraEx_DEPENDENCIES = \
	$(top_builddir)/src/libanasazi.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a 

BlockKrylovSchurEpetraEx_CXXFLAGS = -I$(top_srcdir)/src

#np# List libraries that np is dependent on.
BlockKrylovSchurEpetraEx_LDADD = \
	$(top_builddir)/src/libanasazi.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

BlockKrylovSchurEpetraExSVD_DEPENDENCIES = \
	$(top_builddir)/src/libanasazi.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a 

BlockKrylovSchurEpetraExSVD_CXXFLAGS = -I$(top_srcdir)/src

#np# List libraries that np is dependent on.
BlockKrylovSchurEpetraExSVD_LDADD = $(EXPORT_LIBS)

BlockKrylovSchurEpetraExFile_DEPENDENCIES = \
        $(EXPORT_LIBS)

BlockKrylovSchurEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src

#np# List libraries that np is dependent on.
BlockKrylovSchurEpetraExFile_LDADD = $(EXPORT_LIBS)

BlockKrylovSchurEpetraExGenAztecOO_DEPENDENCIES = \
	$(top_builddir)/src/libanasazi.a \
	$(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../aztecoo/src/libaztecoo.a \
	$(top_builddir)/../ifpack/src/libifpack.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

BlockKrylovSchurEpetraExGenAztecOO_CXXFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/util/ModeLaplace

#np# List libraries that np is dependent on.
BlockKrylovSchurEpetraExGenAztecOO_LDADD = \
	$(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(EXPORT_LIBS)

BlockKrylovSchurEpetraExGenBelos_DEPENDENCIES = \
	$(top_builddir)/src/libanasazi.a \
	$(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../belos/src/libbelos.a \
	$(top_builddir)/../ifpack/src/libifpack.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

BlockKrylovSchurEpetraExGenBelos_CXXFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/util/ModeLaplace

#np# List libraries that np is dependent on.
BlockKrylovSchurEpetraExGenBelos_LDADD = \
	$(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(EXPORT_LIBS)

BlockKrylovSchurEpetraExGenAmesos_DEPENDENCIES = \
	$(top_builddir)/src/libanasazi.a \
	$(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../amesos/src/libamesos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

BlockKrylovSchurEpetraExGenAmesos_CXXFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/util/ModeLaplace

#np# List libraries that np is dependent on.
BlockKrylovSchurEpetraExGenAmesos_LDADD = \
	$(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(EXPORT_LIBS) $(TEUCHOS_MPI_LIBS) $(TEUCHOS_BLAS_LIBS)

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
