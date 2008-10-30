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

#np# Specify locations of include files within Trilinos, but outside np
#np# Add packages that np is dependent on; remove epetra if np does not need it

include $(top_builddir)/Makefile.export.anasazi

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ANASAZI_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ANASAZI_TEST_INCLUDES))
else
EXPORT_LIBS = $(ANASAZI_TEST_LIBS)
EXPORT_INC_PATH = $(ANASAZI_TEST_INCLUDES)
endif

AM_CPPFLAGS = \
	$(EXPORT_INC_PATH)

EXEEXT = .exe

if BUILD_TEUCHOS_COMPLEX
BKS_COMPLEX_LAP_TEST = Tpetra_BlockKrylovSchur_complex_lap_test
if BUILD_TRIUTILS
BKS_COMPLEX_TEST = Tpetra_BlockKrylovSchur_complex_test Dummy
endif
endif

#np# List all executable names here
noinst_PROGRAMS = $(BKS_COMPLEX_TEST) $(BKS_COMPLEX_LAP_TEST)

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/../../../test/MVOPTester/mhd1280a.cua mhd1280a.cua
	cp  $(srcdir)/../../../test/MVOPTester/mhd1280b.cua mhd1280b.cua

#np# One SOURCES primary is needed for EACH exe.  (Notice for all of the following 
#np# primaries that the pattern is "(Exe_name)_(PRIMARY_NAME) = ").  List all source 
#np# files needed to compile the exe.
Tpetra_BlockKrylovSchur_complex_test_SOURCES = $(srcdir)/cxx_main_complex.cpp
Tpetra_BlockKrylovSchur_complex_lap_test_SOURCES = $(srcdir)/cxx_main_complex_lap.cpp

Dummy_SOURCES = Dummy

#np# One DEPENDENCIES primary is needed for each exe.  List any libraries that the 
#np# exe depends on.  (This primary is not absolutely essential.)
Tpetra_BlockKrylovSchur_complex_test_DEPENDENCIES = \
	$(top_builddir)/src/libanasazi.a \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

Tpetra_BlockKrylovSchur_complex_lap_test_DEPENDENCIES = \
	$(top_builddir)/src/libanasazi.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

#np# One CXXFLAGS primary can be defined for each exe.  Include any -I's necessary for 
#np# building the exe.
Tpetra_BlockKrylovSchur_complex_test_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/test/MVOPTester

Tpetra_BlockKrylovSchur_complex_lap_test_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/test/MVOPTester

#np# One LDADD primary can be defined for each exe.  Include any LDFLAGS necessary for 
#np# building the exe.
Tpetra_BlockKrylovSchur_complex_test_LDADD = \
	$(EXPORT_LIBS)

Tpetra_BlockKrylovSchur_complex_lap_test_LDADD = \
	$(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 
