## $Id$ 
## $Source$ 
#@HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
#                 Copyright (2005) Sandia Corporation
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
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
#@HEADER

# ------------------------------------------------------------------------
# FEApp Sources
# ------------------------------------------------------------------------
include $(top_srcdir)/../sacado/example/FEApp/MakefileFEAppStub.am

# Include export Makefiles for correct paths
include $(top_builddir)/Makefile.export.nox

EXEEXT = .exe

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(LOCA_EPETRA_TEST_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_EPETRA_TEST_INCLUDES))
else
TEST_LIBS     = $(LOCA_EPETRA_TEST_LIBS)
TEST_INCLUDES = $(LOCA_EPETRA_TEST_INCLUDES)
endif

if BUILD_THYRA_TESTS
if BUILD_EPETRA_TESTS
if ENABLE_LOCA_SACADO
if BUILD_NOX_EPETRAEXT_SUPPORT
include $(top_builddir)/../sacado/Makefile.export.sacado
noinst_PROGRAMS = tcubed tcubedTP
endif
endif
endif
endif

tcubed_SOURCES = \
	$(srcdir)/tcubed.C 

tcubed_CXXFLAGS = \
	-I$(FEAPP_PATH) \
	$(TEST_INCLUDES) \
	$(SACADO_INCLUDES)

tcubed_LDADD = \
	$(TEST_LIBS) \
	$(SACADO_LIBS)

tcubed_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-thyra/libnoxthyra.a \
	$(top_builddir)/test/utils/libnoxtestutils.a \
	$(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-thyra/liblocathyra.a \
	$(top_builddir)/test/epetra/LOCA_TestProblems/liblocatestproblems.a

tcubedTP_SOURCES = $(srcdir)/tcubedTP.C 

tcubedTP_CXXFLAGS = $(tcubed_CXXFLAGS)

tcubedTP_LDADD = $(tcubed_LDADD)

tcubedTP_DEPENDENCIES = $(tcubed_DEPENDENCIES)