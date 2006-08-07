#@HEADER
# ************************************************************************
# 
#                  LOCA Continuation Algorithm Package
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
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov), Sandia National Laboratories.
# 
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(LOCA_EPETRA_TEST_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_EPETRA_TEST_INCLUDES))
else
TEST_LIBS     = $(LOCA_EPETRA_TEST_LIBS)
TEST_INCLUDES = $(LOCA_EPETRA_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(TEST_INCLUDES)

EXEEXT = .exe

# The Hopf test requries epetraext
if BUILD_NOX_EPETRAEXT_SUPPORT
BRUSSHOPFCONT = brussHopfCont
else
BRUSSHOPFCONT = 
endif

if BUILD_LOCA_EPETRA_TESTS
noinst_PROGRAMS = brussCont $(BRUSSHOPFCONT)
endif

brussCont_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussCont.C

brussCont_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
        $(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

brussCont_LDADD = $(TEST_LIBS)

brussHopfCont_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussHopfCont.C

brussHopfCont_DEPENDENCIES = $(brussCont_DEPENDENCIES)

brussHopfCont_LDADD = $(TEST_LIBS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

