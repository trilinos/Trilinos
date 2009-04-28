#@HEADER
# ************************************************************************
# 
#          Tpetra: Templated Linear Algebra Services Package
#                 Copyright (2008) Sandia Corporation
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

include $(top_builddir)/Makefile.export.tpetra
include $(top_builddir)/../epetra/Makefile.export.epetra
include $(top_builddir)/../triutils/Makefile.export.triutils

EXEEXT = .exe

noinst_PROGRAMS = TpetraOPs

TpetraOPs_SOURCES = $(srcdir)/cxx_main.cpp

TpetraOPs_DEPENDENCIES = $(top_builddir)/src/libtpetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../kokkos/src/libkokkos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../triutils/src/libtriutils.a

_EXPORT_INCLUDES = $(TPETRA_INCLUDES)  $(EPETRA_INCLUDES) $(TRIUTILS_INCLUDES)
_EXPORT_LIBS =  $(TPETRA_LIBS) $(EPETRA_LIBS) $(TRIUTILS_LIBS)

@USING_GNUMAKE_TRUE@EXPORT_INCLUDES  = $(shell @PERL_EXE@ @abs_top_srcdir@/config/strip_dup_incl_paths.pl $(_EXPORT_INCLUDES))
@USING_GNUMAKE_TRUE@EXPORT_LIBS      = $(shell @PERL_EXE@ @abs_top_srcdir@/config/strip_dup_libs.pl $(_EXPORT_LIBS))
@USING_GNUMAKE_FALSE@EXPORT_INCLUDES = $(_EXPORT_INCLUDES)
@USING_GNUMAKE_FALSE@EXPORT_LIBS     = $(_EXPORT_LIBS)

TpetraOPs_CXXFLAGS = $(EXPORT_INCLUDES)

TpetraOPs_LDADD = $(EXPORT_LIBS)
