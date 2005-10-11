#@HEADER
# ************************************************************************
# 
#                     New_Package Example Package
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

#np# This file is an example of a typical Makefile.am to build an example exe. 

EXEEXT = .exe

#
# Actual programs
#

noinst_PROGRAMS = vpl
vpl_SOURCES = $(srcdir)/vpl_main.cpp $(srcdir)/vpl_integrator.cpp $(srcdir)/vpl_interface.cpp $(srcdir)/IntegratorT.cpp $(srcdir)/StiffIntegratorT.cpp $(srcdir)/decsol.cpp
vpl_DEPENDENCIES = $(common_dependencies)
vpl_LDADD = $(common_ldadd) \
	$(top_builddir)/src/libcapo.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../epetra/thyra/src/libepetrathyra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../thyra/src/libthyra.a \
	$(LIBS) $(BLAS_LIBS) $(FLIBS)

#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/src/libcapo.a

include $(top_builddir)/Makefile.export.capo

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(CAPO_EXTRA_INCLUDES) $(CAPO_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(CAPO_EXTRA_LIBS) $(CAPO_LIBS))
else
EXPORT_INCLUDES = $(CAPO_EXTRA_INCLUDES) $(CAPO_INCLUDES)
EXPORT_LIBS = $(CAPO_EXTRA_LIBS) $(CAPO_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES) $(EPETRA_INCLUDES) $(CAPO_EXTRA_INCLUDES) $(CAPO_INCLUDES) $(TEUCHOS_INCLUDES) -I$(srcdir)/../../src

common_ldadd = $(EXPORT_LIBS) $(EPETRA_LIBS) $(CAPO_EXTRA_LIBS) $(CAPO_LIBS) $(TEUCHOS_LIBS) $(THYRA_LIBS) $(NOX_LIBS)

# This is already added part of CAPO_LIBS and therefore automake does not need to add this again!
LIBS =
