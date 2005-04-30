#@HEADER
# ************************************************************************
# 
#               Epetra: Linear Algebra Services Package 
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

include $(top_builddir)/Makefile.export.epetra.thyra

EXEEXT = .exe

#
# Actual programs
#

TEST_EPETRA_ADAPTERS = test_epetra_adapters

noinst_PROGRAMS = $(TEST_EPETRA_ADAPTERS)

test_epetra_adapters_SOURCES = $(srcdir)/test_epetra_adapters.cpp $(common_sources)
test_epetra_adapters_DEPENDENCIES = $(common_dependencies)
test_epetra_adapters_LDADD = $(common_ldadd)

#
# Common stuff
#

common_sources =

common_dependencies = $(top_builddir)/thyra/src/libEpetraThyra.a $(top_builddir)/src/libepetra.a

AM_CPPFLAGS = $(shell perl $(top_srcdir)/config/strip_dup_incl_paths.pl \
	-I$(top_srcdir)/../thyra/example/ana/operator_vector \
	$(EPETRA_THYRA_INCLUDES) \
	)

common_ldadd = $(shell perl $(top_srcdir)/config/strip_dup_libs.pl $(EPETRA_THYRA_LIBS))

# This is already added part of EPETRA_LIBS and therefore automake does not need to add this again!
LIBS =
