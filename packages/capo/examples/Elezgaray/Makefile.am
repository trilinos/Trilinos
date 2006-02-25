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

#np# List the names of all exe's in this directory here
noinst_PROGRAMS = elezgaray

#np# One SOURCES primary is needed for EACH exe listed above.  (Notice for all 
#np# of the  following primaries that the pattern is 
#np# "(Exe_name)_(PRIMARY_NAME) = ").
#np# List all source files needed to compile the exe.
elezgaray_SOURCES = $(srcdir)/EA_main.cpp $(srcdir)/EA_integrator.cpp $(srcdir)/EA_interface.cpp

#np# One DEPENDENCIES primary is needed for each exe listed above.  List any
#np# libraries that the exe depends on.  (This primary is not absolutely
#np# essential.)
elezgaray_DEPENDENCIES = $(top_builddir)/src/libcapo.a \
	$(top_builddir)/../epetra/src/libepetra.a

#np# One CXXFLAGS primary can be defined for each exe listed above.  Include 
#np# any -I's necessary for building the exe.
elezgaray_CXXFLAGS = \
        -I$(srcdir)/.. \
        -I$(top_srcdir)/src \
	-I$(top_srcdir)/../thyra/src \
	-I$(top_srcdir)/../thyra/src/utils \
	-I$(top_srcdir)/../thyra/src/interfaces/operator_vector/ana/fundamental \
	-I$(top_srcdir)/../thyra/src/support/operator_vector/client_support \
	-I$(top_srcdir)/../thyra/src/support/operator_vector/adapter_support \
	-I$(top_srcdir)/../thyra/src/support/operator_vector/adapter_support/mpi \
	-I$(top_srcdir)/../rtop/src \
	-I$(top_srcdir)/../rtop/src/ops_lib \
	-I$(top_srcdir)/../rtop/src/support \
	-I$(top_srcdir)/../rtop/src/interfaces \
	-I$(top_srcdir)/../teuchos/src \
	-I$(top_srcdir)/../epetra/src \
	-I$(top_srcdir)/../epetra/thyra/src \
	-I$(top_builddir)/../thyra/src \
	-I$(top_builddir)/../thyra/src/utils \
	-I$(top_builddir)/../thyra/src/interfaces/operator_vector/ana/fundamental \
	-I$(top_builddir)/../rtop/src \
	-I$(top_builddir)/../teuchos/src \
	-I$(top_builddir)/../epetra/thyra/src \
	-I$(top_builddir)/../epetra/src \
	-I$(top_builddir)/../epetra/thyra/src

#np# One LDADD primary can be defined for each exe listed above.  Include any
#np# LDFLAGS  necessary for building the exe.
elezgaray_LDADD = $(top_builddir)/src/libcapo.a \
	$(top_builddir)/../epetra/thyra/src/libepetrathyra.a \
	$(top_builddir)/../thyra/src/libthyra.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../rtop/src/librtop.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LIBS) $(BLAS_LIBS) $(FLIBS)
