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
noinst_PROGRAMS = rpmtest

#np# One SOURCES primary is needed for EACH exe listed above.  (Notice for all 
#np# of the  following primaries that the pattern is 
#np# "(Exe_name)_(PRIMARY_NAME) = ").
#np# List all source files needed to compile the exe.
rpmtest_SOURCES = $(srcdir)/rpmtest_main.cpp $(srcdir)/rpmtest_integrator.cpp $(srcdir)/rpmtest_interface.cpp

#np# One DEPENDENCIES primary is needed for each exe listed above.  List any
#np# libraries that the exe depends on.  (This primary is not absolutely
#np# essential.)
rpmtest_DEPENDENCIES = $(top_builddir)/src/libcapo.a \
	$(top_builddir)/../epetra/src/libepetra.a

#np# One CXXFLAGS primary can be defined for each exe listed above.  Include 
#np# any -I's necessary for building the exe.
rpmtest_CXXFLAGS = -I$(top_srcdir)/src -I$(srcdir)/.. \
	-I$(top_srcdir)/../epetra/src \
	-I$(top_builddir)/../epetra/src \
	-I$(top_srcdir)/../epetra/thyra/src \
	-I$(top_builddir)/../epetra/thyra/src \
	-I$(top_builddir)/../thyra/src \
	-I$(top_srcdir)/../thyra/src/ana/operator_vector/interfaces \
	-I$(top_srcdir)/../thyra/src/ana/operator_vector/extended_interfaces \
	-I$(top_srcdir)/../thyra/src/ana/operator_vector/ana_support \
	-I$(top_srcdir)/../thyra/src/ana/operator_vector/adapter_support \
	-I$(top_srcdir)/../thyra/src/ana/operator_vector/adapter_support/mpi \
	-I$(top_srcdir)/../thyra/src/ana/operator_vector/adapter_support/serial \
	-I$(top_srcdir)/../thyra/src/ana/rtop \
	-I$(top_srcdir)/../thyra/src/ana/rtop/ops_lib \
	-I$(top_srcdir)/../thyra/src/ana/rtop/base \
	-I$(top_srcdir)/../thyra/src/utils \
	-I$(top_srcdir)/../teuchos/src \
	-I$(top_builddir)/../thyra/src/ana/operator_vector/interfaces \
	-I$(top_builddir)/../thyra/src/ana/operator_vector/extended_interfaces \
	-I$(top_builddir)/../thyra/src/ana/operator_vector/ana_support \
	-I$(top_builddir)/../thyra/src/ana/operator_vector/adapter_support \
	-I$(top_builddir)/../thyra/src/ana/operator_vector/adapter_support/mpi \
	-I$(top_builddir)/../thyra/src/ana/operator_vector/adapter_support/serial \
	-I$(top_builddir)/../thyra/src/ana/rtop \
	-I$(top_builddir)/../thyra/src/ana/rtop/ops_lib \
	-I$(top_builddir)/../thyra/src/ana/rtop/base \
	-I$(top_builddir)/../thyra/src/utils \
	-I$(top_builddir)/../teuchos/src \
	-I$(top_builddir)/../epetra/thyra/src \
	-I$(top_builddir)/../epetra/src \
	-I$(top_srcdir)/../thyra/src			

#np# One LDADD primary can be defined for each exe listed above.  Include any
#np# LDFLAGS  necessary for building the exe.
rpmtest_LDADD = $(top_builddir)/src/libcapo.a \
	$(top_builddir)/../epetra/thyra/src/libepetrathyra.a \
	$(top_builddir)/../thyra/src/libthyra.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	-llapack -lblas -lg2c
