#@HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
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

if BUILD_TSF_EXAMPLES

include $(top_builddir)/../TSFExtended/Makefile.export



noinst_PROGRAMS = kepler ODE

kepler_SOURCES = $(srcdir)/mykeplertest.C 


kepler_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-tsf/libnoxtsf.a 

kepler_CXXFLAGS = $(EX_CXXFLAGS) \
	-I$(top_srcdir)/src -I$(top_srcdir)/src-tsf

kepler_LDADD = \
	$(top_builddir)/src-tsf/libnoxtsf.a \
  $(top_builddir)/src/libnox.a \
	$(TSF_EXPORT_LIBS) \
  $(MPI_LIBS) $(LIBS) 


ODE_SOURCES = $(srcdir)/ODETest.C 


ODE_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-tsf/libnoxtsf.a 

ODE_CXXFLAGS = $(EX_CXXFLAGS) \
	-I$(top_srcdir)/src -I$(top_srcdir)/src-tsf

ODE_LDADD = \
	$(top_builddir)/src-tsf/libnoxtsf.a \
  $(top_builddir)/src/libnox.a \
	$(TSF_EXPORT_LIBS) \
  $(MPI_LIBS) $(LIBS) 



EX_CXXFLAGS = \
	-I$(top_srcdir)/src-loca/src \
	-I$(top_srcdir)/src \
	$(TSFEXTENDED_INCLUDES) \
	$(TSFCORE_INCLUDES)\
	$(AZTECOO_INCLUDES)\
	$(EPETRA_INCLUDES)\
	$(TEUCHOS_INCLUDES)



TSFEXTENDED_INCLUDES =  \
-I$(top_builddir)/../TSFExtended/src \
-I$(top_srcdir)/../TSFExtended/src/Utils \
-I$(top_srcdir)/../TSFExtended/src/Vectors \
-I$(top_srcdir)/../TSFExtended/src/Concrete \
-I$(top_srcdir)/../TSFExtended/src/Operators \
-I$(top_srcdir)/../TSFExtended/src/LinearSolvers \
-I$(top_srcdir)/../TSFExtended/src/Preconditioners \
-I$(top_builddir)/../TSFExtended/src/Utils \
-I$(top_builddir)/../TSFExtended/src/Vectors \
-I$(top_builddir)/../TSFExtended/src/Concrete \
-I$(top_builddir)/../TSFExtended/src/Operators \
-I$(top_builddir)/../TSFExtended/src/LinearSolvers \
-I$(top_builddir)/../TSFExtended/src/Preconditioners 

TSFCORE_INCLUDES =  \
-I$(top_builddir)/../TSFCore/src \
-I$(top_builddir)/../TSFCore/src/interfaces/Core \
-I$(top_builddir)/../TSFCore/src/utilities/Core \
-I$(top_builddir)/../TSFCore/src/adapters/aztecoo \
-I$(top_builddir)/../TSFCore/src/adapters/epetra \
-I$(top_builddir)/../TSFCore/src/adapters/mpi-base \
-I$(top_srcdir)/../TSFCore/src \
-I$(top_srcdir)/../TSFCore/src/interfaces/Core \
-I$(top_srcdir)/../TSFCore/src/utilities/Core \
-I$(top_srcdir)/../TSFCore/src/adapters/aztecoo \
-I$(top_srcdir)/../TSFCore/src/adapters/epetra \
-I$(top_srcdir)/../TSFCore/src/adapters/mpi-base \
-I$(top_srcdir)/../TSFCoreUtils/src \
-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
-I$(top_builddir)/../TSFCoreUtils/src


ML_INCLUDES = \
-I$(top_builddir)/../ml/src/Include\
-I$(top_builddir)/../ml/src\
-I$(top_builddir)/../ml/src/Coarsen\
-I$(top_builddir)/../ml/src/Comm\
-I$(top_builddir)/../ml/src/FEGrid\
-I$(top_builddir)/../ml/src/Krylov\
-I$(top_builddir)/../ml/src/Main\
-I$(top_builddir)/../ml/src/Old\
-I$(top_builddir)/../ml/src/Operator\
-I$(top_builddir)/../ml/src/Runtime\
-I$(top_builddir)/../ml/src/Smoother\
-I$(top_builddir)/../ml/src/Utils\
-I$(top_srcdir)/../ml/src/Include\
-I$(top_srcdir)/../ml/src\
-I$(top_srcdir)/../ml/src/Coarsen\
-I$(top_srcdir)/../ml/src/Comm\
-I$(top_srcdir)/../ml/src/FEGrid\
-I$(top_srcdir)/../ml/src/Krylov\
-I$(top_srcdir)/../ml/src/Main\
-I$(top_srcdir)/../ml/src/Old\
-I$(top_srcdir)/../ml/src/Operator\
-I$(top_srcdir)/../ml/src/Runtime\
-I$(top_srcdir)/../ml/src/Smoother\
-I$(top_srcdir)/../ml/src/Utils

EPETRA_INCLUDES = \
-I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

AZTECOO_INCLUDES = \
-I$(top_srcdir)/../aztecoo/src -I$(top_builddir)/../aztecoo/src

TEUCHOS_INCLUDES = \
-I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src


endif