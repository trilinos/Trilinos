#@HEADER
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
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

EXEEXT = .exe
EXAMPLE = example_AmesosFactory example_AmesosFactory_HB example_AmesosFactory_Tridiag compare_solvers

noinst_PROGRAMS = $(EXAMPLE)

EXTRA_DIST = 662_bus.rsa

#copy these files over to the build area
GetNecessaryFiles: $(EXTRA_DIST)
	        cp -f $(srcdir)/662_bus.rsa .

if HAVE_AMESOS_TRIUTILS
TRIUTILS_INCDIR = -I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src
TRIUTILS_LIBS = $(top_builddir)/../triutils/src/libtriutils.a
else
TRIUTILS_INCDIR = 
TRIUTILS_LIBS = 
endif

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
	      -I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src \
	      ${TRIUTILS_INCDIR}

example_AmesosFactory_SOURCES = $(srcdir)/example_AmesosFactory.cpp

example_AmesosFactory_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	${TRIUTILS_LIBS} \
	GetNecessaryFiles

example_AmesosFactory_CXXFLAGS = -I$(top_srcdir)/src \
				 ${TRIUTILS_INCDIR}

example_AmesosFactory_LDADD = $(top_builddir)/src/libamesos.a \
	$(AMESOS_LIBS) \
	$(TRIUTILS_LIBS) \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)  $(FLIBS)

example_AmesosFactory_Tridiag_SOURCES = $(srcdir)/example_AmesosFactory_Tridiag.cpp

example_AmesosFactory_Tridiag_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	${TRIUTILS_LIBS} \
	GetNecessaryFiles

example_AmesosFactory_Tridiag_CXXFLAGS = -I$(top_srcdir)/src \
				 ${TRIUTILS_INCDIR}

example_AmesosFactory_Tridiag_LDADD = $(top_builddir)/src/libamesos.a \
	$(AMESOS_LIBS) \
	$(TRIUTILS_LIBS) \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)  $(FLIBS)

example_AmesosFactory_HB_SOURCES = $(srcdir)/example_AmesosFactory_HB.cpp

example_AmesosFactory_HB_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	${TRIUTILS_LIBS}

example_AmesosFactory_HB_CXXFLAGS = -I$(top_srcdir)/src \
				 ${TRIUTILS_INCDIR}

example_AmesosFactory_HB_LDADD = $(top_builddir)/src/libamesos.a \
	$(AMESOS_LIBS) \
	$(TRIUTILS_LIBS) \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)  $(FLIBS)

compare_solvers_SOURCES = $(srcdir)/compare_solvers.cpp

compare_solvers_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	${TRIUTILS_LIBS}

compare_solvers_CXXFLAGS = -I$(top_srcdir)/src \
				 ${TRIUTILS_INCDIR}

compare_solvers_LDADD = $(top_builddir)/src/libamesos.a \
	$(AMESOS_LIBS) \
	$(TRIUTILS_LIBS) \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)  $(FLIBS)
