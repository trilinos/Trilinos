## $Id$ 
## $Source$ 
# @HEADER
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
# @HEADER

# ------------------------------------------------------------------------
# Swig interface files
# ------------------------------------------------------------------------

NOX_INTERFACES = \
	NOX_Abstract.i \
	NOX_Parameter.i \
	NOX_StatusTest.i \
	NOX_Solver.i 

if BUILD_LIBNOXEPETRA

EPETRA_INTERFACES = \
	NOX_Epetra.i

else

EPETRA_INTERFACES =

endif

if BUILD_LIBNOXLAPACK

LAPACK_INTERFACES = \
	NOX_LAPACK.i

else

LAPACK_INTERFACES =

endif

INTERFACES = $(NOX_INTERFACES) $(EPETRA_INTERFACES) $(LAPACK_INTERFACES)

# ------------------------------------------------------------------------
# Generated swig wrapper files
# ------------------------------------------------------------------------

NOX_WRAPPERS = \
	$(srcdir)/Abstract.py \
	$(srcdir)/Parameter.py \
	$(srcdir)/StatusTest.py \
	$(srcdir)/Solver.py

if BUILD_LIBNOXEPETRA

EPETRA_WRAPPERS = \
	$(srcdir)/Epetra.py 

else

EPETRA_WRAPPERS =

endif

if BUILD_LIBNOXLAPACK

LAPACK_WRAPPERS = \
	$(srcdir)/LAPACK.py 

else

LAPACK_WRAPPERS =

endif

WRAPPERS = $(NOX_WRAPPERS) $(EPETRA_WRAPPERS) $(LAPACK_WRAPPERS)

# ------------------------------------------------------------------------
# Python extension modules as libtool shared libraries
# ------------------------------------------------------------------------

NOX_EXTENSIONS = \
	_Abstract.la \
	_Parameter.la \
	_StatusTest.la \
	_Solver.la

if BUILD_LIBNOXEPETRA

EPETRA_EXTENSIONS = \
	_Epetra.la 

else

EPETRA_EXTENSIONS =

endif

if BUILD_LIBNOXLAPACK

LAPACK_EXTENSIONS = \
	_LAPACK.la 

else

LAPACK_EXTENSIONS =

endif

EXTENSIONS = $(NOX_EXTENSIONS) $(EPETRA_EXTENSIONS) $(LAPACK_EXTENSIONS)

# ------------------------------------------------------------------------
# Where we want the NOX/Python installed.  We use 'exec' in the name so 
# that it is installed by 'install-exec'
# ------------------------------------------------------------------------

noxPythonexecdir = $(prefix)/python/NOX

if BUILD_PYTHON

# ------------------------------------------------------------------------
# What to build
# ------------------------------------------------------------------------

noxPythonexec_LTLIBRARIES = $(EXTENSIONS)

# ------------------------------------------------------------------------
# Python files to be installed by 'make install-exec'
# ------------------------------------------------------------------------

# __init.py__ is not distributed because it is created by configure
nodist_noxPythonexec_PYTHON = __init__.py

dist_noxPythonexec_PYTHON = $(WRAPPERS)

endif

# ------------------------------------------------------------------------
# Extra files to be installed by 'make dist'
# ------------------------------------------------------------------------

EXTRA_DIST = $(INTERFACES) NOX_LAPACK_import.i NOX_LAPACK.swi

# ------------------------------------------------------------------------
# General definitions
# ------------------------------------------------------------------------

CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS)

CXXFLAGS = -Wall -O2

SWIG_PYTHON_LIBS =

SWIG_NOX_CXXFLAGS = $(TRILINOS_CXXFLAGS) \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/python 

SWIG_EPETRA_CXXFLAGS = \
	-I$(top_srcdir)/src-epetra 

SWIG_LAPACK_CXXFLAGS = \
	-I$(top_srcdir)/src-lapack

SWIG_CXXFLAGS = \
	$(SWIG_NOX_CXXFLAGS) \
	$(SWIG_EPETRA_CXXFLAGS) \
	$(SWIG_LAPACK_CXXFLAGS)

SWIG_LDFLAGS = -avoid-version -module

SWIG_NOX_LIBADD = \
	$(top_builddir)/src/libnox.a \
	$(SWIG_PYTHON_LIBS) $(LIBS) $(XTRALDADD)

SWIG_EPETRA_LIBADD = \
	$(top_builddir)/src-epetra/libnoxepetra.a 

SWIG_LAPACK_LIBADD = \
	$(top_builddir)/src-lapack/libnoxlapack.a 

# ------------------------------------------------------------------------
# Definitions for _Abstract.la
# ------------------------------------------------------------------------

_Abstract_la_SOURCES = $(srcdir)/wrap_NOX_Abstract.C

_Abstract_la_CXXFLAGS = $(SWIG_NOX_CXXFLAGS)

_Abstract_la_LDFLAGS = $(SWIG_LDFLAGS)

_Abstract_la_LIBADD = $(SWIG_NOX_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _Parameter.la
# ------------------------------------------------------------------------

_Parameter_la_SOURCES = \
	$(srcdir)/wrap_NOX_Parameter.C \
	$(srcdir)/Utils_enums.H

_Parameter_la_CXXFLAGS = $(SWIG_NOX_CXXFLAGS)

_Parameter_la_LDFLAGS = $(SWIG_LDFLAGS)

_Parameter_la_LIBADD = $(SWIG_NOX_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _StatusTest.la
# ------------------------------------------------------------------------

_StatusTest_la_SOURCES = $(srcdir)/wrap_NOX_StatusTest.C

_StatusTest_la_CXXFLAGS = $(SWIG_NOX_CXXFLAGS)

_StatusTest_la_LDFLAGS = $(SWIG_LDFLAGS)

_StatusTest_la_LIBADD = $(SWIG_NOX_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _Solver.la
# ------------------------------------------------------------------------

_Solver_la_SOURCES = $(srcdir)/wrap_NOX_Solver.C

_Solver_la_CXXFLAGS = $(SWIG_NOX_CXXFLAGS)

_Solver_la_LDFLAGS = $(SWIG_LDFLAGS)

_Solver_la_LIBADD = $(SWIG_NOX_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _Epetra.la
# ------------------------------------------------------------------------

_Epetra_la_SOURCES = \
	$(srcdir)/wrap_NOX_Epetra.C \
	$(srcdir)/Epetra_VectorHelper.H \
	$(srcdir)/Epetra_VectorHelper.C \
	$(srcdir)/NumPyArray.H \
	$(srcdir)/NumPyArray.C \
	$(srcdir)/NumPyWrapper.H \
	$(srcdir)/NumPyWrapper.C \
	$(srcdir)/Callback.H \
	$(srcdir)/Callback.C \
	$(srcdir)/PyInterface.H \
	$(srcdir)/PyInterface.C

_Epetra_la_CXXFLAGS = $(SWIG_NOX_CXXFLAGS) $(SWIG_EPETRA_CXXFLAGS)

_Epetra_la_LDFLAGS = $(SWIG_LDFLAGS)

_Epetra_la_LIBADD = $(SWIG_EPETRA_LIBADD) $(SWIG_NOX_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _LAPACK.la
# ------------------------------------------------------------------------

_LAPACK_la_SOURCES = $(srcdir)/wrap_NOX_LAPACK.C 

_LAPACK_la_CXXFLAGS = $(SWIG_NOX_CXXFLAGS) $(SWIG_LAPACK_CXXFLAGS)

_LAPACK_la_LDFLAGS = $(SWIG_LDFLAGS)

_LAPACK_la_LIBADD = $(SWIG_LAPACK_LIBADD) $(SWIG_NOX_LIBADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# Include swig rules
include $(top_srcdir)/config/swig.am






