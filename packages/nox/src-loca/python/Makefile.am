## $Id$ 
## $Source$ 
# @HEADER
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
#  CVS Information
#  $Source$
#  $Author$
#  $Date$
#  $Revision$
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Swig interface files
# ------------------------------------------------------------------------

LOCA_INTERFACES = \
	$(srcdir)/LOCA_TopLevel.i \
	$(srcdir)/LOCA_Continuation.i \
	$(srcdir)/LOCA_MultiContinuation.i \
	$(srcdir)/LOCA_Homotopy.i \
	$(srcdir)/LOCA_TimeDependent.i \
	$(srcdir)/LOCA_Bifurcation.i \
	$(srcdir)/LOCA_Abstract.i 

if BUILD_LIBLOCALAPACK

LAPACK_INTERFACES = \
	$(srcdir)/LOCA_LAPACK.i \
	$(srcdir)/LOCA_Chan.i

else

LAPACK_INTERFACES =

endif

INTERFACES = $(LOCA_INTERFACES) $(LAPACK_INTERFACES)

# ------------------------------------------------------------------------
# Generated swig wrapper files
# ------------------------------------------------------------------------

LOCA_WRAPPERS = \
	$(srcdir)/TopLevel.py \
	$(srcdir)/Continuation.py \
	$(srcdir)/MultiContinuation.py \
	$(srcdir)/Homotopy.py \
	$(srcdir)/TimeDependent.py \
	$(srcdir)/Bifurcation.py \
	$(srcdir)/Abstract.py 

if BUILD_LIBLOCALAPACK

LAPACK_WRAPPERS = \
	$(srcdir)/LAPACK.py \
	$(srcdir)/Chan.py 

else

LAPACK_WRAPPERS =

endif

WRAPPERS = $(LOCA_WRAPPERS) $(LAPACK_WRAPPERS)

# ------------------------------------------------------------------------
# Python extension modules as libtool shared libraries
# ------------------------------------------------------------------------

LOCA_EXTENSIONS = \
	_TopLevel.la \
	_Continuation.la \
	_MultiContinuation.la \
	_Homotopy.la \
	_TimeDependent.la \
	_Bifurcation.la \
	_Abstract.la 

if BUILD_LIBLOCALAPACK

LAPACK_EXTENSIONS = \
	_LAPACK.la \
	_Chan.la

else

LAPACK_EXTENSIONS =

endif

EXTENSIONS = $(LOCA_EXTENSIONS) $(LAPACK_EXTENSIONS)

# ------------------------------------------------------------------------
# Where we want the LOCA/Python installed.  We use 'exec' in the name so 
# that it is installed by 'install-exec'
# ------------------------------------------------------------------------

locaPythonexecdir = $(prefix)/python/LOCA

if BUILD_PYTHON

# ------------------------------------------------------------------------
# What to build
# ------------------------------------------------------------------------

locaPythonexec_LTLIBRARIES = $(EXTENSIONS)

# ------------------------------------------------------------------------
# Python files to be installed by 'make install-exec'
# ------------------------------------------------------------------------

# __init.py__ is not distributed because it is created by configure
nodist_locaPythonexec_PYTHON = __init__.py

dist_locaPythonexec_PYTHON = $(WRAPPERS) 

endif

# ------------------------------------------------------------------------
# Extra files to be installed by 'make dist'
# ------------------------------------------------------------------------

EXTRA_DIST = $(INTERFACES)

# ------------------------------------------------------------------------
# General definitions
# ------------------------------------------------------------------------

SWIG_PYTHON_LIBS =

SWIG_LOCA_CXXFLAGS = \
	$(PYTHON_CPPFLAGS) \
	$(SWIG_PYTHON_CPPFLAGS) \
	$(TRILINOS_CXXFLAGS) \
	$(LIBTOOL_COMPILER_FLAGS) \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/python \
	-I$(top_srcdir)/src-loca/src \
	-I$(top_srcdir)/src-loca/python 

SWIG_LAPACK_CXXFLAGS = \
	-I$(top_srcdir)/src-lapack \
	-I$(top_srcdir)/src-loca/src-lapack \
	-I$(top_srcdir)/src-loca/example-lapack/Chan 

SWIG_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS) $(SWIG_LAPACK_CXXFLAGS)

SWIG_LDFLAGS = -avoid-version -module

SWIG_LOCA_LIBADD = \
	$(top_builddir)/src-loca/src/libloca.la \
	$(top_builddir)/src/libnox.a \
	$(CXXREPOS) $(XTRALDADD) $(LIBTOOL_CCLINKER_FLAGS) \
	$(LIBS) $(SWIG_PYTHON_LIBS)

SWIG_LAPACK_LIBADD = \
	$(top_builddir)/src-loca/src-lapack/liblocalapack.la \
	$(top_builddir)/src-lapack/libnoxlapack.a

# ------------------------------------------------------------------------
# Definitions for _TopLevel.la
# ------------------------------------------------------------------------

_TopLevel_la_SOURCES = \
	$(srcdir)/wrap_LOCA_TopLevel.C \
	$(srcdir)/Utils_enums.H

_TopLevel_la_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS)

_TopLevel_la_LDFLAGS = $(SWIG_LDFLAGS)

_TopLevel_la_LIBADD = $(SWIG_LOCA_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _Contination.la
# ------------------------------------------------------------------------

_Continuation_la_SOURCES = $(srcdir)/wrap_LOCA_Continuation.C

_Continuation_la_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS)

_Continuation_la_LDFLAGS = $(SWIG_LDFLAGS)

_Continuation_la_LIBADD = $(SWIG_LOCA_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _MultiContination.la
# ------------------------------------------------------------------------

_MultiContinuation_la_SOURCES = $(srcdir)/wrap_LOCA_MultiContinuation.C

_MultiContinuation_la_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS)

_MultiContinuation_la_LDFLAGS = $(SWIG_LDFLAGS)

_MultiContinuation_la_LIBADD = $(SWIG_LOCA_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _Homotopy.la
# ------------------------------------------------------------------------

_Homotopy_la_SOURCES = $(srcdir)/wrap_LOCA_Homotopy.C

_Homotopy_la_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS)

_Homotopy_la_LDFLAGS = $(SWIG_LDFLAGS)

_Homotopy_la_LIBADD = $(SWIG_LOCA_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _TimeDependent.la
# ------------------------------------------------------------------------

_TimeDependent_la_SOURCES = $(srcdir)/wrap_LOCA_TimeDependent.C

_TimeDependent_la_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS)

_TimeDependent_la_LDFLAGS = $(SWIG_LDFLAGS)

_TimeDependent_la_LIBADD = $(SWIG_LOCA_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _Bifurcation.la
# ------------------------------------------------------------------------

_Bifurcation_la_SOURCES = $(srcdir)/wrap_LOCA_Bifurcation.C

_Bifurcation_la_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS)

_Bifurcation_la_LDFLAGS = $(SWIG_LDFLAGS)

_Bifurcation_la_LIBADD = $(SWIG_LOCA_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _Abstract.la
# ------------------------------------------------------------------------

_Abstract_la_SOURCES = $(srcdir)/wrap_LOCA_Abstract.C

_Abstract_la_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS)

_Abstract_la_LDFLAGS = $(SWIG_LDFLAGS)

_Abstract_la_LIBADD = $(SWIG_LOCA_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _LAPACK.la
# ------------------------------------------------------------------------

_LAPACK_la_SOURCES = $(srcdir)/wrap_LOCA_LAPACK.C

_LAPACK_la_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS) $(SWIG_LAPACK_CXXFLAGS)

_LAPACK_la_LDFLAGS = $(SWIG_LDFLAGS)

_LAPACK_la_LIBADD = $(SWIG_LAPACK_LIBADD) $(SWIG_LOCA_LIBADD)

# ------------------------------------------------------------------------
# Definitions for _Chan.la
# ------------------------------------------------------------------------

_Chan_la_SOURCES = \
	$(srcdir)/wrap_LOCA_Chan.C \
	$(top_srcdir)/src-loca/example-lapack/Chan/ChanProblemInterface.H \
	$(top_srcdir)/src-loca/example-lapack/Chan/ChanProblemInterface.C 

_Chan_la_CXXFLAGS = $(SWIG_LOCA_CXXFLAGS) $(SWIG_LAPACK_CXXFLAGS)

_Chan_la_LDFLAGS = $(SWIG_LDFLAGS)

_Chan_la_LIBADD = $(SWIG_LAPACK_LIBADD) $(SWIG_LOCA_LIBADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

# ------------------------------------------------------------------------
# Include global libtool definitions
# ------------------------------------------------------------------------
include $(top_srcdir)/config/libtool_macros.am

# Include swig rules
include $(top_srcdir)/config/swig.am






