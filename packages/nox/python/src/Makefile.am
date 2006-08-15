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

# Include export Makefiles
include $(top_builddir)/Makefile.export.nox

# Include swig rules
include $(top_srcdir)/config/swig.am
SWIG_CXXFLAGS = $(NOX_INCLUDES) \
		-I$(top_srcdir)/../epetra/python/src

# Default make target
all-am: build

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
if BUILD_NOX_PYTHON
tempdir = $(shell $(PYTHON) -c "import sys,os; from distutils.util import get_platform; print os.path.join('build','temp.%s-%s'%(get_platform(),sys.version[0:3]))")

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@
endif

# Set the appropriate Make macros
if BUILD_NOX_PYTHON
if BUILD_LIBNOXEPETRA
  EPETRA_I=NOX_Epetra.i
  PYTHON_NOX_EPETRA=1
else
  EPETRA_I=
  PYTHON_NOX_EPETRA=0
endif
if BUILD_LIBNOXLAPACK
  LAPACK_I=NOX_LAPACK.i
  PYTHON_NOX_LAPACK=1
else
  LAPACK_I= 
  PYTHON_NOX_LAPACK=0
endif
endif

# Define the SWIG interface files and the resulting wrapper files
INTERFACES  = NOX_TopLevel.i NOX_Abstract.i NOX_Parameter.i NOX_Solver.i \
	      NOX_StatusTest.i $(EPETRA_I) $(LAPACK_I)
WRAPPERS    = $(INTERFACES:.i=_wrap.cpp)
PROXIES_TMP = $(INTERFACES:.i=.py      )
PROXIES     = $(PROXIES_TMP:NOX_=      )

if BUILD_NOX_PYTHON
build: $(WRAPPERS) $(tempdir)
	$(PYTHON) $(srcdir)/setup.py build

install: build
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(PYTHON_PREFIX)

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean
	$(RM) -r ../../packages
else
build:
	echo 'NOX-Python support is disabled!'
endif

CLEANFILES += $(WRAPPERS) $(PROXIES)

# ------------------------------------------------------------------------
# Extra files to be installed by 'make dist'
# ------------------------------------------------------------------------

EXTRA_DIST = \
	$(INTERFACES) \
	setup.py \
	Epetra_VectorHelper.h \
	Epetra_VectorHelper.cpp \
	PyInterface.h \
	PyInterface.cpp \
	NOX_Epetra.i \
	NOX_LAPACK.i \
	Utils_enums.H

