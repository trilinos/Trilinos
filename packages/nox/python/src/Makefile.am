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
# @HEADER

# Include export Makefiles
include $(top_builddir)/Makefile.export.nox

# Include swig rules
include $(top_srcdir)/config/swig.am
SWIG_CXXFLAGS = $(NOX_INCLUDES) \
		-I$(top_srcdir)/../teuchos/python/src \
		-I$(top_srcdir)/../epetra/python/src

# Default make target
all-am: build

if BUILD_NOX_PYTHON

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(PYTHON) -c "import sys,os; from distutils.util import get_platform; print os.path.join('build','temp.%s-%s'%(get_platform(),sys.version[0:3]))")

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@

# Top-level NOX interface, wrapper and proxy variables
NOX_INTERFACES = NOX.__init__.i   \
	         NOX.Abstract.i   \
	         NOX.Solver.i     \
	         NOX.StatusTest.i
NOX_WRAPPERS   = $(NOX_INTERFACES:.i=_wrap.cpp)
NOX_TEMPS      = $(NOX_INTERFACES:NOX.=PyTrilinos/NOX/)
NOX_PROXIES    = $(NOX_PROXIES_TMP:.i=.py)

# NOX.Epetra interface, wrapper and proxy variables
if BUILD_LIBNOXEPETRA
  NOX_EPETRA_INTERFACES = NOX.Epetra.__init__.i \
	                  NOX.Epetra.Interface.i
  NOX_EPETRA_WRAPPERS   = $(NOX_EPETRA_INTERFACES:.i=_wrap.cpp)
  NOX_EPETRA_TEMPS      = $(NOX_EPETRA_INTERFCES:NOX.Epetra.=PyTrilinos/NOX/Epetra/)
  NOX_EPETRA_PROXIES    = $(NOX_EPETRA_TEMPS:.i=.py)
  PYTHON_NOX_EPETRA=1
else
  PYTHON_NOX_EPETRA=0
endif

INTERFACES = $(NOX_INTERFACES) $(NOX_EPETRA_INTERFACES)
WRAPPERS   = $(NOX_WRAPPERS)   $(NOX_EPETRA_WRAPPERS)
PROXIES    = $(NOX_PROXIES)    $(NOX_EPETRA_PROXIES)

build: $(WRAPPERS) $(tempdir) # $(EPETRA_DEP)
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

CLEANFILES += $(WRAPPERS) $(PROXIES) $(EPETRA_DEP)

# ------------------------------------------------------------------------
# Extra files to be installed by 'make dist'
# ------------------------------------------------------------------------

#EXTRA_DIST = \
#	$(INTERFACES) \
#	setup.py \
#	NOX_TopLevel.i \
#	NOX_Abstract.i \
#	NOX_Solver.i \
#	NOX_StatusTest.i \
#	NOX_Epetra.i \
#	NOX_Epetra_Interface.i
