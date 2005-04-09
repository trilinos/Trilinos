# @HEADER
# ************************************************************************
# 
#                  PyTrilinos: Rapid Prototyping Package
#                   Copyright (2005) Sandia Corporation
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
# @HEADER

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
# We will list some source files here because we don't compile them in the
# usual way.
#EXTRA_DIST = \
#	Epetra_VectorHelper.h  NumPyArray.cxx \
#	__init__.py            NumPyArray.h \
#	EpetraExt.i            NumPyWrapper.cxx \
#	Epetra_NumPyVector.cxx NumPyWrapper.h \
#	Epetra_NumPyVector.h \
#	Epetra.py              PyTrilinos.h \
#	Epetra_VectorHelper.cxx \
#	numeric_include.h      RawEpetra.i

## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
#AUX_DIST = config/install-sh config/missing config/mkinstalldirs 
#
#  Again, I hope that AUX_DIR(../../config) eliminates these
#  config/install-sh config/missing config/mkinstalldirs 

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in

# Include paths
EPETRA_INC    = -I$(top_srcdir)/../epetra/src
EPETRAEXT_INC = -I$(top_srcdir)/../epetraext/src/transform \
                -I$(top_srcdir)/../epetraext/src/coloring
NOX_INC       = -I$(top_srcdir)/../nox/src
NOXEPETRA_INC = -I$(top_srcdir)/../nox/src-epetra

# Include swig rules
include $(top_srcdir)/config/swig.am
SWIG_CXXFLAGS = \
	$(EPETRA_INC)    \
	$(EPETRAEXT_INC) \
	$(NOX_INC)       \
	$(NOXEPETRA_INC) \
	-I$(srcdir)

# Default make target
all-am: build

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(PYTHON) $(top_srcdir)/pyLocate --temp)

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@

# Define the SWIG interface files and the resulting wrapper files
INTERFACES = Epetra.i EpetraExt.i NOX/Abstract.i NOX/Epetra.i NOX/Parameter.i \
	     NOX/Solver.i NOX/StatusTest.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cxx)
PROXIES    = $(INTERFACES:.i=.py      )

# File setup.txt is a way to transfer make variables to the setup.py
# python script
setup.txt: $(srcdir)/Makefile.am
	@echo "creating $@"
	@echo "{"                                   >  $@
	@echo "  'top_srcdir'  : '$(top_srcdir)',"  >> $@
	@echo "  'srcdir'      : '$(srcdir)',"      >> $@
	@echo "  'LDFLAGS'     : '$(LDFLAGS)',"     >> $@
	@echo "  'LAPACK_LIBS' : '$(LAPACK_LIBS)'," >> $@
	@echo "  'BLAS_LIBS'   : '$(BLAS_LIBS)',"   >> $@
	@echo "  'FLIBS'       : '$(FLIBS)',"       >> $@
	@echo "  'LIBS'        : '$(LIBS)',"        >> $@
	@echo "  'PYTHON'      : '$(PYTHON)'"       >> $@
	@echo "}"                                   >> $@

build: $(WRAPPERS) setup.txt $(tempdir)
	$(PYTHON) $(srcdir)/setup.py build

install: build
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(PYTHON_PREFIX)

CLEANFILES += setup.txt $(WRAPPERS) $(PROXIES)

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean
	$(RM) -r ../../packages
