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
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Aclocal command (rarely used - don't worry about it)
## #######################################################################
ACLOCAL	= aclocal -I config

#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = \
	Abstract.i    __init__.py  PyInterface.h \
	Callback.cxx  Parameter.i  Solver.i \
	Callback.h    Epetra.i     PyInterface.cxx  StatusTest.i



## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs 
#
#  Again, I hope that AUX_DIR(../../config) eliminates these
#  config/install-sh config/missing config/mkinstalldirs 

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in \
	$(AUX_DIST) 

# SWIG       = swig
# SWIG_FLAGS = -I$(top_srcdir)/src/NOX -Wall -python -c++

EPETRA_INC    = -I$(top_srcdir)/../epetra/src
NOX_INC       = -I$(top_srcdir)/../nox/src
NOXEPETRA_INC = -I$(top_srcdir)/../nox/src-epetra

TARGET     = $(WRAPPERS)

all-am: $(TARGET)

# Include swig rules
include $(top_srcdir)/config/swig.am
SWIG_CXXFLAGS = $(EPETRA_INC) $(NOX_INC) $(NOXEPETRA_INC) -I$(srcdir)

# Default make target
all-am: build

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(PYTHON) $(top_srcdir)/pyLocate --temp)

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@

# Define the SWIG interface files and the resulting wrapper files
INTERFACES = Abstract.i Epetra.i Parameter.i Solver.i StatusTest.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cxx)
PROXIES    = $(INTERFACES:.i=.py      )

build: $(WRAPPERS) $(tempdir)
	$(PYTHON) $(srcdir)/setup.py build

install: build
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(PYTHON_PREFIX)

clean:
	$(PYTHON) $(srcdir)/setup.py clean
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
