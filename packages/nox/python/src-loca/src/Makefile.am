# @HEADER
# ************************************************************************
#
#                PyTrilinos.LOCA: Python Interface to LOCA
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

# Include swig rules
include $(top_srcdir)/config/swig.am
SWIG_CXXFLAGS = \
	-I$(top_srcdir)/src-loca/src        \
	-I$(top_srcdir)/src-loca/src-lapack \
	-I$(top_srcdir)/src-loca/src-epetra \
	-I$(top_srcdir)/src                 \
	-I$(top_srcdir)/src-lapack          \
	-I$(top_srcdir)/src-epetra          \
	-I$(top_srcdir)/python/src          \
	-I$(top_srcdir)/python/src-loca/src

# Default make target
all-am: build

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(PYTHON) -c "import sys,os; from distutils.util import get_platform; print os.path.join('build','temp.%s-%s'%(get_platform(),sys.version[0:3]))")

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@

# Set the appropriate Make macros
if BUILD_LIBLOCAEPETRA
  EPETRA_I=LOCA_Epetra.i
  LOCA_EPETRA=1
else
  LOCA_EPETRA=
endif

if BUILD_LIBLOCALAPACK
  LAPACK_I=LOCA_LAPACK.i
#  CHAN_I=LOCA_Chan.i
  LOCA_LAPACK=1
else
  LOCA_LAPACK=
endif

# Define the SWIG interface files and the resulting wrapper files
INTERFACES = LOCA_TopLevel.i LOCA_Continuation.i LOCA_MultiContinuation.i \
	     LOCA_Homotopy.i LOCA_TimeDependent.i LOCA_Bifurcation.i      \
	     LOCA_Abstract.i $(EPETRA_I) $(LAPACK_I)
#	     $(CHAN_I)
WRAPPERS    = $(INTERFACES:.i=_wrap.cpp)
PROXIES_TMP = $(INTERFACES:.i=.py      )
PROXIES     = $(PROXIES_TMP:LOCA_=     )

# File setup.txt is a way to transfer make variables to the setup.py
# python script
setup.txt: Makefile
	@echo "creating $@"
	@echo "{"                                   >  $@
	@echo "  'top_srcdir' : '$(top_srcdir)',"   >> $@
	@echo "  'srcdir'     : '$(srcdir)',"       >> $@
	@echo "  'LDFLAGS'    : '$(LDFLAGS)',"      >> $@
	@echo "  'LAPACK_LIBS': '$(LAPACK_LIBS)',"  >> $@
	@echo "  'BLAS_LIBS'  : '$(BLAS_LIBS)',"    >> $@
	@echo "  'FLIBS'      : '$(FLIBS)',"        >> $@
	@echo "  'LIBS'       : '$(LIBS)',"         >> $@
	@echo "  'PYTHON'     : '$(PYTHON)',"       >> $@
	@echo "  'LOCA_EPETRA': '$(LOCA_EPETRA)',"  >> $@
	@echo "  'LOCA_LAPACK': '$(LOCA_LAPACK)',"  >> $@
	@echo "}"                                   >> $@

build: $(WRAPPERS) setup.txt $(tempdir)
	$(PYTHON) $(srcdir)/setup.py build

install: build
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(PYTHON_PREFIX)

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean
	$(RM) -r ../../packages

CLEANFILES += setup.txt $(WRAPPERS) $(PROXIES)
