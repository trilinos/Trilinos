# @HEADER
# ************************************************************************
# 
#              PyTrilinos.Epetra: Python Interface to Epetra
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
SWIG_CXXFLAGS = -I$(top_srcdir)/src -I$(srcdir)

# Default make target
all-am: build

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(top_srcdir)/../PyTrilinos/pyLocate --temp)

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@

# Define the SWIG interface files and the resulting wrapper files
INTERFACES = Epetra.i
WRAPPERS   = $(patsubst %.i, $(srcdir)/wrap_%.cpp, $(INTERFACES))

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
	@echo "}"                                   >> $@

build: $(WRAPPERS) setup.txt $(tempdir)
	$(srcdir)/setup.py build

install: build
	$(srcdir)/setup.py install

clean:
	$(srcdir)/setup.py clean
	rm -f setup.txt
	rm -f $(WRAPPERS)

test:
	cd test && make test
