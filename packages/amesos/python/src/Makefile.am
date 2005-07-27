# @HEADER
# ************************************************************************
# 
#              PyTrilinos.Amesos : Python Interface to Amesos 
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
include $(top_builddir)/Makefile.export.amesos

if USING_GNUMAKE
PYTHON_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(AMESOS_LIBS) $(AMESOS_EXTRA_LIBS))
PYTHON_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES))
else
PYTHON_LIBS     = $(AMESOS_LIBS)     $(AMESOS_EXTRA_LIBS)
PYTHON_INCLUDES = $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES)
endif

SWIG_CXXFLAGS = $(PYTHON_INCLUDES) \
		-I$(top_srcdir)/../epetra/python/src

EXTRA_DIST = Amesos.i setup.py

# Default make target
all-am: build

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(PYTHON) $(top_srcdir)/../PyTrilinos/pyLocate --temp)

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@

# Define the SWIG interface files and the resulting wrapper files
INTERFACES = Amesos.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cpp)
PROXIES    = $(INTERFACES:.i=.py      )

# File setup.txt is a way to transfer make variables to the setup.py
# python script
setup.txt: $(srcdir)/Makefile.am
	@echo "creating $@"
	@echo "{"                                   >  $@
	@echo "  'top_builddir' : '$(top_builddir)',"  >> $@
	@echo "  'top_srcdir'   : '$(top_srcdir)',"  >> $@
	@echo "  'srcdir'       : '$(srcdir)',"      >> $@
	@echo "  'CC'           : '$(CC)',"          >> $@
	@echo "  'CXX'          : '$(CXX)',"         >> $@
	@echo "  'CXXFLAGS'     : '$(CXXFLAGS)',"    >> $@
	@echo "  'LDFLAGS'      : '$(PYTHON_LIBS)'," >> $@
	@echo "  'PYTHON'       : '$(PYTHON)'"       >> $@
	@echo "}"                                    >> $@

build: $(WRAPPERS) setup.txt $(tempdir)
	$(PYTHON) $(srcdir)/setup.py build

install: build
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(PYTHON_PREFIX)

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean
	$(RM) -r ../../packages

CLEANFILES += setup.txt $(WRAPPERS) $(PROXIES)
