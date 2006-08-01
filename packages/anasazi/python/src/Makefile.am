# @HEADER
# ************************************************************************
#
#             PyTrilinos.Anasazi: Python Interface to Anasazi
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

# Include the export Makefiles
include $(top_builddir)/Makefile.export.anasazi

# MPI support
@HAVE_MPI_TRUE@MPI_FLAG = -DHAVE_MPI

# Include swig rules
include $(top_srcdir)/config/swig.am
SWIG_CXXFLAGS = $(ANASAZI_INCLUDES) \
		$(MPI_FLAG)

# Default make target
all-am: build

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(PYTHON) -c "import sys,os; from distutils.util import get_platform; print os.path.join('build','temp.%s-%s'%(get_platform(),sys.version[0:3]))")

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@

# Define the SWIG interface files and the resulting wrapper files
INTERFACES = Anasazi.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cpp)
PROXIES    = $(INTERFACES:.i=.py      )

EXTRA_DIST = Anasazi.i setup.py

build: $(WRAPPERS) $(tempdir)
	$(PYTHON) $(srcdir)/setup.py build

install: build
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(PYTHON_PREFIX)

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean
	$(RM) -r ../../packages

CLEANFILES += $(WRAPPERS) $(PROXIES)
