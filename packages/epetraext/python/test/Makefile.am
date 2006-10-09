# @HEADER
# ************************************************************************
# 
#           PyTrilinos.EpetraExt: Python Interface to EpetraExt
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
	setpath.py testColoring.py

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in

# Define the suite of tests
TEST_SUITE = testColoring

SETPATH = setpath.py

all: $(TEST_SUITE)
	if [ `$(PYTHON) -c "import os; print os.path.abspath('$(srcdir)')"` != `pwd` ]; then \
	   /bin/cp $(srcdir)/$(SETPATH) . ; fi

test%:
	if [ `$(PYTHON) -c "import os; print os.path.abspath('$(srcdir)')"` != `pwd` ]; then \
	   /bin/cp $(srcdir)/$@.py . ; fi

CLEANFILES = $(TEST_SUITE:=.py) $(SETPATH) *.pyc
