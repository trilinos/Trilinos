#np# Typically the only line that needs to be edited in this file is the last 
#np# line.  Set 'SUBDIRS= (list all subdirectories that are autotool'ed here)'
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.

# @HEADER
# ************************************************************************
# 
#                     New_Package Example Package
#                 Copyright (2004) Sandia Corporation
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
#EXTRA_DIST = \
#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
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

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################
#  #np# - Most packages should have an example and a full suite 
#  of tests.  At present, new_package only offers src and example.
#  BUILD_TESTS and BUILD_EXAMPLES are defined through the configure.ac file.

SWIG       = swig
SWIG_FLAGS = -I$(top_srcdir)/src/NOX -Wall -c -python -c++

EPETRA_INC    = -I$(top_srcdir)/../epetra/src
NOX_INC       = -I$(top_srcdir)/../nox/src
NOXEPETRA_INC = -I$(top_srcdir)/../nox/src-epetra

INTERFACES = Abstract.i Epetra.i Parameter.i Solver.i StatusTest.i
WRAPPERS   = $(patsubst %.i, %_wrap.cxx, $(INTERFACES))

TARGET     = $(WRAPPERS)

all-am: $(TARGET)

# Generate the Abstract C++ wrapper and proxy file from the Abstract SWIG interface
Abstract_wrap.cxx Abstract.py: Abstract.i
	$(SWIG) $(NOX_INC) $(SWIG_FLAGS) $<

# Generate the Epetra C++ wrapper and proxy file from the Epetra SWIG interface
Epetra_wrap.cxx Epetra.py: Epetra.i
	$(SWIG) $(EPETRA_INC) $(NOX_INC) $(NOXEPETRA_INC) $(SWIG_FLAGS) $<

# Generate the Parameter C++ wrapper and proxy file from the Parameter SWIG interface
Parameter_wrap.cxx Parameter.py: Parameter.i
	$(SWIG) $(NOX_INC) $(SWIG_FLAGS) $<

# Generate the Solver C++ wrapper and proxy file from the Solver SWIG interface
Solver_wrap.cxx Solver.py: Solver.i
	$(SWIG) $(NOX_INC) $(SWIG_FLAGS) $<

# Generate the StatusTest C++ wrapper and proxy file from the StatusTest SWIG interface
StatusTest_wrap.cxx StatusTest.py: StatusTest.i
	$(SWIG) $(NOX_INC) $(SWIG_FLAGS) $<
