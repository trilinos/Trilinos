# @HEADER
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
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

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = bootstrap \
	config/ac_cxx_have_std.m4 config/ac_cxx_have_stl.m4 \
	config/ac_cxx_namespaces.m4 config/ac_cxx_mutable.m4 \
	config/depcomp

## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs 

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	src/common/config-h.in src/common/stamp-h.in \
	$(AUX_DIST) 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################
#SUBDIRS = src
#SUBDIRS = amd src test example 
#SUBDIRS = src test


## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

if BUILD_TESTS
TEST_SUBDIR=test
endif

if BUILD_EXAMPLES
EXAMPLE_SUBDIR=example
endif

SUBDIRS = src $(EXAMPLE_SUBDIR) $(TEST_SUBDIR)

