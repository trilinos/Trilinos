# @HEADER
# ************************************************************************
#
# The copyright to the AMD package belongs to Tim Davis (davis@cse.ufl.edu) 
# Please refer to the copyright notice in src/Makefile.am for details.
#   
# 
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
## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
#EXTRA_DIST = 

#
## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST =
#
## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in $(AUX_DIST) 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

#if BUILD_TESTS
#TEST_SUBDIR=test
#endif

#if BUILD_EXAMPLES
EXAMPLE_SUBDIR=example
#endif

#SUBDIRS = src $(EXAMPLE_SUBDIR) $(TEST_SUBDIR)
SUBDIRS = src $(EXAMPLE_SUBDIR)

