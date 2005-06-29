# @HEADER
# ************************************************************************
# 
#               Epetra: Linear Algebra Services Package 
#                 Copyright (2001) Sandia Corporation
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
#EXTRA_DIST =

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

#The following line helps the test harness recover from build errors.

all-local:
	@echo "Trilinos package epetra built successfully"

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

if BUILD_TESTS
TEST_SUBDIR=test
endif

if BUILD_EXAMPLES
EXAMPLE_SUBDIR=example
endif

if PYTHON_USE
PYTHON_SUBDIR=python
endif

if ENABLE_EPETRA_THYRA
THYRA_SUBDIR=thyra
else
THYRA_SUBDIR=
endif

SUBDIRS = src $(EXAMPLE_SUBDIR) $(TEST_SUBDIR) $(PYTHON_SUBDIR) $(THYRA_SUBDIR)

if ENABLE_EPETRA_THYRA
check:
	$(PERL_EXE) $(top_builddir)/test/scripts/daily/serial/run_epetra_thyra_tests.pl
else
check:
endif

if USING_EXPORT_MAKEFILES

install-exec-hook:	
	cp $(top_builddir)/Makefile.export.epetra $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.epetra \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_srcdir@/thyra/src:@abs_top_builddir@/src:@abs_top_srcdir@/src \
		--my-abs-lib-dirs=@abs_top_builddir@/thyra/src:@abs_top_builddir@/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../thyra
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile EPETRA > $(exec_prefix)/include/Makefile.export.epetra.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.thyra
	rm -f $(exec_prefix)/include/Makefile.export.thyra.macros

else 

install-exec-hook:

uninstall-hook:

endif
