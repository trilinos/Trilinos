#np# Typically the only line that needs to be edited in this file is the last 
#np# line.  Set 'SUBDIRS= (list all subdirectories that are autotool'ed here)'
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.

# @HEADER
# ************************************************************************
# 
# Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms
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
ACLOCAL_AMFLAGS = -I config

#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = bootstrap \
  config/ac_cxx_have_std.m4 config/ac_cxx_have_stl.m4 \
  config/ac_cxx_namespaces.m4 config/ac_cxx_mutable.m4 \
  config/depcomp config/generate-makeoptions.pl \
  config/strip_dup_incl_paths.pl config/strip_dup_libs.pl \
  config/replace-install-prefix.pl config/string-replace.pl \
  config/token-replace.pl


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
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
  configure config.status config.log \
  src/common/config-h.in src/common/stamp-h.in \
  $(AUX_DIST) 

#The following line helps the test harness recover from build errors.

all-local:
	@echo "Trilinos package thyra built successfully"

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though BUILD_TESTS and BUILD_EXAMPLES will never be
#defined, so that the tests and examples are included in the distribution
#tarball.

if SUB_TEST
TEST_SUBDIR=test
endif

if SUB_EXAMPLE
EXAMPLE_SUBDIR=example
endif

if ENABLE_THYRA_EPETRA
BUILD_THYRA_EPETRA_TESTS=cd adapters/epetra/test && $(MAKE)
BUILD_THYRA_EPETRA_EXAMPLES=cd adapters/epetra/example && $(MAKE)
BUILD_THYRA_EPETRA_CLEAN_CMD=clean
BUILD_THYRA_EPETRA_INSTALL_CMD=install
endif

if ENABLE_THYRA_EPETRAEXT
BUILD_THYRA_EPETRAEXT_EXAMPLES=cd adapters/epetraext/example && $(MAKE)
BUILD_THYRA_EPETRAEXT_CLEAN_CMD=clean
BUILD_THYRA_EPETRAEXT_INSTALL_CMD=install
endif

if BUILD_TESTS
tests: examples
	@echo ""
	@echo "Now building thyra tests."
	@echo ""
	cd $(top_builddir)/test && $(MAKE)
	$(BUILD_THYRA_EPETRA_TESTS)
	@echo ""
	@echo "Finished building thyra tests."
	@echo ""
else
tests:
	@echo "thyra tests were disabled at configure time"
endif

if BUILD_EXAMPLES
examples:
	@echo ""
	@echo "Now building thyra examples."
	@echo ""
	cd $(top_builddir)/example && $(MAKE)
	$(BUILD_THYRA_EPETRA_EXAMPLES)
	$(BUILD_THYRA_EPETRAEXT_EXAMPLES)
	@echo ""
	@echo "Finished building thyra examples."
	@echo ""
install-examples:
	cd $(top_builddir)/example && $(MAKE) install
	$(BUILD_THYRA_EPETRA_EXAMPLES) $(BUILD_THYRA_EPETRA_INSTALL_CMD)
	$(BUILD_THYRA_EPETRAEXT_EXAMPLES) $(BUILD_THYRA_EPETRAEXT_INSTALL_CMD)
else
examples:
	@echo "thyra examples were disabled at configure time"

install-examples:
	@echo "thyra examples were disabled at configure time"
endif

clean-tests:
	cd $(top_builddir)/test && $(MAKE) clean
	$(BUILD_THYRA_EPETRA_TESTS) $(BUILD_THYRA_EPETRA_CLEAN_CMD)

clean-examples:
	cd $(top_builddir)/example && $(MAKE) clean
	$(BUILD_THYRA_EPETRA_EXAMPLES) $(BUILD_THYRA_EPETRA_CLEAN_CMD)
	$(BUILD_THYRA_EPETRAEXT_EXAMPLES) $(BUILD_THYRA_EPETRAEXT_CLEAN_CMD)

everything:
	$(MAKE) && $(MAKE) examples && $(MAKE) tests

clean-everything:
	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

install-everything:
	$(MAKE) install && $(MAKE) install-examples

SUBDIRS = src adapters $(EXAMPLE_SUBDIR) $(TEST_SUBDIR)

## #######################################################################
## The below targets allow you to use the new
## testharness to run the test suite as make targets
## #######################################################################

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..
TRILINOS_MPI_MAX_PROC=4
TRILINOS_TEST_CATEGORY=INSTALL

runtests-serial :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=serial \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --category=$(TRILINOS_TEST_CATEGORY) \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=thyra

runtests-mpi :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=mpi \
  --mpi-go=$(TRILINOS_MPI_GO) \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --max-proc=$(TRILINOS_MPI_MAX_PROC) \
  --category=$(TRILINOS_TEST_CATEGORY) \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=thyra

if HAVE_MPI
THYRA_CHECK_COMM=mpi
else
THYRA_CHECK_COMM=serial
endif

## #######################################################################
## Export makefile stuff
## #######################################################################

sample-output:
	cd test ; ../example/generate-sample-output.pl
	cd example/external ; $(MAKE) sample-makefile
	cd $(top_builddir)/adapters ; $(MAKE) sample-output


if USING_EXPORT_MAKEFILES

install-exec-hook:	
	mkdir -p $(DESTDIR)$(includedir)
	cp $(top_builddir)/Makefile.export.thyra $(DESTDIR)$(includedir)/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.thyra \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src:@abs_top_srcdir@/adapters/epetra/src:@abs_top_srcdir@/adapters/tpetra/src:@abs_top_srcdir@/adapters/epetraext/src/model_evaluator \
		--my-abs-lib-dirs=@abs_top_builddir@/src:@abs_top_builddir@/adapters/epetra/src:@abs_top_builddir@/adapters/tpetra/src:@abs_top_builddir@/adapters/epetraext/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../rtop:@abs_top_builddir@/../epetra:@abs_top_builddir@/../tpetra:@abs_top_builddir@/../epetraext
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile THYRA > $(DESTDIR)$(includedir)/Makefile.export.thyra.macros

uninstall-hook:
	rm -f $(includedir)/Makefile.export.thyra
	rm -f $(includedir)/Makefile.export.thyra.macros

else 

install-exec-hook:

uninstall-hook:

endif
