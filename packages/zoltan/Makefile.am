# @HEADER
##############################################################################
# Zoltan Library for Parallel Applications                                   #
# Copyright (c) 2000,2001,2002, Sandia National Laboratories.                #
# This software is distributed under the GNU Lesser General Public License.  #
# For more info, see the README file in the top-level Zoltan directory.      #
##############################################################################
# @HEADER


## #######################################################################
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Aclocal command (rarely used - don't worry about it)
## #######################################################################
ACLOCAL_AMFLAGS = -I config

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = config/generate-makeoptions.pl \
       config/strip_dup_incl_paths.pl config/strip_dup_libs.pl \
       config/replace-install-prefix.pl config/string-replace.pl \
       config/token-replace.pl

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
	@echo ""
	@echo "Trilinos package zoltan built successfully."
	@echo ""

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though BUILD_TESTS and BUILD_EXAMPLES will never be
#defined, so that the tests and examples are included in the distribution
#tarball.

#Add this later
if SUB_TEST
TEST_SUBDIR=
endif

if SUB_EXAMPLE
EXAMPLE_SUBDIR=example
endif

if BUILD_TESTS
tests: examples
	@echo ""
	@echo "Now building zoltan tests."
	@echo ""
	@echo "Zoltan tests use drivers and input filess; no compilation necssary."
	@echo ""
	@echo "Finished building zoltan tests."
	@echo ""
else
tests:
	@echo "zoltan tests were disabled at configure time"
endif

if BUILD_EXAMPLES
examples:
	@echo ""
	@echo "Now building zoltan examples."
	@echo ""
	cd $(top_builddir)/example && $(MAKE)
	@echo ""
	@echo "Finished building zoltan examples."
	@echo ""

install-examples:
	cd $(top_builddir)/example && $(MAKE) install
else
examples:
	@echo "zoltan examples were disabled at configure time"

install-examples:
	@echo "zoltan examples were disabled at configure time"
endif

clean-tests:
	cd $(top_builddir)/test && $(MAKE) clean

clean-examples:
	cd $(top_builddir)/example && $(MAKE) clean

everything:
	$(MAKE) && $(MAKE) examples && $(MAKE) tests

clean-everything:
	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

install-everything:
	$(MAKE) install && $(MAKE) install-examples

SUBDIRS = src $(EXAMPLE_SUBDIR) $(TEST_SUBDIR)

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
  --packages=zoltan

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
  --packages=zoltan

## #######################################################################
## Export makefile stuff
## #######################################################################

if USING_EXPORT_MAKEFILES

install-exec-hook:	
	cp $(top_builddir)/Makefile.export.zoltan $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.zoltan \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src \
		--my-abs-lib-dirs=@abs_top_builddir@/src
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile ZOLTAN > $(exec_prefix)/include/Makefile.export.zoltan.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.zoltan
	rm -f $(exec_prefix)/include/Makefile.export.zoltan.macros

else 

install-exec-hook:

uninstall-hook:

endif

## #######################################################################

