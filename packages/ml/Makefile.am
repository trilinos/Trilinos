# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
# Typically the only line that needs to be edited in this file is the last 
# line.  Set 'SUBDIRS= (list all subdirectories that are autotool'ed here)'
# List all subdirectories even if one or more are configured/built only
# conditionally.

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
EXTRA_DIST = ChangeLog-ML README-ML COPYRIGHT \
       config/generate-makeoptions.pl \
       config/strip_dup_incl_paths.pl config/strip_dup_libs.pl \
       config/replace-install-prefix.pl config/string-replace.pl \
       doc/Doxyfile doc/mainpage.doxygen doc/mlguide.pdf \
       doc/footer.html doc/header.tex doc/ml-image.jpg \
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
	@if test ! -e $(top_builddir)/examples/XML/ml_ParameterList.xml; then cp $(top_srcdir)/examples/XML/ml_ParameterList.xml $(top_builddir)/examples/XML; fi
	@echo "Trilinos package ml built successfully"

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
EXAMPLE_SUBDIR=examples
endif

if BUILD_ML_MATLAB
MATLAB_SUBDIR=matlab
endif

if BUILD_TESTS
tests: examples
	@echo ""
	@echo "Now building ml tests."
	@echo ""
	cd $(top_builddir)/test && $(MAKE)
	@echo ""
	@echo "Finished building ml tests."
	@echo ""
else
tests:
	@echo "ml tests were disabled at configure time"
endif

if ML_BUILD_EXAMPLES
examples: FORCE
	@echo ""
	@echo "Now building ml examples."
	@echo ""
	cd $(top_builddir)/examples && $(MAKE)
	@echo ""
	@echo "Finished building ml examples."
	@echo ""

install-examples:
	cd $(top_builddir)/examples && $(MAKE) install
else
examples: FORCE
	@echo "ml examples were disabled at configure time"

install-examples:
	@echo "ml examples were disabled at configure time"
endif

FORCE:

clean-tests:
	cd $(top_builddir)/test && $(MAKE) clean

clean-examples:
	cd $(top_builddir)/examples && $(MAKE) clean

everything:
	$(MAKE) && $(MAKE) examples && $(MAKE) tests

clean-everything:
	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

install-everything:
	$(MAKE) install && $(MAKE) install-examples

SUBDIRS = src $(EXAMPLE_SUBDIR) $(TEST_SUBDIR) $(MATLAB_SUBDIR)

## #######################################################################
## Export Makefile Installation
## #######################################################################
if USING_EXPORT_MAKEFILES

install-exec-hook:
	mkdir -p $(DESTDIR)$(includedir)
	cp $(top_builddir)/Makefile.export.ml $(DESTDIR)$(includedir)/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.ml \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src:@abs_top_srcdir@/Comm:abs_top_srcdir@/Operator:abs_top_srcdir@/Smoother:@abs_top_srcdir@/Coarsen:@abs_top_srcdir@/Krylov:@abs_top_srcdir@/Main:@abs_top_srcdir@/FEGrid:@abs_top_srcdir@/Utils:@abs_top_srcdir@/MLAPI:@abs_top_srcdir@/Include \
		--my-abs-lib-dirs=@abs_top_builddir@/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../galeri:@abs_top_builddir@/../epetraext:@abs_top_builddir@/../epetra:@abs_top_builddir@/../triutils:@abs_top_builddir@/../amesos:@abs_top_builddir@/../ifpack:@abs_top_builddir@/../teuchos:@abs_top_builddir@/../anasazi:@abs_top_builddir@/../aztecoo
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile ML > $(DESTDIR)$(includedir)/Makefile.export.ml.macros

uninstall-hook:
	rm -f $(includedir)/Makefile.export.ml
	rm -f $(includedir)/Makefile.export.ml.macros

else 

install-exec-hook:

uninstall-hook:

endif

## #######################################################################
## The below targets allow you to use the new
## testharness to run the test suite as make targets
## #######################################################################

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..

if HAVE_MPI
ML_CHECK_COMM=mpi
else
ML_CHECK_COMM=serial
endif

TRILINOS_MPI_MAX_PROC=4
TRILINOS_TEST_CATEGORY=FRAMEWORK

runtests-serial :
	@echo "Running the ML tests."
	@echo "  TRILINOS_TEST_CATEGORY = $(TRILINOS_TEST_CATEGORY)"
	@echo ""
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=serial \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --category=$(TRILINOS_TEST_CATEGORY) \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=ml

runtests-mpi :
	@echo "Running the ML tests."
	@echo "  TRILINOS_TEST_CATEGORY = $(TRILINOS_TEST_CATEGORY)"
	@echo "  TRILINOS_MPI_GO = $(TRILINOS_MPI_GO)"
	@echo ""
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=mpi \
  --mpi-go=$(TRILINOS_MPI_GO) \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --max-proc=$(TRILINOS_MPI_MAX_PROC) \
  --category=$(TRILINOS_TEST_CATEGORY) \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=ml

