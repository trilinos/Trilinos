# @HEADER
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
## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
#np# For a typical package, there is no reason to distribute these files
#np# because users should not have to bootstrap.  We distribute them with
#np# new package so that the files can be used in creating the 
#np# configure script for other packages. 
EXTRA_DIST = \
config/config.guess  config/config.sub  config/install-sh  \
config/token-replace.pl config/missing bootstrap 

## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
#AUX_DIST = config/install-sh config/missing config/mkinstalldirs 
#

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	src/common/config-h.in src/common/stamp-h.in \
	$(AUX_DIST) 

#The following line helps the test harness recover from build errors.
                                                                                
all-local:
	@echo "Trilinos package external built successfully"

## #######################################################################
## Export Makefile Installation
## #######################################################################
#if USING_EXPORT_MAKEFILES

#Fix this as necessary!
#install-exec-hook:
#	cp $(top_builddir)/Makefile.export.new_package $(exec_prefix)/include/.
#	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
#		--exec-prefix=$(exec_prefix) \
#		--my-export-makefile=Makefile.export.new_package \
#		--my-abs-top-srcdir=@abs_top_srcdir@ \
#		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src \
#		--my-abs-lib-dirs=@abs_top_builddir@/src \
#		--dep-package-abs-builddirs=@abs_top_builddir@/../epetra
#	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile \
#		new_package > $(exec_prefix)/include/Makefile.export.new_package.macros

#uninstall-hook:
#	rm -f $(exec_prefix)/include/Makefile.export.new_package
#	rm -f $(exec_prefix)/include/Makefile.export.new_package.macros

#else

#install-exec-hook:

#uninstall-hook:

#endif

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though BUILD_TESTS and BUILD_EXAMPLES will never be
#defined, so that the tests and examples are included in the distribution
#tarball.

#if SUB_TEST
#TEST_SUBDIR=test
#endif

#if SUB_EXAMPLE
#EXAMPLE_SUBDIR=example
#endif

#if BUILD_TESTS
#tests:
#	cd $(top_builddir)/test && $(MAKE)
#else
#tests:
#	@echo "sample tests were disabled at configure time"
#endif

#if BUILD_EXAMPLES
#examples:
#	cd $(top_builddir)/example && $(MAKE)

#install-examples:
#	cd $(top_builddir)/example && $(MAKE) install
#else
#examples:
#	@echo "sample examples were disabled at configure time"

#install-examples:
#	@echo "sample examples were disabled at configure time"
#endif

#clean-tests:
#	cd $(top_builddir)/test && $(MAKE) clean

#clean-examples:
#	cd $(top_builddir)/example && $(MAKE) clean

#everything:
#	$(MAKE) && $(MAKE) examples && $(MAKE) tests

#clean-everything:
#	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

#install-everything:
#	$(MAKE) install && $(MAKE) install-examples

#TEST_SUBDIR=EXTERNAL_PACKAGE_NAME

#SUBDIRS = $(TEST_SUBDIR)

#TRILINOS_HOME_DIR=@abs_top_srcdir@/../../../../
#TRILINOS_BUILD_DIR=@abs_top_builddir@/../../../../

