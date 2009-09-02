# @HEADER
# ************************************************************************
# 
#                Komplex: Complex Linear Solver Package
#                 Copyright (2002) Sandia Corporation
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

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = config/strip_dup_incl_paths.pl config/strip_dup_libs.pl \
             config/token-replace.pl config/replace-install-prefix.pl \
             config/generate-makeoptions.pl \
             doc/build_docs doc/DoxyfileWeb \
             doc/index.doc \
             example/simple/Makefile \
             example/simple/README \
             example/simple/cxx_main.cpp \
             example/simple/invoke-configure

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
	@echo "Trilinos package komplex built successfully"

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though BUILD_TESTS and BUILD_EXAMPLES will never be
#defined, so that the tests and examples are included in the distribution
#tarball.  (Komplex does not have any tests at this time.)

if SUB_EXAMPLE
EXAMPLE_SUBDIR=example
endif

tests:
	@echo "komplex does not have any tests yet"

if BUILD_EXAMPLES
examples:
	@echo ""
	@echo "Now building komplex examples."
	@echo ""
	cd $(top_builddir)/example && $(MAKE)
	@echo ""
	@echo "Finished building komplex examples."
	@echo ""

install-examples:
	cd $(top_builddir)/example && $(MAKE) install
else
examples:
	@echo "komplex examples were disabled at configure time"

install-examples:
	@echo "komplex examples were disabled at configure time"
endif

clean-tests:
	@echo "komplex does not have any tests yet"

clean-examples:
	cd $(top_builddir)/example && $(MAKE) clean

everything:
	$(MAKE) && $(MAKE) examples && $(MAKE) tests

clean-everything:
	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

install-everything:
	$(MAKE) install && $(MAKE) install-examples

SUBDIRS = src $(EXAMPLE_SUBDIR)

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
  --packages=komplex

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
  --packages=komplex

## #######################################################################
## Export makefile stuff
## #######################################################################

if USING_EXPORT_MAKEFILES

install-exec-hook:	
	mkdir -p $(DESTDIR)$(includedir)
	cp $(top_builddir)/Makefile.export.komplex $(DESTDIR)$(includedir)/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.komplex \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src \
		--my-abs-lib-dirs=@abs_top_builddir@/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../teuchos
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile KOMPLEX > $(DESTDIR)$(includedir)/Makefile.export.komplex.macros

uninstall-hook:
	rm -f $(includedir)/Makefile.export.komplex
	rm -f $(includedir)/Makefile.export.komplex.macros

else 

install-exec-hook:

uninstall-hook:

endif

## #######################################################################
