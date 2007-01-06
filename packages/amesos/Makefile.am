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
EXTRA_DIST = README-Amesos ChangeLog-Amesos \
       config/generate-makeoptions.pl \
       config/strip_dup_incl_paths.pl config/strip_dup_libs.pl \
       config/replace-install-prefix.pl config/string-replace.pl \
       doc/Doxyfile doc/AmesosLogo.eps doc/AmesosLogo.jpeg \
       doc/Doxyfile doc/footer.html doc/index.doc \
       doc/header.tex doc/AmesosReferenceGuide/AmesosReferenceGuide.pdf


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

#The following line helps the test harness recover from build errors.

all-local:
	@cp -f $(srcdir)/example/In_662_bus.rsa $(top_builddir)/example/662_bus.rsa
	@cp -f $(srcdir)/example/_stratimikos_amesos_lapack.xml $(top_builddir)/example/stratimikos_amesos_lapack.xml
	@cp -f $(srcdir)/example/_FourByFour.mtx $(top_builddir)/example/FourByFour.mtx
	@echo 
	@echo "        +----------------------------------------------------+"
	@echo "        | Trilinos package amesos built successfully         |"
	@echo "        |                                                    |"
	@echo "        | You might want to type \"make run-examples\" to      |"
	@echo "        | verify that all the Amesos examples run            |"
	@echo "        | successfully. For more details, please refer to    |"
	@echo "        | Trilinos/doc/check_installation.txt.               |"
	@echo "        |                                                    |"
	@echo "        | The Amesos Users' Guide is located in              |"
	@echo "        | Trilinos/packages/amesos/doc/AmesosReferenceGuide. |"
	@echo "        | Also you can consult the Amesos web page.          |"
	@echo "        |                                                    |"
	@echo "        | Please send questions/comments to:                 |"
	@echo "        | - amesos-users@software.sandia.gov                 |"
	@echo "        | - amesos-developers@software.sandia.gov            |"
	@echo "        |                                                    |"
	@echo "        | Have fun!                                          |"
	@echo "        +----------------------------------------------------+"
	@echo 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

if BUILD_TESTS
TEST_SUBDIR=test
endif

if BUILD_EXAMPLES
EXAMPLE_SUBDIR=example
endif

if ENABLE_AMESOS_THYRA
THYRA_SUBDIR=thyra
else
THYRA_SUBDIR=
endif

SUBDIRS = src $(THYRA_SUBDIR) $(EXAMPLE_SUBDIR) $(TEST_SUBDIR)

## #######################################################################
## Target for circular dependencies
## #######################################################################

if ENABLE_AMESOS_THYRA
BUILD_THYRA_SUBLIB=cd $(top_builddir)/thyra/src ; $(MAKE)
else
BUILD_THYRA_SUBLIB=
endif

all-libs:
	cd $(top_builddir)/src ; $(MAKE) all-libs
	$(BUILD_THYRA_SUBLIB)

libs:
	cd $(top_builddir)/src ; $(MAKE)
	$(BUILD_THYRA_SUBLIB)

## #######################################################################
## Export Makefile Installation
## #######################################################################
if USING_EXPORT_MAKEFILES

install-exec-hook:	
	cp $(top_builddir)/Makefile.export.amesos $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.amesos \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src:@abs_top_srcdir@/src/klu_btf:@abs_top_srcdir@/src/klu_kernel:@abs_top_srcdir@/src/btf:@abs_top_srcdir@/src/AMD:@abs_top_srcdir@/src/AMD/Include:@abs_top_srcdir@/src/AMD/Source:@abs_top_srcdir@/src/colamd:@abs_top_srcdir@/src/ksparse \
		--my-abs-lib-dirs=@abs_top_builddir@/src:@abs_top_builddir@/thyra/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../epetraext:@abs_top_builddir@/../epetra:@abs_top_builddir@/../thyra:@abs_top_builddir@/../teuchos:@abs_top_builddir@/../triutils:@abs_top_builddir@/../galeri
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile AMESOS > $(exec_prefix)/include/Makefile.export.amesos.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.amesos
	rm -f $(exec_prefix)/include/Makefile.export.amesos.macros

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

runtests-serial :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=serial \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --category=INSTALL \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=amesos

if HAVE_MPI
AMESOS_CHECK_COMM=mpi
else
AMESOS_CHECK_COMM=serial
endif

run-examples:
	@echo
	@echo "=============================="
	@echo "Running the Amesos examples..."
	@echo "=============================="
	@echo
	@echo "None of these examples should take more than a minute to execute on a typical compter"
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(AMESOS_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=AmesosExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=amesos

run-tests:
	@echo
	@echo "==========================="
	@echo "Running the Amesos tests..."
	@echo "==========================="
	@echo
	@echo "None of these tests should take more than a minute or two to execute on a multi-processor computer"
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo "make run-full-tests perform a more complete set of tests"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(AMESOS_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=INSTALL \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=amesos

run-full-tests:
	@echo
	@echo "==========================="
	@echo "Running ALL Amesos tests..."
	@echo "==========================="
	@echo
	@echo "None of these tests should take more than an hour to execute on a typical computer"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(AMESOS_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=FRAMEWORK \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=amesos
