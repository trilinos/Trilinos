## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
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
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov), Sandia National Laboratories.
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
EXTRA_DIST = LICENSE config/generate-makeoptions.pl \
       config/strip_dup_incl_paths.pl config/strip_dup_libs.pl \
       config/replace-install-prefix.pl config/string-replace.pl


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
	@echo "Trilinos package nox built successfully"

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

SUBDIRS = doc src src-lapack src-epetra src-petsc src-belos \
	  src-loca test examples python

## #######################################################################
## Documentation
## #######################################################################

dox:
	cd doc;	$(MAKE) dox; 

## #######################################################################
## Target for circular dependencies
## #######################################################################
all-libs:
	cd src ; $(MAKE) all-libs
	cd ../src-lapack; $(MAKE)
	cd ../src-epetra; $(MAKE)
	cd ../src-petsc; $(MAKE)
	cd ../src-belos; $(MAKE)
	cd ../src-loca/src; $(MAKE)
	cd ../src-loca/src-lapack; $(MAKE)
	cd ../src-loca/src-epetra; $(MAKE)
	cd ../src-loca/src-mf; $(MAKE)

libs:
	cd src ; $(MAKE)
	cd ../src-lapack; $(MAKE)
	cd ../src-epetra; $(MAKE)
	cd ../src-petsc; $(MAKE)
	cd ../src-belos; $(MAKE)
	cd ../src-loca/src; $(MAKE)
	cd ../src-loca/src-lapack; $(MAKE)
	cd ../src-loca/src-epetra; $(MAKE)
	cd ../src-loca/src-mf; $(MAKE)

## #######################################################################
## Export Makefile Installation
## #######################################################################
if USING_EXPORT_MAKEFILES

install-exec-hook:	
	cp $(top_builddir)/Makefile.export.nox $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.nox \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src:@abs_top_srcdir@/src-epetra:@abs_top_srcdir@/src-lapack:@abs_top_srcdir@/src-petsc:@abs_top_srcdir@/src-belos:@abs_top_srcdir@/src-loca/src:@abs_top_srcdir@/src-loca/src-lapack:@abs_top_srcdir@/src-loca/src-epetra:@abs_top_srcdir@/src-loca/src-mf \
		--my-abs-lib-dirs=@abs_top_builddir@/src:@abs_top_builddir@/src-epetra:@abs_top_builddir@/src-lapack:@abs_top_builddir@/src-petsc:@abs_top_builddir@/src-belos:@abs_top_builddir@/src-loca/src:@abs_top_builddir@/src-loca/src-epetra:@abs_top_builddir@/src-loca/src-lapack:@abs_top_builddir@/src-loca/src-mf \
		--dep-package-abs-builddirs=@abs_top_builddir@/../epetraext:@abs_top_builddir@/../aztecoo:@abs_top_builddir@/../teuchos:@abs_top_builddir@/../ifpack:@abs_top_builddir@/../epetra:@abs_top_builddir@/../ml:@abs_top_builddir@/../belos:@abs_top_builddir@/../anasazi:@abs_top_builddir@/../TSFExtended
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile nox > $(exec_prefix)/include/Makefile.export.nox.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.nox
	rm -f $(exec_prefix)/include/Makefile.export.nox.macros

else 

install-exec-hook:

uninstall-hook:

endif

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..

if HAVE_MPI
NOX_CHECK_COMM=mpi
else
NOX_CHECK_COMM=serial
endif

run-examples:
	@echo
	@echo "==========================="
	@echo "Running the NOX examples..."
	@echo "==========================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(NOX_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=NOXExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=nox

run-tests:
	@echo
	@echo "========================"
	@echo "Running the NOX tests..."
	@echo "========================"
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(NOX_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=FRAMEWORK \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=nox

run-pyexamples:
	@echo
	@echo "=================================="
	@echo "Running the Python/NOX examples..."
	@echo "=================================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(NOX_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=NOXPyExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=nox

run-pytests:
	@echo
	@echo "==============================="
	@echo "Running the Python/NOX tests..."
	@echo "==============================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(NOX_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=NOXPyTests \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=nox
