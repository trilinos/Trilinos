# @HEADER
# ************************************************************************
# 
#               Isorropia: Partitioning and Load Balancing Package
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact Alan Williams (william@sandia.gov) 
#                 or Erik Boman    (egboman@sandia.gov)
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
EXTRA_DIST = config/generate-makeoptions.pl \
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
	@echo ""
	@echo "Trilinos package isorropia built successfully."
	@echo "For more details, please refer to the instructions"
	@echo "contained in Trilinos/doc/check_installation.txt"
	@echo ""
#	@echo "To check the installation, simply type \"$(MAKE) run-examples\"."

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

if BUILD_TESTS
TEST_SUBDIR=test
endif

if BUILD_EXAMPLES
EXAMPLE_SUBDIR=example
endif

# if PYTHON_USE
# PYTHON_SUBDIR=python
# endif

SUBDIRS = src $(EXAMPLE_SUBDIR) $(TEST_SUBDIR) $(PYTHON_SUBDIR)

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
  --packages=isorropia

runtests-mpi :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=mpi \
  --mpi-go=$(TRILINOS_MPI_GO) \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --category=INSTALL \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=isorropia

## #######################################################################
## Export makefile stuff
## #######################################################################

all-libs:
	cd $(top_builddir)/src ; $(MAKE) all-libs

libs:
	cd $(top_builddir)/src ; $(MAKE)

if USING_EXPORT_MAKEFILES

install-exec-hook:	
	cp $(top_builddir)/Makefile.export.isorropia $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.isorropia \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src \
		--my-abs-lib-dirs=@abs_top_builddir@/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../teuchos
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile ISORROPIA > $(exec_prefix)/include/Makefile.export.isorropia.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.isorropia
	rm -f $(exec_prefix)/include/Makefile.export.isorropia.macros

else 

install-exec-hook:

uninstall-hook:

endif

## #######################################################################

if HAVE_MPI
ISORROPIA_COMM=mpi
else
ISORROPIA_COMM=serial
endif

run-examples:
	@echo
	@echo "=============================="
	@echo "Running the Isorropia examples..."
	@echo "=============================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ISOROPIA_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=IsorropiaExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=isorropia

run-tests:
	@echo
	@echo "==========================="
	@echo "Running the Isorropia tests..."
	@echo "==========================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ISORROPIA_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=FRAMEWORK \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=isorropia

run-pyexamples:
	@echo
	@echo "====================================="
	@echo "Running the Python/Isorropia examples..."
	@echo "====================================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ISORROPIA_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=IsorropiaPyExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=isorropia

run-pytests:
	@echo
	@echo "=================================="
	@echo "Running the Python/Isorropia tests..."
	@echo "=================================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ISORROPIA_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=IsorropiaPyTests \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=isorropia
