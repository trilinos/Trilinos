# Typically the only line that needs to be edited in this file is the last 
# line.  Set 'SUBDIRS= (list all subdirectories that are autotool'ed here)'
# List all subdirectories even if one or more are configured/built only
# conditionally.

# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001,2002,2003) Sandia Corporation
#

# Copyright (2001,2002,2003) Sandia Corportation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in ths data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from

# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
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
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = ChangeLog-ML README-ML COPYRIGHT \
       config/generate-makeoptions.pl \
       config/strip_dup_incl_paths.pl config/strip_dup_libs.pl \
       config/replace-install-prefix.pl config/string-replace.pl


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
	@cp $(top_srcdir)/examples/BasicExamples/ml_ParameterList.xml $(top_builddir)/examples/BasicExamples
	@echo "Trilinos package ml built successfully"
	@echo 
	@echo "        +----------------------------------------------------+"
	@echo "        | Trilinos package ML built successfully             |"
	@echo "        |                                                    |"
	@echo "        | You might want to type \"$(MAKE) run-examples\" to      |"
	@echo "        | verify that all the ML examples run                |"
	@echo "        | successfully. For more details, please refer to    |"
	@echo "        | Trilinos/doc/check_installation.txt.               |"
	@echo "        |                                                    |"
	@echo "        | The ML Users' Guide is located in                  |"
	@echo "        | Trilinos/packages/ml/doc/mlguide.pdf. You can also |"
	@echo "        | consult the ML web site for more details.          |"
	@echo "        |                                                    |"
	@echo "        | Please send questions/comments to:                 |"
	@echo "        | - ml-users@software.sandia.gov                     |"
	@echo "        | - ml-developers@software.sandia.gov                |"
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

if ML_BUILD_EXAMPLES
EXAMPLE_SUBDIR=examples
endif

if PYTHON_USE
if BUILD_ML_PYTHON
PYTHON_SUBDIR=python
endif
endif

if BUILD_ML_MATLAB
MATLAB_SUBDIR=matlab
endif


if ENABLE_ML_THYRA
THYRA_SUBDIR=thyra
else
THYRA_SUBDIR=
endif

SUBDIRS = src $(THYRA_SUBDIR) $(EXAMPLE_SUBDIR) $(TEST_SUBDIR) doc $(PYTHON_SUBDIR) $(MATLAB_SUBDIR)


## #######################################################################
## Target for circular dependencies
## #######################################################################

if ENABLE_ML_THYRA
BUILD_THYRA_SUBLIB=cd $(top_builddir)/thyra/src ; $(MAKE)
else
BUILD_THYRA_SUBLIB=
endif

all-libs:
	cd src ; $(MAKE) all-libs
	$(BUILD_THYRA_SUBLIB)

libs:
	cd src ; $(MAKE)
	$(BUILD_THYRA_SUBLIB)

## #######################################################################
## Export Makefile Installation
## #######################################################################
if USING_EXPORT_MAKEFILES

install-exec-hook:
	cp $(top_builddir)/Makefile.export.ml $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.ml \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src:@abs_top_srcdir@/Comm:abs_top_srcdir@/Operator:abs_top_srcdir@/Smoother:@abs_top_srcdir@/Coarsen:@abs_top_srcdir@/Krylov:@abs_top_srcdir@/Main:@abs_top_srcdir@/FEGrid:@abs_top_srcdir@/Utils:@abs_top_srcdir@/MLAPI:@abs_top_srcdir@/Include \
		--my-abs-lib-dirs=@abs_top_builddir@/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../epetraext:@abs_top_builddir@/../epetra:@abs_top_builddir@/../triutils:@abs_top_builddir@/../amesos:@abs_top_builddir@/../ifpack:@abs_top_builddir@/../teuchos:@abs_top_builddir@/../anasazi:@abs_top_builddir@/../aztecoo
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile ML > $(exec_prefix)/include/Makefile.export.ml.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.ml
	rm -f $(exec_prefix)/include/Makefile.export.ml.macros

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

TRILINOS_TEST_CATEGORY=FRAMEWORK

run-examples:
	@echo
	@echo "=========================="
	@echo "Running the ML examples..."
	@echo "=========================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ML_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=MLExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=ml

run-tests:
	@echo
	@echo "====================================================="
	@echo "Running the ML tests in category ${TRILINOS_TEST_CATEGORY}"
	@echo "====================================================="
	@echo
	@echo "To change categories, type"
	@echo "  make run-tests TRILINOS_TEST_CATEGORY=some_test_category"
	@echo "See ml/test/definition for your choices."
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ML_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=FRAMEWORK \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=ml

run-pyexamples:
	@echo
	@echo "================================="
	@echo "Running the Python/ML examples..."
	@echo "================================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ML_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=MLPyExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=ml

run-pytests:
	@echo
	@echo "=============================="
	@echo "Running the Python/ML tests..."
	@echo "=============================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ML_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=MLPyTests \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=ml
