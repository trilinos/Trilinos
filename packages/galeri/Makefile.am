#np# Typically the only line that needs to be edited in this file is the last 
#np# line.  Set 'SUBDIRS= (list all subdirectories that are autotool'ed here)'
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.

# @HEADER
# ************************************************************************
#
#                  Galeri Matrix Generation Package
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

#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
#np# For a typical package, there is no reason to distribute these files
#np# because users should not have to bootstrap.  We distribute them with
#np# new package so that the files can be used in creating the 
#np# configure script for other packages. 
EXTRA_DIST = \
config/ac_cxx_bool.m4                   config/depcomp \
config/ac_cxx_have_std.m4               config/install-sh \
config/ac_cxx_have_stl.m4               config/missing \
config/ac_cxx_mutable.m4                config/mkinstalldirs \
config/ac_cxx_namespaces.m4             config/tac_arg_check_mpi.m4 \
config/ac_cxx_new_for_scoping.m4        config/tac_arg_config_mpi.m4 \
config/ac_cxx_reinterpret_cast.m4 \
config/tac_arg_enable_compiler_options.m4 \
config/ac_cxx_std_ios_base_fmtflags.m4  config/tac_arg_enable_feature.m4 \
config/ac_cxx_std_sprintf.m4            config/tac_arg_enable_option.m4 \
config/ac_prompt_user.m4                config/tac_arg_with_ar.m4 \
config/ac_prompt_user_no_define.m4      config/tac_arg_with_blaslib.m4 \
config/acx_blas.m4                      config/tac_arg_with_compiler.m4 \
config/acx_check_dos_filesys.m4         config/tac_arg_with_flags.m4 \
config/acx_lapack.m4                    config/tac_arg_with_incdirs.m4 \
config/acx_pthread.m4                   config/tac_arg_with_lapacklib.m4 \
config/config.guess                     config/tac_arg_with_libdirs.m4 \
config/config.sub                       config/tac_arg_with_package.m4 \
config/strip_dup_incl_paths.pl		config/strip_dup_libs.pl \
config/generate-makeoptions.pl		config/replace-install-prefix.pl \
config/string-replace.pl		bootstrap \
README-Galeri \
ChangeLog-Galeri \
TODO-Galeri \
data/Square.grid \
data/TwoSquares.grid \
data/VerySimple.grid \
data/Circle.grid \
data/SquareInSquare.grid \
data/TwoSquares.m \
data/convert.m

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
	@echo 
	@echo "        +----------------------------------------------------+"
	@echo "        | Trilinos package Galeri built successfully         |"
	@echo "        |                                                    |"
	@echo "        | You might want to type \"$(MAKE) run-examples\" to      |"
	@echo "        | verify that all the Galeri examples run            |"
	@echo "        | successfully. For more details, please refer to    |"
	@echo "        | Trilinos/doc/check_installation.txt.               |"
	@echo "        |                                                    |"
	@echo "        | Please send questions/comments to:                 |"
	@echo "        | - galeri-users@software.sandia.gov                 |"
	@echo "        | - galeri-developers@software.sandia.gov            |"
	@echo "        |                                                    |"
	@echo "        | Have fun!                                          |"
	@echo "        +----------------------------------------------------+"
	@echo 

all-libs:
	cd $(top_builddir)/src ; $(MAKE) all-libs

libs:
	cd $(top_builddir)/src ; $(MAKE)

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

if BUILD_TESTS
TEST_SUBDIR=test
endif

if BUILD_EXAMPLES
EXAMPLE_SUBDIR=example example-fem
endif

if BUILD_PYTHON
if BUILD_GALERI_PYTHON
PYTHON_SUBDIR=python
endif
endif

SUBDIRS = src $(TEST_SUBDIR) $(EXAMPLE_SUBDIR) $(PYTHON_SUBDIR)

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..

## #######################################################################
## Export makefile stuff
## #######################################################################

all-libs:
	cd $(top_builddir)/src ; $(MAKE) all-libs

libs:
	cd $(top_builddir)/src ; $(MAKE)

if USING_EXPORT_MAKEFILES

install-exec-hook:
	cp $(top_builddir)/Makefile.export.galeri $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
                --exec-prefix=$(exec_prefix) \
                --my-export-makefile=Makefile.export.galeri \
                --my-abs-top-srcdir=@abs_top_srcdir@ \
                --my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src:@abs_top_srcdir@/Maps:@abs_top_srcdir@/FiniteElements:@abs_top_srcdir@/CrsMatrices \
                --my-abs-lib-dirs=@abs_top_builddir@/src \
                --dep-package-abs-builddirs=@abs_top_builddir@/../teuchos:@abs_top_builddir@/../epetra
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile GALERI > $(exec_prefix)/include/Makefile.export.galeri.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.galeri
	rm -f $(exec_prefix)/include/Makefile.export.galeri.macros

else

install-exec-hook:

uninstall-hook:

endif

## #######################################################################

if HAVE_MPI
GALERI_CHECK_COMM=mpi
else
GALERI_CHECK_COMM=serial
endif

run-examples:
	@echo
	@echo "=============================="
	@echo "Running the Galeri examples..."
	@echo "=============================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(GALERI_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=GaleriExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=galeri

run-tests:
	@echo
	@echo "==========================="
	@echo "Running the Galeri tests..."
	@echo "==========================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(GALERI_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=FRAMEWORK \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=galeri

run-pyexamples:
	@echo
	@echo "====================================="
	@echo "Running the Python/Galeri examples..."
	@echo "====================================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(GALERI_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=GaleriPyExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=galeri

run-pytests:
	@echo
	@echo "=================================="
	@echo "Running the Python/Galeri tests..."
	@echo "=================================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(GALERI_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=GaleriPyTests \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=galeri
