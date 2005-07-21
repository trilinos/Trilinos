# @HEADER
# ************************************************************************
# 
#       Ifpack: Object-Oriented Algebraic Preconditioner Package
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
ACLOCAL	= aclocal -I config

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
#EXTRA_DIST =


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
	@echo 
	@echo "        +----------------------------------------------------+"
	@echo "        | Trilinos package IFPACK built successfully         |"
	@echo "        |                                                    |"
	@echo "        | You might want to type \"make check\" to verify that |"
	@echo "        | all the IFPACK examples run successfully, or type  |"
	@echo "        | \"make check-all\" to run the test suite of IFPACK   |"
	@echo "        | (this may take some time).                         |"
	@echo "        |                                                    |"
	@echo "        | The IFPACK Users' Guide is located in              |"
	@echo "        | Trilinos/packages/ifpack/doc/UsersGuide, or        |"
	@echo "        | consult the IFPACK web site for more details.      |"
	@echo "        |                                                    |"
	@echo "        | You can send questions/comments to:                |"
	@echo "        | - ifpack-users@software.sandia.gov                 |"
	@echo "        | - ifpack-developers@software.sandia.gov            |"
	@echo "        |                                                    |"
	@echo "        | Have fun!                                          |"
	@echo "        +----------------------------------------------------+"
	@echo 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################
#SUBDIRS = src example test doc

if BUILD_TESTS
TEST_SUBDIR=test
endif

if BUILD_EXAMPLES
EXAMPLE_SUBDIR=example
endif 

if PYTHON_USE
PYTHON_SUBDIR=python
endif

SUBDIRS = src $(EXAMPLE_SUBDIR) $(TEST_SUBDIR) doc $(PYTHON_SUBDIR)

## #######################################################################
## Target for circular dependencies
## #######################################################################
all-libs:
	cd src ; make all-libs

## #######################################################################
## Export Makefile Installation
## #######################################################################
if USING_EXPORT_MAKEFILES

install-exec-hook:	
	cp $(top_builddir)/Makefile.export.ifpack $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.ifpack \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src \
		--my-abs-lib-dirs=@abs_top_builddir@/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../epetraext:@abs_top_builddir@/../aztecoo:@abs_top_builddir@/../teuchos:@abs_top_builddir@/../amesos:@abs_top_builddir@/../epetra:@abs_top_builddir@/../triutils
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile ifpack > $(exec_prefix)/include/Makefile.export.ifpack.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.ifpack
	rm -f $(exec_prefix)/include/Makefile.export.ifpack.macros

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

check-local:
	@cd example; make
	@echo
	@echo "============================="
	@echo "Testing the IFPACK package..."
	@echo "============================="
	@echo
	@echo "You are now running (most of) the executables in the"
	@echo "example subdirectory."
	@echo
if HAVE_MPI
	@echo "The MPI command is defaulted to \"mpirun -np \""
	@echo "You can change it by setting the environmental variable"
	@echo "TRILINOS_MPI_GO, for example in bash you might need"
	@echo "> export TRILINOS_MPI_GO=\"prun -n \""
	@echo
endif
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ML_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=IFPACKExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=ifpack

check-all:
	@cd test; make
	@echo
	@echo "============================="
	@echo "Testing the IFPACK package..."
	@echo "============================="
	@echo
	@echo "You are now running the complete IFPACK test suite (this may"
	@echo "test the executables in the examples subdirectory)."
	@echo
if HAVE_MPI
	@echo "By default, the MPI command is set to \"mpirun -np \"."
	@echo "You can change it by setting the environmental variable"
	@echo "TRILINOS_MPI_GO, for example in bash you might need"
	@echo "> export TRILINOS_MPI_GO=\"prun -n \""
	@echo
endif
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(ML_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=FRAMEWORK \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=ifpack
