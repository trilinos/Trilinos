# @HEADER
# ************************************************************************
# 
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2005) Sandia Corporation
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
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
EXTRA_DIST = doc/Doxygen/Teuchos_dox.i                  \
	     doc/Doxygen/Epetra_dox.i                   \
	     doc/Doxygen/TriUtils_dox.i                 \
	     doc/Doxygen/EpetraExt_dox.i                \
	     doc/Doxygen/AztecOO_dox.i                  \
	     doc/Doxygen/Galeri_dox.i                   \
	     doc/Doxygen/Amesos_dox.i                   \
	     doc/Doxygen/IFPACK_dox.i                   \
	     doc/Doxygen/Anasazi_dox.i                  \
	     doc/Doxygen/ML_dox.i                       \
	     doc/Doxygen/NOX_dox.i                      \
	     doc/build_docs                             \
	     doc/DoxyfileWeb                            \
	     doc/DevelopersGuide/AddingModules.txt      \
	     doc/DevelopersGuide/AddingModules.html     \
	     doc/DevelopersGuide/BuildSystem.txt        \
	     doc/DevelopersGuide/BuildSystem.html       \
	     doc/DevelopersGuide/Documentation.txt      \
	     doc/DevelopersGuide/Documentation.html     \
	     doc/DevelopersGuide/ExceptionHandling.txt  \
	     doc/DevelopersGuide/ExceptionHandling.html \
	     doc/DevelopersGuide/Testing.txt            \
	     doc/DevelopersGuide/Testing.html           \
	     INSTALL
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
	$(AUX_DIST) 

#The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package PyTrilinos built successfully"
	@echo ""

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though BUILD_TESTS and BUILD_EXAMPLES will never be
#defined, so that the tests and examples are included in the distribution
#tarball.

if BUILD_PYTRILINOS
SRC_SUBDIR=src
endif

if SUB_TEST
TEST_SUBDIR=test
endif

if SUB_EXAMPLE
EXAMPLE_SUBDIR=example
endif

if BUILD_PYTRILINOS
if BUILD_TESTS
tests:
	@echo ""
	@echo "Now building PyTrilinos tests."
	@echo ""
	cd $(top_builddir)/test && $(MAKE)
	@echo ""
	@echo "Finished building PyTrilinos tests."
	@echo ""
else
tests:
	@echo "PyTrilinos tests were disabled at configure time"
endif

if BUILD_EXAMPLES
examples:
	@echo ""
	@echo "Now building PyTrilinos examples."
	@echo ""
	cd $(top_builddir)/example && $(MAKE)
	@echo ""
	@echo "Finished building PyTrilinos examples."
	@echo ""

install-examples:
	cd $(top_builddir)/example && $(MAKE) install
else
examples:
	@echo "PyTrilinos examples were disabled at configure time"

install-examples:
	@echo "PyTrilinos examples were disabled at configure time"
endif
else
# BUILD_PYTRILINOS = no
tests:
	@echo "PyTrilinos was disabled at configure time"

examples:
	@echo "PyTrilinos was disabled at configure time"

install-examples:
	@echo "PyTrilinos was disabled at configure time"
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

SUBDIRS = shared $(SRC_SUBDIR) $(TEST_SUBDIR) $(EXAMPLE_SUBDIR)

## #######################################################################
## The below targets allow you to use the new
## testharness to run the test suite as make targets
## #######################################################################

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..
TRILINOS_TEST_CATEGORY=INSTALL

runtests-serial :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=serial \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=$(TRILINOS_TEST_CATEGORY) \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=PyTrilinos

runtests-mpi :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=mpi \
	--mpi-go=$(TRILINOS_MPI_GO) \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=$(TRILINOS_TEST_CATEGORY) \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=PyTrilinos
