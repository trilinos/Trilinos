# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#

# Copyright (2001) Sandia Corportation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in ths data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
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
ACLOCAL_AMFLAGS = -I config

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = \
	Trilinos_version.h \
	commonTools/buildTools/strip_dup_incl_paths.pl \
	commonTools/buildTools/strip_dup_libs.pl \
	commonTools/test/coverage/README \
	commonTools/test/coverage/gcovmerge \
	commonTools/test/harness/config \
	commonTools/test/harness/README \
	commonTools/test/harness/runharness \
	commonTools/test/harness/invoke-configures/linux-beowulf \
	commonTools/test/harness/invoke-configures/linux-beowulf-mpi \
	commonTools/test/utilities/compareTestOutput \
	commonTools/test/utilities/dependencies \
	commonTools/test/utilities/packages \
	commonTools/test/utilities/README \
	commonTools/test/utilities/report \
	commonTools/test/utilities/runbuild \
	commonTools/test/utilities/runconfigure \
	commonTools/test/utilities/runmake \
	commonTools/test/utilities/runtests \
	sampleScripts/aix_mpi_xlC sampleScripts/aix_mpi_xlC2 \
	sampleScripts/compaqCluster_mpi_stratus \
	sampleScripts/compaqCluster_serial_stratus sampleScripts/cplant_mpi_ross \
	sampleScripts/cplant_serial_ross \
	sampleScripts/cygwin_mpi sampleScripts/cygwin_serial_gcc \
	sampleScripts/cygwin_serial_icl sampleScripts/cygwingcc_amesos \
	sampleScripts/irix_mpi_sierra sampleScripts/janus_mpi \
	sampleScripts/janus_mpi_stlport \
	sampleScripts/linux_mpi_herouxsmp sampleScripts/linux_mpi_instCluster \
	sampleScripts/macosx_mpi_paul sampleScripts/macosx_serial_paul \
	sampleScripts/macosx_serial_samt5980 \
	sampleScripts/macosx_mpi_marzio \
	sampleScripts/macosx_serial_marzio \
	sampleScripts/linux_mpi_marzio \
	sampleScripts/kythira_serial_marzio \
	sampleScripts/dec_mpi_marzio \
	sampleScripts/cray_xt3_marzio \
	sampleScripts/scyld_mpi_beowulf sampleScripts/sgi64_mpi_atlantis \
	sampleScripts/sgi64_serial_atlantis \
	sampleScripts/solaris_mpi_conejo sampleScripts/solaris_serial_paul \
	sampleScripts/README config/AutotoolsVersionInfo \
	config/makeall \
	sampleScripts/cygwin_icl.pl sampleScripts/stratus_ar.pl \
	doc/Tutorial/README RELEASE_NOTES \
	doc/build_docs.pl \
	commonTools/refactoring/copy-files-list.pl \
	commonTools/refactoring/source_tree_stats.pl \
	commonTools/refactoring/c++-to-c-comments \
	commonTools/refactoring/string-replace-list.pl \
	commonTools/refactoring/string-replace-list-r \
	commonTools/refactoring/fix-dos-text-file.pl \
	commonTools/refactoring/string-replace.pl \
	commonTools/refactoring/fix-dos-text-file-r \
	commonTools/refactoring/string-replace-r \
	commonTools/refactoring/fix-doxygen-documentation.pl \
	commonTools/refactoring/token-replace-list.pl \
	commonTools/refactoring/fix-doxygen-documentation-r \
	commonTools/refactoring/token-replace-list-r \
	commonTools/refactoring/README \
	commonTools/refactoring/token-replace.pl \
	commonTools/refactoring/replace-tabs-with-spaces.pl \
	commonTools/refactoring/token-replace-r \
	commonTools/refactoring/replace-tabs-with-spaces-r \
	commonTools/refactoring/update-copyright-header.pl \
	commonTools/refactoring/source_file_stats.py \
	CMakeLists.txt \
	cmake/Trilinos_Add_Executable.cmake \
	cmake/TrilinosConfig.cmake.in \
	cmake/ConfigureMPI.cmake

## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs config/makeall

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
SUBDIRS = packages

## #######################################################################
## Install the perl scrpts used for external makefiles
## #######################################################################

install-exec-hook:	
	cp $(top_srcdir)/Trilinos_version.h $(DESTDIR)$(includedir)/. ;
	cp $(top_srcdir)/commonTools/buildTools/strip_dup_incl_paths.pl $(DESTDIR)$(includedir)/. ;
	cp $(top_srcdir)/commonTools/buildTools/strip_dup_libs.pl $(DESTDIR)$(includedir)/.

uninstall-hook:
	rm -f $(includedir)/Trilinos_version.h
	rm -f $(includedir)/strip_dup_incl_paths.pl
	rm -f $(includedir)/strip_dup_libs.pl


#
# The below targets allow you to use the new
# testharness to run the test suite as make targets
#

TRILINOS_HOME_DIR=@abs_top_srcdir@
TRILINOS_BUILD_DIR=@abs_top_builddir@
TRILINOS_MPI_MAX_PROC=4
TRILINOS_TEST_CATEGORY=INSTALL

runtests-serial :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=serial \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --category=$(TRILINOS_TEST_CATEGORY) \
  --output-dir=@abs_top_builddir@/runtests-results

runtests-mpi :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=mpi \
  --mpi-go=$(TRILINOS_MPI_GO) \
  --max-proc=$(TRILINOS_MPI_MAX_PROC) \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --category=$(TRILINOS_TEST_CATEGORY) \  --output-dir=@abs_top_builddir@/runtests-results

run-pytests: all-am
	@if [ -d packages/teuchos     ]; then cd packages/teuchos;     $(MAKE) run-pytests; fi
	@if [ -d packages/thyra       ]; then cd packages/thyra;       $(MAKE) run-pytests; fi
	@if [ -d packages/epetra      ]; then cd packages/epetra;      $(MAKE) run-pytests; fi
	@if [ -d packages/triutils    ]; then cd packages/triutils;    $(MAKE) run-pytests; fi
	@if [ -d packages/epetraext   ]; then cd packages/epetraext;   $(MAKE) run-pytests; fi
	@if [ -d packages/aztecoo     ]; then cd packages/aztecoo;     $(MAKE) run-pytests; fi
	@if [ -d packages/galeri      ]; then cd packages/galeri;      $(MAKE) run-pytests; fi
	@if [ -d packages/amesos      ]; then cd packages/amesos;      $(MAKE) run-pytests; fi
	@if [ -d packages/ifpack      ]; then cd packages/ifpack;      $(MAKE) run-pytests; fi
	@if [ -d packages/anasazi     ]; then cd packages/anasazi;     $(MAKE) run-pytests; fi
	@if [ -d packages/ml          ]; then cd packages/ml;          $(MAKE) run-pytests; fi
	@if [ -d packages/nox         ]; then cd packages/nox;         $(MAKE) run-pytests; fi
	@if [ -d packages/new_package ]; then cd packages/new_package; $(MAKE) run-pytests; fi

run-pyexamples: all-am
	@if [ -d packages/teuchos     ]; then cd packages/teuchos;     $(MAKE) run-pyexamples; fi
	@if [ -d packages/thyra       ]; then cd packages/thyra;       $(MAKE) run-pyexamples; fi
	@if [ -d packages/epetra      ]; then cd packages/epetra;      $(MAKE) run-pyexamples; fi
	@if [ -d packages/triutils    ]; then cd packages/triutils;    $(MAKE) run-pyexamples; fi
	@if [ -d packages/epetraext   ]; then cd packages/epetraext;   $(MAKE) run-pyexamples; fi
	@if [ -d packages/aztecoo     ]; then cd packages/aztecoo;     $(MAKE) run-pyexamples; fi
	@if [ -d packages/galeri      ]; then cd packages/galeri;      $(MAKE) run-pyexamples; fi
	@if [ -d packages/amesos      ]; then cd packages/amesos;      $(MAKE) run-pyexamples; fi
	@if [ -d packages/ifpack      ]; then cd packages/ifpack;      $(MAKE) run-pyexamples; fi
	@if [ -d packages/anasazi     ]; then cd packages/anasazi;     $(MAKE) run-pyexamples; fi
	@if [ -d packages/ml          ]; then cd packages/ml;          $(MAKE) run-pyexamples; fi
	@if [ -d packages/nox         ]; then cd packages/nox;         $(MAKE) run-pyexamples; fi
	@if [ -d packages/new_package ]; then cd packages/new_package; $(MAKE) run-pyexamples; fi

clean-python: all-am
	@if [ -d packages/PyTrilinos         ]; then cd packages/PyTrilinos;         $(MAKE) clean; fi
	@if [ -d packages/teuchos/python     ]; then cd packages/teuchos/python;     $(MAKE) clean; fi
	@if [ -d packages/thyra/python       ]; then cd packages/thyra/python;       $(MAKE) clean; fi
	@if [ -d packages/epetra/python      ]; then cd packages/epetra/python;      $(MAKE) clean; fi
	@if [ -d packages/triutils/python    ]; then cd packages/triutils/python;    $(MAKE) clean; fi
	@if [ -d packages/epetraext/python   ]; then cd packages/epetraext/python;   $(MAKE) clean; fi
	@if [ -d packages/aztecoo/python     ]; then cd packages/aztecoo/python;     $(MAKE) clean; fi
	@if [ -d packages/galeri/python      ]; then cd packages/galeri/python;      $(MAKE) clean; fi
	@if [ -d packages/amesos/python      ]; then cd packages/amesos/python;      $(MAKE) clean; fi
	@if [ -d packages/ifpack/python      ]; then cd packages/ifpack/python;      $(MAKE) clean; fi
	@if [ -d packages/anasazi/python     ]; then cd packages/anasazi/python;     $(MAKE) clean; fi
	@if [ -d packages/ml/python          ]; then cd packages/ml/python;          $(MAKE) clean; fi
	@if [ -d packages/nox/python         ]; then cd packages/nox/python;         $(MAKE) clean; fi
	@if [ -d packages/new_package/python ]; then cd packages/new_package/python; $(MAKE) clean; fi

#
# General targets to build tests and examples etc. since the default target
# now just builds the libraries.
#

examples:
	cd packages && $(MAKE) $@

install-examples:
	cd packages && $(MAKE) $@

clean-examples:
	cd packages && $(MAKE) $@

tests:
	cd packages && $(MAKE) $@

clean-tests:
	cd packages && $(MAKE) $@

everything:
	cd packages && $(MAKE) $@

install-everything:
	cd packages && $(MAKE) $@

clean-everything:
	cd packages && $(MAKE) $@

partial-clean:
	cd packages && $(MAKE) $@

