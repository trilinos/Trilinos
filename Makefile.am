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
	sampleScripts/cygwin_icl.pl sampleScripts/stratus_ar.pl \
	doc/Tutorial/README RELEASE_NOTES \
	doc/build_docs.pl \
	doc/check_installation.txt

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
SUBDIRS = packages

all-libs:
	cd $(top_builddir)/packages ; $(MAKE) all-libs

libs:
	cd $(top_builddir)/packages ; $(MAKE) libs

## #######################################################################
## Install the perl scrpts used for external makefiles
## #######################################################################

install-exec-hook:	
	cp $(top_srcdir)/commonTools/buildTools/strip_dup_incl_paths.pl $(exec_prefix)/include/. ;
	cp $(top_srcdir)/commonTools/buildTools/strip_dup_libs.pl $(exec_prefix)/include/.

uninstall-hook:
	rm -f $(exec_prefix)/include/strip_dup_incl_paths.pl
	rm -f $(exec_prefix)/include/strip_dup_libs.pl


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
  --category=$(TRILINOS_TEST_CATEGORY) \
  --output-dir=@abs_top_builddir@/runtests-results

run-examples: all-am
	@if [ -d packages/epetra    ]; then cd packages/epetra;    $(MAKE) run-examples; fi
	@if [ -d packages/epetraext ]; then cd packages/epetraext; $(MAKE) run-examples; fi
	@if [ -d packages/aztecoo   ]; then cd packages/aztecoo;   $(MAKE) run-examples; fi
	@if [ -d packages/ifpack    ]; then cd packages/ifpack;    $(MAKE) run-examples; fi
	@if [ -d packages/amesos    ]; then cd packages/amesos;    $(MAKE) run-examples; fi
	@if [ -d packages/ml        ]; then cd packages/ml;        $(MAKE) run-examples; fi
	@if [ -d packages/nox       ]; then cd packages/nox;       $(MAKE) run-examples; fi
	@if [ -d packages/didasko   ]; then cd packages/didasko;   $(MAKE) run-examples; fi

run-tests: all-am
	@if [ -d packages/epetra    ]; then cd packages/epetra;    $(MAKE) run-tests; fi
	@if [ -d packages/epetraext ]; then cd packages/epetraext; $(MAKE) run-tests; fi
	@if [ -d packages/aztecoo   ]; then cd packages/aztecoo;   $(MAKE) run-tests; fi
	@if [ -d packages/ifpack    ]; then cd packages/ifpack;    $(MAKE) run-tests; fi
	@if [ -d packages/amesos    ]; then cd packages/amesos;    $(MAKE) run-tests; fi
	@if [ -d packages/ml        ]; then cd packages/ml;        $(MAKE) run-tests; fi
	@if [ -d packages/nox       ]; then cd packages/nox;       $(MAKE) run-tests; fi
	@if [ -d packages/didasko   ]; then cd packages/didasko;   $(MAKE) run-tests; fi

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

examples:
	@if [ -d packages/teuchos     ]; then cd packages/teuchos;     $(MAKE) examples; fi
	@if [ -d packages/rtop        ]; then cd packages/rtop;        $(MAKE) examples; fi
	@if [ -d packages/kokkos      ]; then cd packages/kokkos;      $(MAKE) examples; fi
	@if [ -d packages/epetra      ]; then cd packages/epetra;      $(MAKE) examples; fi
	@if [ -d packages/triutils    ]; then cd packages/triutils;    $(MAKE) examples; fi
	@if [ -d packages/tpetra      ]; then cd packages/tpetra;      $(MAKE) examples; fi
	@if [ -d packages/epetraext   ]; then cd packages/epetraext;   $(MAKE) examples; fi
	@if [ -d packages/thyra       ]; then cd packages/thyra;       $(MAKE) examples; fi
	@if [ -d packages/isorropia   ]; then cd packages/isorropia;   $(MAKE) examples; fi
	@if [ -d packages/aztecoo     ]; then cd packages/aztecoo;     $(MAKE) examples; fi
	@if [ -d packages/galeri      ]; then cd packages/galeri;      $(MAKE) examples; fi
	@if [ -d packages/amesos      ]; then cd packages/amesos;      $(MAKE) examples; fi
	@if [ -d packages/ifpack      ]; then cd packages/ifpack;      $(MAKE) examples; fi
	@if [ -d packages/komplex     ]; then cd packages/komplex;     $(MAKE) examples; fi
	@if [ -d packages/claps       ]; then cd packages/claps;       $(MAKE) examples; fi
	@if [ -d packages/belos       ]; then cd packages/belos;       $(MAKE) examples; fi
	@if [ -d packages/anasazi     ]; then cd packages/anasazi;     $(MAKE) examples; fi
	@if [ -d packages/rbgen       ]; then cd packages/rbgen;       $(MAKE) examples; fi
	@if [ -d packages/pliris      ]; then cd packages/pliris;      $(MAKE) examples; fi
	@if [ -d packages/ml          ]; then cd packages/ml;          $(MAKE) examples; fi
	@if [ -d packages/aristos     ]; then cd packages/aristos;     $(MAKE) examples; fi
	@if [ -d packages/moertel     ]; then cd packages/moertel;     $(MAKE) examples; fi
	@if [ -d packages/stratimikos ]; then cd packages/stratimikos; $(MAKE) examples; fi
	@if [ -d packages/meros       ]; then cd packages/meros;       $(MAKE) examples; fi
	@if [ -d packages/sacado      ]; then cd packages/sacado;      $(MAKE) examples; fi
	@if [ -d packages/nox         ]; then cd packages/nox;         $(MAKE) examples; fi
	@if [ -d packages/capo        ]; then cd packages/capo;        $(MAKE) examples; fi
	@if [ -d packages/rythmos     ]; then cd packages/rythmos;     $(MAKE) examples; fi
	@if [ -d packages/trilinoscouplings ]; then cd packages/trilinoscouplings; $(MAKE) examples; fi
	@if [ -d packages/new_package ]; then cd packages/new_package; $(MAKE) examples; fi
	@if [ -d packages/moocho      ]; then cd packages/moocho;      $(MAKE) examples; fi
	@if [ -d packages/PyTrilinos  ]; then cd packages/PyTrilinos;  $(MAKE) examples; fi
	@if [ -d packages/didasko     ]; then cd packages/didasko;     $(MAKE) examples; fi
	@if [ -d packages/WebTrilinos ]; then cd packages/WebTrilinos; $(MAKE) examples; fi
	@if [ -d packages/external    ]; then cd packages/external;    $(MAKE) examples; fi

install-examples:
	@if [ -d packages/teuchos     ]; then cd packages/teuchos;     $(MAKE) install-examples; fi
	@if [ -d packages/rtop        ]; then cd packages/rtop;        $(MAKE) install-examples; fi
	@if [ -d packages/kokkos      ]; then cd packages/kokkos;      $(MAKE) install-examples; fi
	@if [ -d packages/epetra      ]; then cd packages/epetra;      $(MAKE) install-examples; fi
	@if [ -d packages/triutils    ]; then cd packages/triutils;    $(MAKE) install-examples; fi
	@if [ -d packages/tpetra      ]; then cd packages/tpetra;      $(MAKE) install-examples; fi
	@if [ -d packages/epetraext   ]; then cd packages/epetraext;   $(MAKE) install-examples; fi
	@if [ -d packages/thyra       ]; then cd packages/thyra;       $(MAKE) install-examples; fi
	@if [ -d packages/isorropia   ]; then cd packages/isorropia;   $(MAKE) install-examples; fi
	@if [ -d packages/aztecoo     ]; then cd packages/aztecoo;     $(MAKE) install-examples; fi
	@if [ -d packages/galeri      ]; then cd packages/galeri;      $(MAKE) install-examples; fi
	@if [ -d packages/amesos      ]; then cd packages/amesos;      $(MAKE) install-examples; fi
	@if [ -d packages/ifpack      ]; then cd packages/ifpack;      $(MAKE) install-examples; fi
	@if [ -d packages/komplex     ]; then cd packages/komplex;     $(MAKE) install-examples; fi
	@if [ -d packages/claps       ]; then cd packages/claps;       $(MAKE) install-examples; fi
	@if [ -d packages/belos       ]; then cd packages/belos;       $(MAKE) install-examples; fi
	@if [ -d packages/anasazi     ]; then cd packages/anasazi;     $(MAKE) install-examples; fi
	@if [ -d packages/rbgen       ]; then cd packages/rbgen;       $(MAKE) install-examples; fi
	@if [ -d packages/pliris      ]; then cd packages/pliris;      $(MAKE) install-examples; fi
	@if [ -d packages/ml          ]; then cd packages/ml;          $(MAKE) install-examples; fi
	@if [ -d packages/aristos     ]; then cd packages/aristos;     $(MAKE) install-examples; fi
	@if [ -d packages/moertel     ]; then cd packages/moertel;     $(MAKE) install-examples; fi
	@if [ -d packages/stratimikos ]; then cd packages/stratimikos; $(MAKE) install-examples; fi
	@if [ -d packages/meros       ]; then cd packages/meros;       $(MAKE) install-examples; fi
	@if [ -d packages/sacado      ]; then cd packages/sacado;      $(MAKE) install-examples; fi
	@if [ -d packages/nox         ]; then cd packages/nox;         $(MAKE) install-examples; fi
	@if [ -d packages/capo        ]; then cd packages/capo;        $(MAKE) install-examples; fi
	@if [ -d packages/rythmos     ]; then cd packages/rythmos;     $(MAKE) install-examples; fi
	@if [ -d packages/trilinoscouplings ]; then cd packages/trilinoscouplings; $(MAKE) install-examples; fi
	@if [ -d packages/new_package ]; then cd packages/new_package; $(MAKE) install-examples; fi
	@if [ -d packages/moocho      ]; then cd packages/moocho;      $(MAKE) install-examples; fi
	@if [ -d packages/PyTrilinos  ]; then cd packages/PyTrilinos;  $(MAKE) install-examples; fi
	@if [ -d packages/didasko     ]; then cd packages/didasko;     $(MAKE) install-examples; fi
	@if [ -d packages/WebTrilinos ]; then cd packages/WebTrilinos; $(MAKE) install-examples; fi
	@if [ -d packages/external    ]; then cd packages/external;    $(MAKE) install-examples; fi

clean-examples:
	@if [ -d packages/teuchos     ]; then cd packages/teuchos;     $(MAKE) clean-examples; fi
	@if [ -d packages/rtop        ]; then cd packages/rtop;        $(MAKE) clean-examples; fi
	@if [ -d packages/kokkos      ]; then cd packages/kokkos;      $(MAKE) clean-examples; fi
	@if [ -d packages/epetra      ]; then cd packages/epetra;      $(MAKE) clean-examples; fi
	@if [ -d packages/triutils    ]; then cd packages/triutils;    $(MAKE) clean-examples; fi
	@if [ -d packages/tpetra      ]; then cd packages/tpetra;      $(MAKE) clean-examples; fi
	@if [ -d packages/epetraext   ]; then cd packages/epetraext;   $(MAKE) clean-examples; fi
	@if [ -d packages/thyra       ]; then cd packages/thyra;       $(MAKE) clean-examples; fi
	@if [ -d packages/isorropia   ]; then cd packages/isorropia;   $(MAKE) clean-examples; fi
	@if [ -d packages/aztecoo     ]; then cd packages/aztecoo;     $(MAKE) clean-examples; fi
	@if [ -d packages/galeri      ]; then cd packages/galeri;      $(MAKE) clean-examples; fi
	@if [ -d packages/amesos      ]; then cd packages/amesos;      $(MAKE) clean-examples; fi
	@if [ -d packages/ifpack      ]; then cd packages/ifpack;      $(MAKE) clean-examples; fi
	@if [ -d packages/komplex     ]; then cd packages/komplex;     $(MAKE) clean-examples; fi
	@if [ -d packages/claps       ]; then cd packages/claps;       $(MAKE) clean-examples; fi
	@if [ -d packages/belos       ]; then cd packages/belos;       $(MAKE) clean-examples; fi
	@if [ -d packages/anasazi     ]; then cd packages/anasazi;     $(MAKE) clean-examples; fi
	@if [ -d packages/rbgen       ]; then cd packages/rbgen;       $(MAKE) clean-examples; fi
	@if [ -d packages/pliris      ]; then cd packages/pliris;      $(MAKE) clean-examples; fi
	@if [ -d packages/ml          ]; then cd packages/ml;          $(MAKE) clean-examples; fi
	@if [ -d packages/aristos     ]; then cd packages/aristos;     $(MAKE) clean-examples; fi
	@if [ -d packages/moertel     ]; then cd packages/moertel;     $(MAKE) clean-examples; fi
	@if [ -d packages/stratimikos ]; then cd packages/stratimikos; $(MAKE) clean-examples; fi
	@if [ -d packages/meros       ]; then cd packages/meros;       $(MAKE) clean-examples; fi
	@if [ -d packages/sacado      ]; then cd packages/sacado;      $(MAKE) clean-examples; fi
	@if [ -d packages/nox         ]; then cd packages/nox;         $(MAKE) clean-examples; fi
	@if [ -d packages/capo        ]; then cd packages/capo;        $(MAKE) clean-examples; fi
	@if [ -d packages/rythmos     ]; then cd packages/rythmos;     $(MAKE) clean-examples; fi
	@if [ -d packages/trilinoscouplings ]; then cd packages/trilinoscouplings; $(MAKE) clean-examples; fi
	@if [ -d packages/new_package ]; then cd packages/new_package; $(MAKE) clean-examples; fi
	@if [ -d packages/moocho      ]; then cd packages/moocho;      $(MAKE) clean-examples; fi
	@if [ -d packages/PyTrilinos  ]; then cd packages/PyTrilinos;  $(MAKE) clean-examples; fi
	@if [ -d packages/didasko     ]; then cd packages/didasko;     $(MAKE) clean-examples; fi
	@if [ -d packages/WebTrilinos ]; then cd packages/WebTrilinos; $(MAKE) clean-examples; fi
	@if [ -d packages/external    ]; then cd packages/external;    $(MAKE) clean-examples; fi

tests:
	@if [ -d packages/teuchos     ]; then cd packages/teuchos;     $(MAKE) tests; fi
	@if [ -d packages/rtop        ]; then cd packages/rtop;        $(MAKE) tests; fi
	@if [ -d packages/kokkos      ]; then cd packages/kokkos;      $(MAKE) tests; fi
	@if [ -d packages/epetra      ]; then cd packages/epetra;      $(MAKE) tests; fi
	@if [ -d packages/triutils    ]; then cd packages/triutils;    $(MAKE) tests; fi
	@if [ -d packages/tpetra      ]; then cd packages/tpetra;      $(MAKE) tests; fi
	@if [ -d packages/epetraext   ]; then cd packages/epetraext;   $(MAKE) tests; fi
	@if [ -d packages/thyra       ]; then cd packages/thyra;       $(MAKE) tests; fi
	@if [ -d packages/isorropia   ]; then cd packages/isorropia;   $(MAKE) tests; fi
	@if [ -d packages/aztecoo     ]; then cd packages/aztecoo;     $(MAKE) tests; fi
	@if [ -d packages/galeri      ]; then cd packages/galeri;      $(MAKE) tests; fi
	@if [ -d packages/amesos      ]; then cd packages/amesos;      $(MAKE) tests; fi
	@if [ -d packages/ifpack      ]; then cd packages/ifpack;      $(MAKE) tests; fi
	@if [ -d packages/komplex     ]; then cd packages/komplex;     $(MAKE) tests; fi
	@if [ -d packages/claps       ]; then cd packages/claps;       $(MAKE) tests; fi
	@if [ -d packages/belos       ]; then cd packages/belos;       $(MAKE) tests; fi
	@if [ -d packages/anasazi     ]; then cd packages/anasazi;     $(MAKE) tests; fi
	@if [ -d packages/rbgen       ]; then cd packages/rbgen;       $(MAKE) tests; fi
	@if [ -d packages/pliris      ]; then cd packages/pliris;      $(MAKE) tests; fi
	@if [ -d packages/ml          ]; then cd packages/ml;          $(MAKE) tests; fi
	@if [ -d packages/aristos     ]; then cd packages/aristos;     $(MAKE) tests; fi
	@if [ -d packages/moertel     ]; then cd packages/moertel;     $(MAKE) tests; fi
	@if [ -d packages/stratimikos ]; then cd packages/stratimikos; $(MAKE) tests; fi
	@if [ -d packages/meros       ]; then cd packages/meros;       $(MAKE) tests; fi
	@if [ -d packages/sacado      ]; then cd packages/sacado;      $(MAKE) tests; fi
	@if [ -d packages/nox         ]; then cd packages/nox;         $(MAKE) tests; fi
	@if [ -d packages/capo        ]; then cd packages/capo;        $(MAKE) tests; fi
	@if [ -d packages/rythmos     ]; then cd packages/rythmos;     $(MAKE) tests; fi
	@if [ -d packages/trilinoscouplings ]; then cd packages/trilinoscouplings; $(MAKE) tests; fi
	@if [ -d packages/new_package ]; then cd packages/new_package; $(MAKE) tests; fi
	@if [ -d packages/moocho      ]; then cd packages/moocho;      $(MAKE) tests; fi
	@if [ -d packages/PyTrilinos  ]; then cd packages/PyTrilinos;  $(MAKE) tests; fi
	@if [ -d packages/didasko     ]; then cd packages/didasko;     $(MAKE) tests; fi
	@if [ -d packages/WebTrilinos ]; then cd packages/WebTrilinos; $(MAKE) tests; fi
	@if [ -d packages/external    ]; then cd packages/external;    $(MAKE) tests; fi

clean-tests:
	@if [ -d packages/teuchos     ]; then cd packages/teuchos;     $(MAKE) clean-tests; fi
	@if [ -d packages/rtop        ]; then cd packages/rtop;        $(MAKE) clean-tests; fi
	@if [ -d packages/kokkos      ]; then cd packages/kokkos;      $(MAKE) clean-tests; fi
	@if [ -d packages/epetra      ]; then cd packages/epetra;      $(MAKE) clean-tests; fi
	@if [ -d packages/triutils    ]; then cd packages/triutils;    $(MAKE) clean-tests; fi
	@if [ -d packages/tpetra      ]; then cd packages/tpetra;      $(MAKE) clean-tests; fi
	@if [ -d packages/epetraext   ]; then cd packages/epetraext;   $(MAKE) clean-tests; fi
	@if [ -d packages/thyra       ]; then cd packages/thyra;       $(MAKE) clean-tests; fi
	@if [ -d packages/isorropia   ]; then cd packages/isorropia;   $(MAKE) clean-tests; fi
	@if [ -d packages/aztecoo     ]; then cd packages/aztecoo;     $(MAKE) clean-tests; fi
	@if [ -d packages/galeri      ]; then cd packages/galeri;      $(MAKE) clean-tests; fi
	@if [ -d packages/amesos      ]; then cd packages/amesos;      $(MAKE) clean-tests; fi
	@if [ -d packages/ifpack      ]; then cd packages/ifpack;      $(MAKE) clean-tests; fi
	@if [ -d packages/komplex     ]; then cd packages/komplex;     $(MAKE) clean-tests; fi
	@if [ -d packages/claps       ]; then cd packages/claps;       $(MAKE) clean-tests; fi
	@if [ -d packages/belos       ]; then cd packages/belos;       $(MAKE) clean-tests; fi
	@if [ -d packages/anasazi     ]; then cd packages/anasazi;     $(MAKE) clean-tests; fi
	@if [ -d packages/rbgen       ]; then cd packages/rbgen;       $(MAKE) clean-tests; fi
	@if [ -d packages/pliris      ]; then cd packages/pliris;      $(MAKE) clean-tests; fi
	@if [ -d packages/ml          ]; then cd packages/ml;          $(MAKE) clean-tests; fi
	@if [ -d packages/aristos     ]; then cd packages/aristos;     $(MAKE) clean-tests; fi
	@if [ -d packages/moertel     ]; then cd packages/moertel;     $(MAKE) clean-tests; fi
	@if [ -d packages/stratimikos ]; then cd packages/stratimikos; $(MAKE) clean-tests; fi
	@if [ -d packages/meros       ]; then cd packages/meros;       $(MAKE) clean-tests; fi
	@if [ -d packages/sacado      ]; then cd packages/sacado;      $(MAKE) clean-tests; fi
	@if [ -d packages/nox         ]; then cd packages/nox;         $(MAKE) clean-tests; fi
	@if [ -d packages/capo        ]; then cd packages/capo;        $(MAKE) clean-tests; fi
	@if [ -d packages/rythmos     ]; then cd packages/rythmos;     $(MAKE) clean-tests; fi
	@if [ -d packages/trilinoscouplings ]; then cd packages/trilinoscouplings; $(MAKE) clean-tests; fi
	@if [ -d packages/new_package ]; then cd packages/new_package; $(MAKE) clean-tests; fi
	@if [ -d packages/moocho      ]; then cd packages/moocho;      $(MAKE) clean-tests; fi
	@if [ -d packages/PyTrilinos  ]; then cd packages/PyTrilinos;  $(MAKE) clean-tests; fi
	@if [ -d packages/didasko     ]; then cd packages/didasko;     $(MAKE) clean-tests; fi
	@if [ -d packages/WebTrilinos ]; then cd packages/WebTrilinos; $(MAKE) clean-tests; fi
	@if [ -d packages/external    ]; then cd packages/external;    $(MAKE) clean-tests; fi

everything:
	 $(MAKE) ; $(MAKE) examples ; $(MAKE) tests

install-everything:
	$(MAKE) install; $(MAKE) install-examples

clean-everything:
	$(MAKE) clean-examples ; $(MAKE) clean-tests ; $(MAKE) clean

