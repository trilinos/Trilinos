# @HEADER
# ************************************************************************
# 
#	   Finite Element Interface to Linear Solvers (FEI)
#		 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# Questions? Contact Alan Williams (william@sandia.gov) 
# 
# ************************************************************************
# @HEADER

## #######################################################################
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log $(AUX_DIST) 

CORE_H = \
	$(top_builddir)/FEI_config.h \
	$(top_srcdir)/fei_base.hpp \
	$(top_srcdir)/fei_fwd.hpp

include_HEADERS = $(CORE_H)

TRILINOS_DIR=support-Trilinos

if HAVE_FEI_PETSC
PETSC_DIR=support-PETSc
else
PETSC_DIR=
endif

if SUB_EXAMPLE
EXAMPLES=example
else
EXAMPLES=
endif

if SUB_TEST
TESTS=test utest_cases utest \
	test/elemDOF_Aztec \
	test/remEssBCs \
	test/remEssBCs_Trilinos \
	test/cr_lagrange_20quad \
	test/cr_slave_20quad \
	test/cr_slave_20quad_old \
	test/cr_lagrange_20quad_AZLSC \
	test/cr_lagrange_20quad_old \
	test/cr_multifield_slave
else
TESTS=
endif

SUBDIRS=base $(TRILINOS_DIR) $(PETSC_DIR) test_utils $(TESTS) $(EXAMPLES)

EXTRA_DIST = \
	RELEASE_NOTES.txt

##
## The dist-hook target below gets run when 'make dist' is used,
## and deletes files that we don't want included in fei distributions.
##

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..
TRILINOS_MPI_MAX_PROC=4
TRILINOS_TEST_CATEGORY=INSTALL

DONT_DISTRIBUTE = bootstrap configure.ac FEI_config.h

if HAVE_MPI
FEI_CHECK_COMM=mpi
else
FEI_CHECK_COMM=serial
endif

if BUILD_TESTS
tests:
	@echo ""
	@echo "Now building fei tests."
	@echo ""
	cd $(top_builddir)/test && $(MAKE)
	cd $(top_builddir)/utest_cases && $(MAKE)
	cd $(top_builddir)/utest && $(MAKE)
	cd $(top_builddir)/test/elemDOF_Aztec && $(MAKE)
	cd $(top_builddir)/test/remEssBCs && $(MAKE)
	cd $(top_builddir)/test/remEssBCs_Trilinos && $(MAKE)
	cd $(top_builddir)/test/cr_lagrange_20quad && $(MAKE)
	cd $(top_builddir)/test/cr_slave_20quad && $(MAKE)
	cd $(top_builddir)/test/cr_slave_20quad_old && $(MAKE)
	cd $(top_builddir)/test/cr_lagrange_20quad_AZLSC && $(MAKE)
	cd $(top_builddir)/test/cr_lagrange_20quad_old && $(MAKE)
	cd $(top_builddir)/test/cr_multifield_slave && $(MAKE)
	@echo ""
	@echo "Finished building fei tests."
	@echo ""
else
tests:
	@echo "fei tests were disabled at configure time"
endif

if BUILD_EXAMPLES
examples:
	@echo ""
	@echo "Now building fei examples."
	@echo ""
	cd $(top_builddir)/example && $(MAKE)
	@echo ""
	@echo "Finished building fei examples."
	@echo ""

install-examples:
	cd $(top_builddir)/example && $(MAKE) install
else
examples:
	@echo "fei examples were disabled at configure time"

install-examples:
	@echo "fei examples were disabled at configure time"
endif

clean-tests:
	cd $(top_builddir)/test && $(MAKE) clean
	cd $(top_builddir)/utest_cases && $(MAKE) clean
	cd $(top_builddir)/utest && $(MAKE) clean
	cd $(top_builddir)/test/elemDOF_Aztec && $(MAKE) clean
	cd $(top_builddir)/test/remEssBCs && $(MAKE) clean
	cd $(top_builddir)/test/remEssBCs_Trilinos && $(MAKE) clean
	cd $(top_builddir)/test/cr_lagrange_20quad && $(MAKE) clean
	cd $(top_builddir)/test/cr_slave_20quad && $(MAKE) clean
	cd $(top_builddir)/test/cr_slave_20quad_old && $(MAKE) clean
	cd $(top_builddir)/test/cr_lagrange_20quad_AZLSC && $(MAKE) clean
	cd $(top_builddir)/test/cr_lagrange_20quad_old && $(MAKE) clean
	cd $(top_builddir)/test/cr_multifield_slave && $(MAKE) clean

clean-examples:
	cd $(top_builddir)/example && $(MAKE) clean

everything:
	$(MAKE) && $(MAKE) examples && $(MAKE) tests

clean-everything:
	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

install-everything:
	$(MAKE) install && $(MAKE) install-examples

dist-hook:
	for file in $(DONT_DISTRIBUTE) ; do \
		rm -f $(distdir)/$$file ; \
	done

run-examples:
	@echo
	@echo "=============================="
	@echo "Running FEI examples..."
	@echo "=============================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(FEI_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=FeiExamples \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=fei

run-tests:
	@echo
	@echo "==========================="
	@echo "Running FEI tests..."
	@echo "==========================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(FEI_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=FeiTests \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=fei

