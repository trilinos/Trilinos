#!/bin/bash
# Simple scripts that executes all the ML tests, as defined in the ML/test/testharness
# subdirectory. This file is supposed to be executed as a cron job.
#
# This file jumps into the Trilinos test harness, and executes the perl script
# with a variety of different options. Note that the test-XXX files differ only
# in the specification of the elements files, and the directory names.
#
# WARNING: still to be completed...
#
# MS, last modified on Nov-19
#
TRILINOS_HOME=${HOME}/Trilinos
cd ${TRILINOS_HOME}/testharness
#
# ML as a stand-alone package
#
#perl test-harness.plx -f ../packages/ml/test/testharness/test-ml-standalone
#
# ML with basic Trilinos packages (Epetra, AztecOO, triutils)
#
#perl test-harness.plx -f ../packages/ml/test/testharness/test-ml-basic
#
# This is the standard Trilinos usage of ML (with Teuchos)
#
perl test-harness.plx -f ../packages/ml/test/testharness/test-ml-trilinos
#
# as ml-trilinos, with METIS and ParMETIS support
#
perl test-harness.plx -f ../packages/ml/test/testharness/test-ml-trilinos-parmetis
