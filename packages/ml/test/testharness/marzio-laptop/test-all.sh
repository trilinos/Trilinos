#!/bin/bash
# Simple scripts that executes all the ML tests, as defined in the ML/test/testharness
# subdirectory. This file is supposed to be executed as a cron job.
#
# This file jumps into the Trilinos test harness, and executes the perl script
# with a variety of different options. 
#
# MS, last modified on Nov-21
#
TRILINOS_HOME=${HOME}/Trilinos
TEST_HOME=${HOME}/Trilinos/packages/ml/test/testharness/marzio-laptop

# -f  to configure/build/run
# -tf to run only
OPT=-f
cd ${TRILINOS_HOME}/testharness

# standalone    ML alone (nothing else)
# epetra        ML with epetra, but without AztecOO
# basic         only basic trilinos packages (no Teuchos)
# trilinos      this should be the basic usage of ML within Trilinos
# parmetis      Epetra + Teuchos + ParMETIS (only MPI)
# all           My script, with all I have on my machien (only MPI)
TEST="standalone epetra basic trilinos parmetis all"
StartTime=`date`
for t in $TEST
do
  echo "--------------------------------"
  echo "Testing $t..."
  echo "--------------------------------"
  perl test-harness.plx $OPT ${TEST_HOME}/$t-input
  /bin/rm -rf ${TRILINOS_HOME}/$t-results
  cp -r ${TRILINOS_HOME}/testharness/results ${TRILINOS_HOME}/$t-results/
  echo 
  echo "Analyzing results:"
  echo "1) build error files:"
  ls -1 $t-results/ | grep -i error
  echo "2) test failed files:"
  ls -1 $t-results/ | grep -i failed
done
EndTime=`date`
echo "Test started on $StartTime,"
echo "ended on $EndTime"
