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
TEST_HOME=${HOME}/Trilinos/packages/ml/test/testharness/s861769

# -f  to configure/build/run
# -tf to run only
OPT=-f
cd ${TRILINOS_HOME}/testharness

# delete previous directories
/bin/rm -rf ${TRILINOS_HOME}/ml-test/
mkdir ${TRILINOS_HOME}/ml-test

TEST="basic"
StartTime=`date`
lamboot
for t in $TEST
do
  echo "========================="
  echo "Testing $t..."
  echo "========================="
  perl test-harness.plx $OPT ${TEST_HOME}/$t-input
  cp -r ${TRILINOS_HOME}/testharness/results ${TRILINOS_HOME}/ml-test/$t-results/
  echo 
  echo "Analyzing results:"
  echo "1) build error files (none is fine):"
  ls -1 ${TRILINOS_HOME}/ml-test/$t-results/ | grep -i error
  echo "2) test failed files (none is fine):"
  ls -1 ${TRILINOS_HOME}/ml-test/$t-results/ | grep -i failed
done
EndTime=`date`
echo "Test started on $StartTime,"
echo "ended on $EndTime"
lamhalt
cd -
