#!/bin/tcsh
# Simple scripts that executes all the ML tests, as defined in the ml/test/testharness
# subdirectory. This file is supposed to be executed as a cron job.
#
# This file jumps into the Trilinos test harness, and executes the perl script
# with a variety of different options. 
#
# Modification of Marzio Sala's script.
#
#  Current revision: $Revision$
#  Branch:           $Branch$
#  Last modified:    $Date$
#  Modified by:      $Author$
#
set TRILINOS_HOME="${HOME}/Trilinos/development-branch/Trilinos"
set TEST_HOME="${TRILINOS_HOME}/packages/ml/test/testharness/jhu-aphrodite"
set TEST_SAVE="${TRILINOS_HOME}/ml-test"

# -f  to configure/build/run
# -tf to run only
set OPT=-f
cd ${TRILINOS_HOME}/testharness

#set TEST="basic teuchos partition"
set TEST="partition"
set StartTime=`date`
echo ""
echo "*** Tests to run: ${TEST}"
/bin/rm -rf ${TRILINOS_HOME}/ml-test/
echo "*** Results will be copied into ${TEST_SAVE}"
mkdir ${TEST_SAVE}
foreach t ($TEST)
  echo "========================="
  echo "Testing $t..."
  echo "========================="
  perl test-harness.plx $OPT ${TEST_HOME}/$t-input
  cp -rp ${TRILINOS_HOME}/testharness/results ${TEST_SAVE}/$t-results/
  echo 
  echo "Analyzing results:"
  echo "1) build error files (none is fine):"
  ls -1 ${TEST_SAVE}/$t-results/ | grep -i error
  echo "2) test failed files (none is fine):"
  ls -1 ${TEST_SAVE}/$t-results/ | grep -i failed
end
set EndTime=`date`
echo ""
echo "Test started on $StartTime"
echo "Test ended on   $EndTime"
