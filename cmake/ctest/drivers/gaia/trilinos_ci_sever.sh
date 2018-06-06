#!/bin/bash
#
# This is the main driver script for a poor-man's CI server for the Trilinos
# pre-push CI build.
#

DRIVER_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
echo "DRIVER_SCRIPT_DIR = '$DRIVER_SCRIPT_DIR'"

#TRILINOS_DIR=`readlink -f ${DRIVER_SCRIPT_DIR}/../../../..`
TRILINOS_DIR=${DRIVER_SCRIPT_DIR}/../../../..
echo "TRILINOS_DIR='${TRILINOS_DIR}'"

BASE_COMMAND=$DRIVER_SCRIPT_DIR/single_ci_iter.sh
LOOP_INTERVAL=3m
TODAY_RUN_TILL=18:00:00
PAUSE_FILE=$PWD/pause_ci_server.txt

# Un-comment to test out locally
#BASE_COMMAND="echo mock_command"
#LOOP_INTERVAL=1s

echo "BASE_COMMAND = '$BASE_COMMAND'"
echo "LOOP_INTERVAL = '$LOOP_INTERVAL'"
echo "TODAY_RUN_TILL = '$TODAY_RUN_TILL'"
echo "PAUSE_FILE = '$PAUSE_FILE'"

#
# Execute
#

# Have to source this or you don't even get the module command!
source /etc/bashrc

HOSTNAME=`hostname`
MAILMSG=$TRILINOS_DIR/cmake/tribits/python_utils/mailmsg.py

if [ "$TRILINOS_SKIP_CI_EMAILS" == "" ] ; then
  $MAILMSG "Starting Trilinos standrad CI testing server on '$HOSTNAME'"
fi

# A) Run the first build starting from scratch

env CI_FIRST_ITERATION=1 \
  CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE \
  CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=OFF \
  $BASE_COMMAND

# B) Run a CI loop that will terminate on time

# B.1) Define the inner loop command

CI_COMMAND="env CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=OFF CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON $BASE_COMMAND"
echo "CI_COMMAND = $CI_COMMAND"

# B.2) Run the CI loop

if [ "$TRILINOS_SKIP_CI_ITERATION" == "" ] ; then
  $TRILINOS_DIR/cmake/tribits/python_utils/generic-looping-demon.py \
  --command="$CI_COMMAND" \
  --today-run-till=$TODAY_RUN_TILL \
  --loop-interval=$LOOP_INTERVAL \
  --pause-file=$PAUSE_FILE
else
  echo "Skipping CI iterations because TRILINOS_SKIP_CI_ITERATION = '$TRILINOS_SKIP_CI_ITERATION' != ''"
fi

if [ "$TRILINOS_SKIP_CI_EMAILS" == "" ] ; then
  $MAILMSG "Ending Trilinos standrad CI testing server on '$HOSTNAME'"
fi
