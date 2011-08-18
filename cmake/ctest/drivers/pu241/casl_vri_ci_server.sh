#!/bin/bash

# This script runs a CI server that builds the CASL VRI packages in a
# continuous integration loop and sends results to the CASL VRI CDash
# server.

#
# Set up
#

TRILINOS_BASE=/home/casl-vri-admin/CIDashboards/Trilinos
echo "TRILINOS_BASE = '$TRILINOS_BASE'"

# Must always start from this dir!
DASHBOARD_DRIVER_DIR=/home/casl-vri-admin/CIDashboards/Trilinos/cmake/ctest/drivers/pu241
cd $DASHBOARD_DRIVER_DIR

BASE_COMMAND=$TRILINOS_BASE/cmake/ctest/drivers/pu241/cron_driver_single_ci_iter_pu241.sh
LOOP_INTERVAL=1m
TODAY_RUN_TILL=23:59:00

echo "BASE_COMMAND = '$BASE_COMMAND'"
echo "LOOP_INTERVAL = '$LOOP_INTERVAL'"
echo "TODAY_RUN_TILL = '$TODAY_RUN_TILL'"

#BASE_COMMAND="echo mock_command"
#LOOP_INTERVAL=1s


#
# Executate
#

export PATH=/opt/trilinos-toolset/bin:$PATH

$TRILINOS_BASE/cmake/python/mailmsg.py "Starting CASL VRI CI testing server on pu241"

echo
echo "***"
echo "*** Starting CI server that tests changes that can affect CASL"
echo "***"
echo

# A) Run the first build all SS packages from scratch

env CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=OFF $BASE_COMMAND

# B) Run a CI loop that will terminate on time

# B.1) Define the inner loop command

CI_COMMAND="env TDD_FORCE_CMAKE_INSTALL=0 TDD_DO_SUBMIT=OFF CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=OFF CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON $BASE_COMMAND"
echo "CI_COMMAND = $CI_COMMAND"

# B.2) Run the loop

$TRILINOS_BASE/cmake/python/generic-looping-demon.py \
--command="$CI_COMMAND" \
--today-run-till=$TODAY_RUN_TILL \
--loop-interval=$LOOP_INTERVAL

$TRILINOS_BASE/cmake/python/mailmsg.py "Finished CASL VRI CI testing server on pu241"

# NOTES:
#
# On the first CI iteration, we let it update CMake/CTest and do the
# outer submit to the TDD dashboard.  This provides some detail on how
# the outter configure is being done and should provide some useful
# info.  However, for now, we turn off CDash submits to the outer TDD
# dashboard because currently it would submit worthless no-op results
# every LOOP_INTERVAL (e.g. 1 minute).  To be able to do the submit of
# the outer results we would have to make some tricky modifications to
# the CTest driver script so that the outter TDD driver could detect
# when at least one build case was done and if it was then do the
# submit in one shot at the very end.  I think you can do this but it
# would require some work.
#
# Above, we only install the offical version of CMake/CTest once on
# the first CI iteration.  We do not reinstall CMake over and over
# again to save time and there is little value in doing so;5B.
