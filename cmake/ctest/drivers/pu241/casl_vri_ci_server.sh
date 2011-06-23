#!/bin/bash


# This script runs a CI server that builds the CASL VRI packages in a
# continuous integration loop and sends results to the CASL VRI CDash
# server.

#
# Set up
#

# Override the pointer to Trilinos if you need to.
if ["$TRILINOS_BASE" == ""] ; then
  TRILINOS_BASE=/home/casl-vri-admin/Dashbaords/Trilinos
fi
echo "TRILINOS_BASE = '$TRILINOS_BASE'"

# Must always start from this dir!
DASHBOARD_DRIVER_DIR=/home/casl-vri-admin/Dashbaords/Trilinos/cmake/ctest/drivers/pu241
cd $DASHBOARD_DRIVER_DIR

BASE_COMMAND=$TRILINOS_BASE/cmake/ctest/drivers/pu241/cron_driver_single_ci_iter_pu241.sh
LOOP_INTERVAL=1m
TODAY_RUN_TILL=20:00:00

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

# ToDo: Set TDD_CONTINUOUS_INTEGRATION_MODE on the env and have it
# picked up in the CMakeLists.txt file and have it run the CI test script(s)!

echo
echo "***"
echo "*** Starting CI server that tests changes that can affect CASL"
echo "***"
echo

# Run the first build from scratch and build everything
env \
CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=ON \
CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=OFF \
$BASE_COMMAND

# Run a CI loop that will terminate on time

CI_COMMAND=env \
CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=OFF \
CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON \
$BASE_COMMAND

$TRILINOS_BASE/cmake/python/generic-looping-demon.py \
--command="$CI_COMMAND" \
--today-run-till=$TODAY_RUN_TILL \
--loop-interval=$LOOP_INTERVAL

$TRILINOS_BASE/cmake/python/mailmsg.py "Finished CASL VRI CI testing server on pu241"
