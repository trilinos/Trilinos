#!/bin/bash

#
# This is the script that I use to checkin to Trilinos on s903186.sandia.gov.
# You can copy this script and adapt it to your own machine.
#
# Options to run with:
#
# To automatically invoke a remote pull/test/push on trilinos-test2, you can
# use the argument:
#
#    "--execute-on-ready-to-push=\"ssh -q trilinos-test2 '~/Trilinos-remote-testing/checkin-test-godel-remote.sh --extra-pull-from s903186:master' &\""
#
# NOTE: You will need the funny quotes when passing through a shell script.
# You also need to leave the '=' off of the --extra-pull-from argument because
# the python scripting code can't handle an '=' inside of an quoted argument
# (yet).
#

#
# Allow command-line arguments to pass through to cmake configure!
#

EXTRA_ARGS=$@

#
# Set up configuration files
#

echo "-DBUILD_SHARED_LIBS:BOOL=ON" > COMMON.config

echo "
-DMPI_BASE_DIR:PATH=/Users/jmwille/install
" > MPI_DEBUG.config

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
" > SERIAL_DEBUG.config


#
# Run the standard checkin testing script with my specializations
#

../../Trilinos/checkin-test.py \
-j8 \
--ctest-timeout=300 \
$EXTRA_ARGS  
