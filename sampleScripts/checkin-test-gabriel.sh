#!/bin/bash

#
# This is the script that I use to checkin to Trilinos on gabriel.sandia.gov.
# You can copy this script and adapt it to your own machine.
#
# Options to run with:
#
#  For the extra builds pass in (some subset of):
#
#    --extra-builds=SERIAL_DEBUG_BOOST_TRACE,SERIAL_DEBUG_TRACE,SERIAL_DEBUG_BOOST,MPI_DEBUG_INT
#
# If you want to automatically invoke a remote pull/test/push on godel, you can
# use the argument:
#
#    "--execute-on-ready-to-push=\"ssh -q godel '~/PROJECTS/Trilinos.base.checkin/checkin-test-godel-remote.sh --extra-pull-from gabriel:master' &\""
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
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTPL_ENABLE_MPI:BOOL=ON
-DTPL_ENABLE_Boost:BOOL=ON
-DBoost_INCLUDE_DIRS:PATH=$HOME/PROJECTS/install/boost-1.40.0/include
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DTeuchos_ENABLE_DEBUG_RCP_NODE_TRACING:BOOL=ON
" > MPI_DEBUG_BOOST_TRACE.config

echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTPL_ENABLE_MPI:BOOL=ON
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DTeuchos_ORDINAL_TYPE:STRIRNG=int
" > MPI_DEBUG_INT.config

echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTPL_ENABLE_Boost:BOOL=ON
-DBoost_INCLUDE_DIRS:PATH=$HOME/PROJECTS/install/boost-1.40.0/include
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DTeuchos_ENABLE_DEBUG_RCP_NODE_TRACING:BOOL=ON
" > SERIAL_DEBUG_BOOST_TRACE.config

echo "
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DTeuchos_ENABLE_DEBUG_RCP_NODE_TRACING:BOOL=ON
" > SERIAL_DEBUG_TRACE.config

echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTPL_ENABLE_Boost:BOOL=ON
-DBoost_INCLUDE_DIRS:PATH=$HOME/PROJECTS/install/boost-1.40.0/include
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
" > SERIAL_DEBUG_BOOST.config


#
# Run the standard checkin testing script with my specializations
#

../../Trilinos/checkin-test.py \
--make-options=-j8 \
--ctest-options=-j8 \
--ctest-timeout=180 \
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS  

#
# NOTES:
#
# (*) Enabling shared libaries makes relinks go *much* faster and massively
# speeds up the checkin testing process when rebulding!
#
# (*) Be sure to set --make-options="-jN" to speed up building.  It makes a
# *big* difference.
#
# (*) Passing -jN to ctest with --ctest-options="-jN" can speed up running
# the tests but I have not seen very good speedup in general and some people
# have reported no speedup at all.  You should experiment with ctest -jN to
# see what number N works well on your machine.
#
