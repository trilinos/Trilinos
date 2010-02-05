#!/bin/bash

#
# This is the script that I use to checkin to Trilinos on s903186.sandia.gov.
# You can copy this script and adapt it to your own machine.
#
# Options to run with:
#
#  For the extra builds pass in (some subset of):
#
#    --extra-builds=SERIAL_DEBUG_BOOST_TRACE,SERIAL_DEBUG_TRACE,SERIAL_DEBUG_BOOST,MPI_DEBUG_INT,MPI_OPT
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
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTPL_ENABLE_MPI:BOOL=ON
-DMPI_BASE_DIR:PATH=/Users/jmwille/install
-DTPL_ENABLE_Boost:BOOL=ON
-DBoost_INCLUDE_DIRS:PATH=/Users/jmwille/install/boost_1_41_0
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DTeuchos_ENABLE_DEBUG_RCP_NODE_TRACING:BOOL=ON
" > MPI_DEBUG_BOOST_TRACE.config

echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTPL_ENABLE_MPI:BOOL=ON
-DMPI_BASE_DIR:PATH=/Users/jmwille/install
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DTeuchos_ORDINAL_TYPE:STRIRNG=int
" > MPI_DEBUG_INT.config

echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTPL_ENABLE_MPI:BOOL=ON
-DMPI_BASE_DIR:PATH=/Users/jmwille/install
-DTrilinos_ENABLE_DEBUG:BOOL=OFF \
-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF \
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF \
" > MPI_OPT.config

echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTPL_ENABLE_Boost:BOOL=ON
-DBoost_INCLUDE_DIRS:PATH=/Users/jmwille/install/boost_1_41_0
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DTeuchos_ENABLE_DEBUG_RCP_NODE_TRACING:BOOL=ON
-DCMAKE_CXX_COMPILER:FILEPATH=/Users/jmwille/install/bin/g++
-DCMAKE_C_COMPILER:FILEPATH=/Users/jmwille/install/bin/gcc
-DCMAKE_Fortran_COMPILER:FILEPATH=/Users/jmwille/install/bin/gfortran
" > SERIAL_DEBUG_BOOST_TRACE.config

echo "
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DTeuchos_ENABLE_DEBUG_RCP_NODE_TRACING:BOOL=ON
-DCMAKE_CXX_COMPILER:FILEPATH=/Users/jmwille/install/bin/g++ 
-DCMAKE_C_COMPILER:FILEPATH=/Users/jmwille/install/bin/gcc 
-DCMAKE_Fortran_COMPILER:FILEPATH=/Users/jmwille/install/bin/gfortran 
" > SERIAL_DEBUG_TRACE.config

echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=ON
-DTPL_ENABLE_Boost:BOOL=ON
-DBoost_INCLUDE_DIRS:PATH=/Users/jmwille/install/boost_1_41_0
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-DCMAKE_CXX_COMPILER:FILEPATH=/Users/jmwille/install/bin/g++ 
-DCMAKE_C_COMPILER:FILEPATH=/Users/jmwille/install/bin/gcc 
-DCMAKE_Fortran_COMPILER:FILEPATH=/Users/jmwille/install/bin/gfortran 
" > SERIAL_DEBUG_BOOST.config


#
# Run the standard checkin testing script with my specializations
#

../../Trilinos/checkin-test.py \
--make-options=-j8 \
--ctest-options=-j2 \
--ctest-timeout=300 \
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
