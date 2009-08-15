#!/bin/bash

#
# This is the script that I used to checkin to Trilinos on gabriel.sandia.gov.
# You can copy this script and adapt it to your own machine.
#

EXTRA_ARGS=$@

if [ "$TRILINOS_HOME" == "" ] ; then
  TRILINOS_HOME=/home/rabartl/PROJECTS/Trilinos.base/Trilinos
fi

echo "-DBUILD_SHARED:BOOL=ON" > COMMON.config

$TRILINOS_HOME/cmake/python/checkin-test.py \
--make-options="-j4" \
--ctest-options="-W 100 -j4" \
$EXTRA_ARGS  
