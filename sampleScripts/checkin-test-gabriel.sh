#!/bin/bash

#
# This is the script that I used to checkin to Trilinos on gabriel.sandia.gov.
# You can copy this script and adapt it to your own machine.
#

EXTRA_ARGS=$@

if [ "$TRILINOS_HOME" == "" ] ; then
  TRILINOS_HOME=/home/rabartl/PROJECTS/Trilinos.base/Trilinos
fi

echo "-DBUILD_SHARED_LIBS:BOOL=ON" > COMMON.config

$TRILINOS_HOME/cmake/python/checkin-test.py \
--make-options="-j4" \
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS  
