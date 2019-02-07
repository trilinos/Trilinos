#!/bin/bash

#
# Gives the CDash testing day YYYY-MM-DD that that correpsonds to a CDash
# build start time right now.  This takes into account the current Trilinos
# CDash project testing day start time (in UTC).
#
# NOTE: The argument --cdash-project-start-time below must match what is set
# in the Trilinos CDash project settings!
#

# Get the location of the base Trilinos directory and tribits/ci_support/

if [ "$ATDM_TRILINOS_DIR" == "" ] ; then
  # Grab from the symlink (only works on Linux)
  _ABS_FILE_PATH=`readlink -f $0` || \
   echo "Could not follow symlink to set TRILINOS_DIR!"
  if [ "$_ABS_FILE_PATH" != "" ] ; then
    export STD_ATDM_DIR=`dirname $_ABS_FILE_PATH`
    export ATDM_TRILINOS_DIR=`readlink -f $STD_ATDM_DIR/../../..`
  fi
fi

TRIBITS_CI_SUPPORT_DIR="${ATDM_TRILINOS_DIR}/cmake/tribits/ci_support"

# Run the script and return the CDash testing day YYYY-MM-DD to the STDOUT

${TRIBITS_CI_SUPPORT_DIR}/cdash_build_testing_date.py \
  --cdash-project-start-time="04:01"  # UTC
