#!/bin/bash

#
# Test ATDM configurations of Trilinos
#
# Usage:
#
#   ./checkin-test-atdm.sh \
#     <job-name-keys0> <job-name-keys1> ... <job-name-keysn> \
#     [other checkin-test options]
#
# If just 'all' is passsed in for the <job-name-keys> list, then the list of
# all of of the supported jobs for that current system will be loaded from
# <system_name>/all_supported_builds.sh.
#
# To use this script, just symlink this into any desired directory where to
# run from like:
#
#   cd <some-base-dir>/
#   ln -s $TRILNOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
#
# then run it as above.  For example, to locally test a few builds for just
# Kokkos use:
#
#   ./checkin-test-atdm.sh \
#     gnu-opt-openmp intel-debug-serial \
#     --no-enable-fwd-packages --enable-packages=Kokkos \
#     --local-do-all
#
# This will only send email for the final check of all of the builds
# specified.
#

echo
echo "***"
echo "*** $0 " "$@"
echo "***"
echo

if [ "$ATDM_TRILINOS_DIR" == "" ] ; then
  # Grab from the symlink (only works on Linux)
  _ABS_FILE_PATH=`readlink -f $0` || \
   echo "Could not follow symlink to set TRILINOS_DIR!"
  if [ "$_ABS_FILE_PATH" != "" ] ; then
    export STD_ATDM_DIR=`dirname $_ABS_FILE_PATH`
    export ATDM_TRILINOS_DIR=$STD_ATDM_DIR/../../..
  fi
fi

echo "ATDM_TRILINOS_DIR = '$ATDM_TRILINOS_DIR'"

if [ "$ATDM_TRILINOS_DIR" == "" ] ; then
  echo "ERROR: Cannot determine TRILINOS_DIR (you must be on a non-Linux system or you must have copied the script instead of symlinking it as per instructions)."
  exit 1
fi

echo
source $STD_ATDM_DIR/utils/get_known_system_name.sh

#
# A) Parse the arguments
#

# A.1) Pull off the initial <job-name-keysi> arguments

ATDM_JOB_NAME_KEYS_LIST=
while [[ ! "$1" == "--"* ]] && [[ ! "$1" == "" ]] ; do
  if [[ "$ATDM_JOB_NAME_KEYS_LIST" == "" ]] ; then
    ATDM_JOB_NAME_KEYS_LIST="$1"
  else
    ATDM_JOB_NAME_KEYS_LIST="$ATDM_JOB_NAME_KEYS_LIST $1"
  fi
  shift
done

if [[ "$ATDM_JOB_NAME_KEYS_LIST" == "" ]] ; then
  echo "Error, at least one <job-name-keys> (e.g. gnu-opt-openmp) argument is required!"
  exit 1
fi

if [[ "$ATDM_JOB_NAME_KEYS_LIST" == "all" ]] ; then
  export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=
  source $STD_ATDM_DIR/$ATDM_CONFIG_KNOWN_SYSTEM_NAME/all_supported_builds.sh
  ATDM_JOB_NAME_KEYS_LIST="$ATDM_CONFIG_ALL_SUPPORTED_BUILDS"
fi

ATDM_JOB_NAME_KEYS_COMMA_LIST=
ATDM_NUM_BULDS=0
for ATDM_JOB_NAME_KEYS in $ATDM_JOB_NAME_KEYS_LIST ; do
  if [[ "$ATDM_JOB_NAME_KEYS_COMMA_LIST" == "" ]] ; then
    ATDM_JOB_NAME_KEYS_COMMA_LIST="$ATDM_JOB_NAME_KEYS"
  else
    ATDM_JOB_NAME_KEYS_COMMA_LIST="$ATDM_JOB_NAME_KEYS_COMMA_LIST,$ATDM_JOB_NAME_KEYS"
  fi
  ATDM_NUM_BULDS=$((ATDM_NUM_BULDS+1))
done

# A.2) Look for --pull and --push arguments and pull them out

ATDM_CHT_FOUND_PULL=0
ATDM_CHT_FOUND_PUSH=0
ARG_ITX=1

for ATDM_CHT_CURENT_ARG in "$@" ; do
  if [[ "$ATDM_CHT_CURENT_ARG" == "--pull" ]] ; then
    echo "Found --pull"
    ATDM_CHT_FOUND_PULL=1
  fi
  if [[ "ATDM_CHT_CURENT_ARG" == "--push" ]] ; then
    echo "Found --push"
    ATDM_CHT_FOUND_PUSH=1
  fi
done

#
# B) Do an initial pull
#

# ToDo: Implement

#
# C) Loop over individual builds and run them
#

echo
echo "Running configure, build, and/or testing for $ATDM_NUM_BULDS builds: $ATDM_JOB_NAME_KEYS_LIST"
echo

ATDM_CHT_BUILD_CASE_IDX=0
for ATDM_JOB_NAME_KEYS in $ATDM_JOB_NAME_KEYS_LIST ; do
  echo
  echo "***"
  echo "*** $ATDM_CHT_BUILD_CASE_IDX) Process build case $ATDM_JOB_NAME_KEYS"
  echo "***"
  echo
  $STD_ATDM_DIR/utils/checkin-test-atdm-single.sh $ATDM_JOB_NAME_KEYS "$@" \
    --default-builds= --allow-no-pull --send-email-to=
  if [[ "$?" == "0" ]] ; then
    echo "$ATDM_JOB_NAME_KEYS: PASSED!"
  else
    echo "$ATDM_JOB_NAME_KEYS: FAILED!"
  fi
  ATDM_CHT_BUILD_CASE_IDX=$((ATDM_CHT_BUILD_CASE_IDX+1))
done

#
# D) Collect the results from all the builds
#

echo
echo "Collect and report final results:"
echo
echo "  ==> See output file checkin-test.final.out" 
echo

$ATDM_TRILINOS_DIR/cmake/tribits/ci_support/checkin-test.py \
  --allow-no-pull --default-builds= --st-extra-builds=$ATDM_JOB_NAME_KEYS_COMMA_LIST \
  &> checkin-test.final.out

ATDM_CHT_RETURN_CODE=$?

# NOTE The return value will be 0 if everything passed!

# Print final status
echo
grep -A 1000 "Commit status email being sent" checkin-test.final.out \
  | grep -B 1000 "Commits for repo" \
  | grep -v "Commit status email being sent" \
  | grep -v "Commits for repo"
echo
grep "REQUESTED ACTIONS" checkin-test.final.out

# ToDo: Add logic to --push if requested!

exit $ATDM_CHT_RETURN_CODE
