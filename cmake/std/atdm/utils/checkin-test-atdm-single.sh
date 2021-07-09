#!/bin/bash

# NOTE: This script should only be run by checkin-test-atdm.sh.  It does not
# have enough info to run on its own.

ATDM_BUILD_NAME_KEYS=$1 ; shift

source $STD_ATDM_DIR/load-env.sh $ATDM_BUILD_NAME_KEYS

if [ "${ATDM_CONFIG_USE_NINJA}" == "ON" ] ; then
  echo "-GNinja" > $ATDM_BUILD_NAME_KEYS.config
  CHECKIN_TEST_USE_NINJA_ARG=--use-ninja
else
  echo > $ATDM_BUILD_NAME_KEYS.config
  CHECKIN_TEST_USE_NINJA_ARG=
fi

echo "-DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake
-DTrilinos_TRACE_ADD_TEST=ON" \
>> $ATDM_BUILD_NAME_KEYS.config

echo
echo "Running: checkin-test.py --st-extra-builds=$ATDM_BUILD_NAME_KEYS ..."
echo
echo "  ==> See output file checkin-test.$ATDM_BUILD_NAME_KEYS.out"
echo

if [ "$ATDM_CONFIG_BUILD_COUNT" -gt "0" ] ; then
  make_options="-j $ATDM_CONFIG_BUILD_COUNT"
else
  make_options=
fi

set -x

${ATDM_TRIBITS_DIR}/ci_support/checkin-test.py \
  --src-dir=$ATDM_TRILINOS_DIR \
  --make-options="${make_options}" \
  --ctest-options="-j $ATDM_CONFIG_CTEST_PARALLEL_LEVEL" \
  --st-extra-builds=$ATDM_BUILD_NAME_KEYS "$@" \
  $CHECKIN_TEST_USE_NINJA_ARG \
  --log-file=checkin-test.$ATDM_BUILD_NAME_KEYS.out \
  &> /dev/null
ATDM_CHT_SINGLE_RETURN_CODE=$?

set +x

echo "source $STD_ATDM_DIR/load-env.sh $ATDM_BUILD_NAME_KEYS" \
 > $ATDM_BUILD_NAME_KEYS/load-env.sh

exit $ATDM_CHT_SINGLE_RETURN_CODE
