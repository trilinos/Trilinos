#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/..`

#
# Test compiler parsing
#


testAll() {

  ATDM_CONFIG_SYSTEM_DIR=${ATDM_CONFIG_SCRIPT_DIR}/van1-tx2

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-arm-20.0-openmpi-4.0.2_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-arm-20.0_openmpi-4.0.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before_arm-20.0-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before_arm-20-afert
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-arm-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=BEFORE-ARM-AFTER
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  # Make sure 'arms' does not match 'arm'
  ATDM_CONFIG_BUILD_NAME=anything-arms
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} DEFAULT

}


#
# Run the unit tests
#

. ${ATDM_CONFIG_SCRIPT_DIR}/test/shunit2/shunit2
