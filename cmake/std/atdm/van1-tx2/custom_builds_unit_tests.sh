#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/..`

#
# Test compiler parsing
#


testAll() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/van1-tx2/custom_builds.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-arm-20.0-openmpi-4.0.2_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/van1-tx2/custom_builds.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-arm-20.0_openmpi-4.0.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/van1-tx2/custom_builds.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=arm-20.0
  . ${ATDM_CONFIG_SCRIPT_DIR}/van1-tx2/custom_builds.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=arm-20
  . ${ATDM_CONFIG_SCRIPT_DIR}/van1-tx2/custom_builds.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=arm
  . ${ATDM_CONFIG_SCRIPT_DIR}/van1-tx2/custom_builds.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} ARM-20.0_OPENMPI-4.0.2

}


#
# Run the unit tests
#

. ${ATDM_CONFIG_SCRIPT_DIR}/test/shunit2/shunit2
