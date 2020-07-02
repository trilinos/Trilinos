#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/..`

#
# Test compiler parsing
#


testAll() {

  ATDM_CONFIG_SYSTEM_DIR=${ATDM_CONFIG_SCRIPT_DIR}/van1-tx2

  #### TEST ARM-20.0 ####
  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ARM-20.0_OPENMPI-4.0.2 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before-arm-20.0-openmpi-4.0.2_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ARM-20.0_OPENMPI-4.0.2 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before-arm-20.0_openmpi-4.0.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ARM-20.0_OPENMPI-4.0.2 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_arm-20.0-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ARM-20.0_OPENMPI-4.0.2 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_arm-20-afert
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ARM-20.0_OPENMPI-4.0.2 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before-arm-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ARM-20.0_OPENMPI-4.0.2 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=BEFORE-ARM-AFTER
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ARM-20.0_OPENMPI-4.0.2 ${ATDM_CONFIG_COMPILER}

  # Make sure 'arms' does not match 'arm'
  ATDM_CONFIG_BUILD_NAME=anything-arms
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} DEFAULT ${ATDM_CONFIG_COMPILER}

  #### Test ARM-20.1 ####
  ATDM_CONFIG_BUILD_NAME=before-arm-20.1-openmpi-4.0.3_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ARM-20.1_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before-arm-20.1_openmpi-4.0.3-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ARM-20.1_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_arm-20.1-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_}  ARM-20.1_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

}


#
# Run the unit tests
#

SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
. ${SHUNIT2_DIR}/shunit2
