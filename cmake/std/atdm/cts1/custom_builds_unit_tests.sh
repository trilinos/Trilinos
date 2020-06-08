#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/..`

#
# Test compiler parsing
#

testAll() {

  ATDM_CONFIG_SYSTEM_DIR=${ATDM_CONFIG_SCRIPT_DIR}/cts1

  ATDM_CONFIG_BUILD_NAME=before-intel-18.0.2-openmpi-4.0.3_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-18.0.2_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before-intel-18.0.2_openmpi-4.0.3-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-18.0.2_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_intel-18.0.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-18.0.2_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-18.0.2_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before-intel-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-18.0.2_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=BEFORE-INTEL-AFTER
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-18.0.2_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before-intel-19.0.4-openmpi-4.0.3_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-19.0.4_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before-intel-19.0.4_openmpi-4.0.3-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-19.0.4_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_intel-19.0.4-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-19.0.4_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_intel-19-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} INTEL-19.0.4_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=anything-intell
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} DEFAULT ${ATDM_CONFIG_COMPILER}

}


#
# Run the unit tests
#

SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
. ${SHUNIT2_DIR}/shunit2
