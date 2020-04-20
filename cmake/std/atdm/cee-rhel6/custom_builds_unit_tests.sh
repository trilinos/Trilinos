#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/..`

#
# Test compiler parsing
#


testAll() {

  ATDM_CONFIG_SYSTEM_DIR=${ATDM_CONFIG_SCRIPT_DIR}/cee-rhel6

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-clang-9.0.1-openmpi-4.0.2_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-clang-9.0_1-openmpi-4.0.2_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-clang-9.0_1-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-clang-9-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before-clang-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before_gnu-7.2.0-openmpi-4.0.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before_gnu-7.2.0_openmpi-4.0.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before_gnu-7.2.0-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before_gnu-7-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before_gnu-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.2

  ATDM_CONFIG_BUILD_NAME=before_intel-18.0.2-mpich2-3.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-18.0.2_MPICH2-3.2

  ATDM_CONFIG_BUILD_NAME=before_intel-18.0.2_mpich2-3.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-18.0.2_MPICH2-3.2

  ATDM_CONFIG_BUILD_NAME=before_intel-18.0.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-18.0.2_MPICH2-3.2

  ATDM_CONFIG_BUILD_NAME=before_intel-18-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-18.0.2_MPICH2-3.2

  ATDM_CONFIG_BUILD_NAME=before-intel-19.0.3-intelmpi-2018.4_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-19.0.3_INTELMPI-2018.4

  ATDM_CONFIG_BUILD_NAME=before-intel-19.0.3_intelmpi-2018.4_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-19.0.3_INTELMPI-2018.4

  ATDM_CONFIG_BUILD_NAME=before-intel-19.0.3_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-19.0.3_INTELMPI-2018.4

  ATDM_CONFIG_BUILD_NAME=before-intel-19_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-19.0.3_INTELMPI-2018.4

  ATDM_CONFIG_BUILD_NAME=before-intel_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-19.0.3_INTELMPI-2018.4

  # Run the error case
  ATDM_CONFIG_BUILD_NAME=anything-unsupported-compiler
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} DEFAULT

}


#
# Run the unit tests
#

. ${ATDM_CONFIG_SCRIPT_DIR}/test/shunit2/shunit2
