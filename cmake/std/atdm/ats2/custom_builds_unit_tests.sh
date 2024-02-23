#!/bin/bash

if [[ "$(uname)" == "Darwin" ]]; then
  ATDM_CONFIG_SCRIPT_DIR=".."
  source "${ATDM_CONFIG_SCRIPT_DIR}/utils/define_atdm_match_keyword.sh"
  ATDM_CONFIG_SYSTEM_DIR=${ATDM_CONFIG_SCRIPT_DIR}/ats2
  SHUNIT2_DIR="${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2"
else
  CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
  ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/..`
  ATDM_CONFIG_SYSTEM_DIR=${ATDM_CONFIG_SCRIPT_DIR}/ats2
  SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
fi

#
# Test compiler parsing
#

testGnu() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-11.2.152-gcc-8.3.1-spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-11.2.152_gcc-8.3.1_spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-11.2.152-gcc-8.3.1_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-11.2.152_gcc-8.3.1_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-11.2.152-gcc-8_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-11.2.152_gcc-8_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-11.2.152_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-11_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_gnu-7.3.1-spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_gnu-7.3.1_spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_gnu-7.3.1_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_gnu-7_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_gnu_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  # This should not match anything and should be an error!
  ATDM_CONFIG_BUILD_NAME=anything-cudas-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} DEFAULT ${ATDM_CONFIG_COMPILER}

}


#
# Run the unit tests
#

. ${SHUNIT2_DIR}/shunit2
