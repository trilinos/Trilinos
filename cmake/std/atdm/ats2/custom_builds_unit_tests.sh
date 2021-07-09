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

testXl() {
  ATDM_CONFIG_BUILD_NAME=default-xl
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_xl-2020.03.18-spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_xl-2020.03.18_spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_xl-2020.03.18_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_xl-2020_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_xl_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  # This should not match anything and should be an error!
  ATDM_CONFIG_BUILD_NAME=anything-xls-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} DEFAULT ${ATDM_CONFIG_COMPILER}
}

testXlCuda() {
  ATDM_CONFIG_BUILD_NAME=default-cuda-xl
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243-xl-2020.03.18-spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243-xl-2020.03.18_spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243-xl-2020.03.18_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243-xl-2020_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-xl_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243_xl-2020.03.18-spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243_xl-2020.03.18_spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243_xl-2020.03.18_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243_xl-2020_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda_xl_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  # This should not match anything and should be an error!
  ATDM_CONFIG_BUILD_NAME=anything-cudas-xls-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} DEFAULT ${ATDM_CONFIG_COMPILER}
}

#
# Test compiler parsing
#

testGnu() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243-gnu-7.3.1-spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243_gnu-7.3.1_spmpi-rolling_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243-gnu-7.3.1_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243_gnu-7.3.1_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243-gnu-7_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243_gnu-7_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING ${ATDM_CONFIG_COMPILER}

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
