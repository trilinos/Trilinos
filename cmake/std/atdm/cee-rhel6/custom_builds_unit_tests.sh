#!/bin/bash

if [[ "$(uname)" == "Darwin" ]]; then
  ATDM_CONFIG_SCRIPT_DIR=..
  source "${ATDM_CONFIG_SCRIPT_DIR}/utils/define_atdm_match_keyword.sh"
  SHUNIT2_DIR="${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2"
else
  CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
  ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/..`
  SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
fi

ATDM_CONFIG_SYSTEM_DIR=${ATDM_CONFIG_SCRIPT_DIR}/cee-rhel6

#
# Test compiler parsing
#

testCuda() {
  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243-gcc-7.2.0_openmpi-4.0.3-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GCC-7.2.0_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243_gcc-7.2.0_openmpi-4.0.3-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GCC-7.2.0_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10.1.243-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GCC-7.2.0_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-10-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GCC-7.2.0_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}

  ATDM_CONFIG_BUILD_NAME=before_cuda-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} CUDA-10.1.243_GCC-7.2.0_OPENMPI-4.0.3 ${ATDM_CONFIG_COMPILER}
  
  # Check that 'somethingcuda' does not match 'cuda'! (Shows true keyword
  # matching is working)
  ATDM_CONFIG_BUILD_NAME=somethingcuda
  STDOUT=$(. ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh 2>&1)
  GREP_LINE=$(echo "${STDOUT}" | grep "ERROR: A supported compiler was not selected")
  GREP_LINE_EXPECTED="*** ERROR: A supported compiler was not selected for 'cee-rhel6' env in buildname 'somethingcuda'"
  #echo "GREP_LINE = [${GREP_LINE}]"
  echo "GREP_LINE_EXPECTED = [${GREP_LINE_EXPECTED}]"
  ${_ASSERT_EQUALS_} '"${GREP_LINE}"' '"${GREP_LINE_EXPECTED}"'

  # Check that missing compiler prints right error messagematch 'gnu'!
  ATDM_CONFIG_BUILD_NAME=help
  STDOUT=$(. ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh 2>&1)
  GREP_LINE=$(echo "${STDOUT}" | grep "ERROR: A supported compiler was not selected")
  GREP_LINE_EXPECTED="*** ERROR: A supported compiler was not selected for 'cee-rhel6' env in buildname 'help'"
  #echo "GREP_LINE = [${GREP_LINE}]"
  echo "GREP_LINE_EXPECTED = [${GREP_LINE_EXPECTED}]"
  ${_ASSERT_EQUALS_} '"${GREP_LINE}"' '"${GREP_LINE_EXPECTED}"'
}

testAll() {
  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=default-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before-clang-9.0.1-openmpi-4.0.3_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before-clang-9.0_1-openmpi-4.0.3_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before-clang-9.0_1-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before-clang-9-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before-clang-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-9.0.1_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before_gnu-7.2.0-openmpi-4.0.3-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before_gnu-7.2.0_openmpi-4.0.3-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before_gnu-7.2.0-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before_gnu-7-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.3

  ATDM_CONFIG_BUILD_NAME=before_gnu-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} GNU-7.2.0_OPENMPI-4.0.3

  # Check that 'gnus' does not match 'gnu'! (Shows true keyword matching is
  # working)
  ATDM_CONFIG_BUILD_NAME=before_gnus-after
  STDOUT=$(. ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh 2>&1)
  GREP_LINE=$(echo "${STDOUT}" | grep "ERROR: A supported compiler was not selected")
  GREP_LINE_EXPECTED="*** ERROR: A supported compiler was not selected for 'cee-rhel6' env in buildname 'before_gnus-after'"
  #echo "GREP_LINE = [${GREP_LINE}]"
  echo "GREP_LINE_EXPECTED = [${GREP_LINE_EXPECTED}]"
  ${_ASSERT_EQUALS_} '"${GREP_LINE}"' '"${GREP_LINE_EXPECTED}"'

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

  # Check that 'somethingintel' does not match 'intl'! (Shows true keyword
  # matching is working)
  ATDM_CONFIG_BUILD_NAME=somethingintel
  STDOUT=$(. ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh 2>&1)
  GREP_LINE=$(echo "${STDOUT}" | grep "ERROR: A supported compiler was not selected")
  GREP_LINE_EXPECTED="*** ERROR: A supported compiler was not selected for 'cee-rhel6' env in buildname 'somethingintel'"
  #echo "GREP_LINE = [${GREP_LINE}]"
  echo "GREP_LINE_EXPECTED = [${GREP_LINE_EXPECTED}]"
  ${_ASSERT_EQUALS_} '"${GREP_LINE}"' '"${GREP_LINE_EXPECTED}"'

  # Check that missing compiler prints right error messagematch 'gnu'!
  ATDM_CONFIG_BUILD_NAME=help
  STDOUT=$(. ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh 2>&1)
  GREP_LINE=$(echo "${STDOUT}" | grep "ERROR: A supported compiler was not selected")
  GREP_LINE_EXPECTED="*** ERROR: A supported compiler was not selected for 'cee-rhel6' env in buildname 'help'"
  #echo "GREP_LINE = [${GREP_LINE}]"
  echo "GREP_LINE_EXPECTED = [${GREP_LINE_EXPECTED}]"
  ${_ASSERT_EQUALS_} '"${GREP_LINE}"' '"${GREP_LINE_EXPECTED}"'

}


#
# Run the unit tests
#

. ${SHUNIT2_DIR}/shunit2
