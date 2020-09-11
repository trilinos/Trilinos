#!/bin/bash

#
# Get location of directories
#

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`

if [[ "$(uname)" == "Darwin" ]]; then
  ATDM_CONFIG_SCRIPT_DIR="../.."
  ATDM_UTIL_SCRIPT_GET_KNOWN_SYSTEM_INFO="${ATDM_CONFIG_SCRIPT_DIR}/utils/get_known_system_info.sh"
  SHUNIT2_DIR="${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2"
else
  ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/../..`
  ATDM_UTIL_SCRIPT_GET_KNOWN_SYSTEM_INFO=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_known_system_info.sh`
  SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
fi

ATDM_UNIT_TESTS_DIR=${ATDM_CONFIG_SCRIPT_DIR}/test/unit_tests


#
# Unit test helper functions
# 


# Create an intial default test env
function init_system_info_test_env() {
  # Args
  realhostname=$1
  # Body
  SNLSYSTEM=
  SEMS_PLATFORM=
  ATDM_SYSTEM_NAME=
  SNLCLUSTER=
  HOST=
  ATDM_CONFIG_SEMS_GET_PLATFORM=/fake/path/for/unit/testing/
  ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING=${realhostname}
  ATDM_CONFIG_BUILD_NAME=default
  ATDM_CONFIG_DISABLE_WARNINGS=ON
  unset ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG
  unset ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR
  unset ATDM_CONFIG_SYSTEM_NAME
  unset ATDM_CONFIG_REAL_HOSTNAME
  unset ATDM_CONFIG_CDASH_HOSTNAME
  unset ATDM_CONFIG_SYSTEM_NAME
  unset ATDM_CONFIG_SYSTEM_DIR
  unset ATDM_CONFIG_CUSTOM_CONFIG_DIR
}


function assert_system_info_output() {
  # Args
  realhostname=$1
  expected_cdash_hostname=$2
  expected_system_name=$3
  # Body
  ${_ASSERT_EQUALS_} '"${realhostname}"'  '"${ATDM_CONFIG_REAL_HOSTNAME}"'
  ${_ASSERT_EQUALS_} '"${expected_cdash_hostname}"'  '"${ATDM_CONFIG_CDASH_HOSTNAME}"'
  ${_ASSERT_EQUALS_} '"${expected_system_name}"'  '"${ATDM_CONFIG_SYSTEM_NAME}"'
  ${_ASSERT_EQUALS_} '"${ATDM_CONFIG_SCRIPT_DIR}/${expected_system_name}"' \
    '"${ATDM_CONFIG_SYSTEM_DIR}"'
  ${_ASSERT_EQUALS_} '""' '"${ATDM_CONFIG_CUSTOM_SYSTEM_DIR}"'
}


function assert_custom_system_info_output() {
  # Args
  realhostname=$1
  expected_cdash_hostname=$2
  expected_system_name=$3
  expected_system_dir=$4
  # Body
  assertEquals ${expected_system_name}  "${ATDM_CONFIG_SYSTEM_NAME}"
  # NOTE: Above works even if arg is empty to help show this was not set!
  ${_ASSERT_EQUALS_} '"${realhostname}"'  '"${ATDM_CONFIG_REAL_HOSTNAME}"'
  ${_ASSERT_EQUALS_} '"${expected_cdash_hostname}"'  '"${ATDM_CONFIG_CDASH_HOSTNAME}"'
  ${_ASSERT_EQUALS_} '"${expected_system_name}"'  '"${ATDM_CONFIG_SYSTEM_NAME}"'
  ${_ASSERT_EQUALS_} '"${expected_system_dir}"'  '"${ATDM_CONFIG_SYSTEM_DIR}"'
  ${_ASSERT_EQUALS_} '"${expected_system_dir}"'  '"${ATDM_CONFIG_CUSTOM_CONFIG_DIR}"'
}


function known_system_by_hostname_test_helper() {
  # Args
  realhostname=$1
  expected_cdash_hostname=$2
  expected_system_name=$3
  # Body
  init_system_info_test_env "${realhostname}" 
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output "${realhostname}" "${expected_cdash_hostname}" \
    "${expected_system_name}"
}


#
# Test specific some known configurations and basic use cases
#

function test_ride_default() {
  known_system_by_hostname_test_helper ride11 ride ride
}

function test_white_default() {
  known_system_by_hostname_test_helper white22 white ride
}

function test_vortex_default() {
  known_system_by_hostname_test_helper vortex60 vortex ats2
}

function test_ats1_by_hostname() {
  init_system_info_test_env mutrino
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output mutrino mutrino ats1
}

function test_ats1_by_hostname() {
  init_system_info_test_env mtnode52  # Hostname does not matter, just has not to match
  HOST=mutrino
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output mtnode52 mutrino ats1
}

function test_van1_tx2() {
  init_system_info_test_env astra22
  SNLSYSTEM=astra45
  SNLCLUSTER=astra-test
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output astra22 astra-test van1-tx2
}

function test_cts1empire() {
  init_system_info_test_env dummy45  # hostname ignored for matching
  SNLSYSTEM=cts1
  SNLCLUSTER=chama-test
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy45 chama-test cts1empire
}

function test_cts1() {
  init_system_info_test_env dummy45  # hostname ignored for matching
  SNLSYSTEM=cts1
  SNLCLUSTER=chama-test
  ATDM_CONFIG_BUILD_NAME=cts1-default  # Have to put 'cts1' in build name!
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy45 chama-test cts1
}

function test_tlcc2() {
  init_system_info_test_env dummy22  # hostname ignored for matching
  SNLSYSTEM=tlcc2-92
  SNLCLUSTER=dummy-cluster
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy22 dummy-cluster tlcc2
}

function test_sems_rhel6_sems_platform() {
  init_system_info_test_env dummy-rhel6
  SEMS_PLATFORM=rhel6-x86_64
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel6 sems-rhel6 sems-rhel6
}

function test_sems_rhel7_sems_platform() {
  init_system_info_test_env dummy-rhel7
  SEMS_PLATFORM=rhel7-x86_64
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel7 sems-rhel7 sems-rhel7
}

function test_sems_rhel6_sems_get_platform() {
  init_system_info_test_env dummy-rhel6
  ATDM_CONFIG_SEMS_GET_PLATFORM=${ATDM_UNIT_TESTS_DIR}/sems_rhel6_get_platform.sh
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel6 sems-rhel6 sems-rhel6
}

function test_sems_rhel7_sems_get_platform() {
  init_system_info_test_env dummy-rhel7
  ATDM_CONFIG_SEMS_GET_PLATFORM=${ATDM_UNIT_TESTS_DIR}/sems_rhel7_get_platform.sh
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel7 sems-rhel7 sems-rhel7
}

function test_cee_rhel6_linux_rh6() {
  init_system_info_test_env dummy-rhel6
  SNLSYSTEM=cee
  SNLCLUSTER=linux_rh6
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel6 cee-rhel6 cee-rhel6
}

function test_cee_rhel6_linux_rh7() {
  init_system_info_test_env dummy-rhel7
  SNLSYSTEM=cee
  SNLCLUSTER=linux_rh7
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel7 cee-rhel6 cee-rhel6
}


#
# Test cases where multiple known systems are supported and test cases for
# selecting which one
#

function test_sems_rhel7_cee_rhel6_select_sems_rhel7() {
  init_system_info_test_env dummy-rhel7
  SEMS_PLATFORM=rhel7-x86_64
  SNLSYSTEM=cee
  SNLCLUSTER=linux_rh7
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel7 sems-rhel7 sems-rhel7
}

function test_sems_rhel7_cee_rhel6_select_cee_rhel6() {
  init_system_info_test_env dummy-rhel7
  SEMS_PLATFORM=rhel7-x86_64
  SNLSYSTEM=cee
  SNLCLUSTER=linux_rh7
  ATDM_CONFIG_BUILD_NAME=cee-rhel6-default  # Have to put 'cee-rhel6' in name!
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel7 cee-rhel6 cee-rhel6
}

function test_cts1_sems_sems_rhel7_select_cts1emprire() {
  init_system_info_test_env dummy-cts1
  SEMS_PLATFORM=rhel7-x86_64
  SNLSYSTEM=cts1
  SNLCLUSTER=chama-test
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-cts1 chama-test cts1empire
}

function test_cts1_sems_sems_rhel7_select_sems_rhel7() {
  init_system_info_test_env dummy-cts1
  SEMS_PLATFORM=rhel7-x86_64
  SNLSYSTEM=cts1
  SNLCLUSTER=chama-test
  ATDM_CONFIG_BUILD_NAME=sems-rhel7-default   # Must have 'sems-rhel7' in name
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-cts1 sems-rhel7 sems-rhel7
}


#
# Test use cases for custom builds
#

function test_custom_system_arg_with_default_known_system() {
  thishost=$(hostname)
  init_system_info_test_env ${thishost}
  SEMS_PLATFORM=rhel7-x86_64  # Makes sure there is a known default system
  ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG=${ATDM_UNIT_TESTS_DIR}/dummy_custom_system
  ATDM_CONFIG_BUILD_NAME=dummy_custom_system-default # Must be in build name!
  # Check the STDOUT
  output=$(source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh)
  assertEquals \
    "Selecting custom system configuration 'dummy_custom_system'" \
    "${output}"
  # Check the vars set
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_custom_system_info_output ${thishost} ${thishost} dummy_custom_system \
    ${ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG}
}

function test_custom_system_reg_with_default_known_system() {
  thishost=$(hostname)
  init_system_info_test_env ${thishost}
  SEMS_PLATFORM=rhel7-x86_64  # Makes sure there is a known default system
  ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR=${ATDM_UNIT_TESTS_DIR}/dummy_custom_system
  ATDM_CONFIG_BUILD_NAME=dummy_custom_system-default # Must be in build name!
  # Check the STDOUT
  output=$(source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh)
  assertEquals \
    "Selecting custom system configuration 'dummy_custom_system'" \
    "${output}"
  # Check the vars set
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_custom_system_info_output ${thishost} ${thishost} dummy_custom_system \
    ${ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR}
}

# ToDo: Add add test for ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR with no
# default known system!


#
# Test some failure error conditions
#

# ToDo: Implement!


#
# Run the unit tests
#

. ${SHUNIT2_DIR}/shunit2
