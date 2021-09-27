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

function test_sems_rhel7_sems_platform() {
  init_system_info_test_env dummy-rhel7
  SEMS_PLATFORM=rhel7-x86_64
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel7 sems-rhel7 sems-rhel7
}

function test_sems_rhel7_sems_get_platform() {
  init_system_info_test_env dummy-rhel7
  ATDM_CONFIG_SEMS_GET_PLATFORM=${ATDM_UNIT_TESTS_DIR}/sems_rhel7_get_platform.sh
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel7 sems-rhel7 sems-rhel7
}

function test_cee_rhel7_linux_rh7() {
  init_system_info_test_env dummy-rhel7
  SNLSYSTEM=cee
  SNLCLUSTER=linux_rh7
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel7 cee-rhel7 cee-rhel7
}


#
# Test cases where multiple known systems are supported and test cases for
# selecting which one
#

function test_sems_rhel7_cee_rhel7_select_sems_rhel7() {
  init_system_info_test_env dummy-rhel7
  SEMS_PLATFORM=rhel7-x86_64
  SNLSYSTEM=cee
  SNLCLUSTER=linux_rh7
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-rhel7 sems-rhel7 sems-rhel7
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
  thishost=dummy-machine
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
  thishost=dummy-machine
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

function test_custom_system_reg_no_default_known_system() {
  thishost=dummy-machine
  init_system_info_test_env ${thishost}
  ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR=${ATDM_UNIT_TESTS_DIR}/dummy_custom_system
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


#
# Test some failure error conditions
#

function test_no_supported_system() {
  init_system_info_test_env dummy-machine
  ATDM_CONFIG_BUILD_NAME=default
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_REAL_HOSTNAME}"'
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_CDASH_HOSTNAME}"'
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_SYSTEM_NAME}"'
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_SYSTEM_DIR}"'
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_CUSTOM_CONFIG_DIR}"'
  # NOTE: The error message for no supported system is generated by
  # the load-env.sh script, not the get_system_info.sh script.
}

function test_no_supported_system_try_ats2() {
  init_system_info_test_env dummy-machine
  ATDM_CONFIG_BUILD_NAME=ats2-default
  # Capture the STDOUT and check out
  cmnd_stdout=$(source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh)
  cmnd_stdout_expected="Hostname 'dummy-machine' matches known ATDM host 'dummy-machine' and system 'ats2'"
  assertEquals \
    "${cmnd_stdout_expected}" \
    "${cmnd_stdout}"
  # Run the command and check the vars
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  assert_system_info_output dummy-machine "" ats2
}
# NOTE: The above use case makes no sense.  This should be a failure
# to try to select a system thta is not supported.  We need to fix this so that
# this generates an error.  This current test just pins down this (bad) behavior.

function test_sems_rhel7_cee_rhel7_try_ats2_fail() {
  init_system_info_test_env dummy-rhel7
  SEMS_PLATFORM=rhel7-x86_64
  SNLSYSTEM=cee
  SNLCLUSTER=linux_rh7
  ATDM_CONFIG_BUILD_NAME=ats2-default
  # Capture the STDOUT and check out
  cmnd_stdout=$(source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh)
  cmnd_stdout_expected="
***
*** Error, the system name 'ats2' given in the build name:
***
***   ats2-default
***
*** does not match the current host:
***
***   dummy-rhel7
***
**** and does not match one of the supported system types on this machine which includes:
***
***   (sems-rhel7 cee-rhel7)
***
*** To address this, either remove 'ats2' from the build name
*** or change it to one of the above suppported system types.
***"
  assertEquals \
    "${cmnd_stdout_expected}" \
    "${cmnd_stdout}"
  # Run the command and check the vars
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_REAL_HOSTNAME}"'
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_CDASH_HOSTNAME}"'
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_SYSTEM_NAME}"'
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_SYSTEM_DIR}"'
  ${_ASSERT_EQUALS_} '""'  '"${ATDM_CONFIG_CUSTOM_CONFIG_DIR}"'
}


#
# Run the unit tests
#

. ${SHUNIT2_DIR}/shunit2
