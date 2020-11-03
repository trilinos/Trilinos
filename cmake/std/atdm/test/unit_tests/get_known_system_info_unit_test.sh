#!/bin/bash

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

#
# Test atdm get_known_system_info script
#
test_atdm_get_known_system_info() {
  # Tweak these values for testing on other platforms
  # TODO: Test branch D.2 in get_known_system_info.sh
  SNLSYSTEM=
  SEMS_PLATFORM=
  ATDM_SYSTEM_NAME=
  SNLCLUSTER=
  HOST=
  ATDM_CONFIG_SEMS_GET_PLATFORM=/fake/path/for/unit/testing/
  ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING=unit_test
  ATDM_CONFIG_BUILD_NAME=$ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING
  ATDM_CONFIG_DISABLE_WARNINGS=ON
  # Populate ATDM_KNOWN_SYSTEM_NAMES_LIST
  source ${ATDM_UTIL_SCRIPT_GET_KNOWN_SYSTEM_INFO}

  # Check that all known system names pass
  for ATDM_CONFIG_BUILD_NAME in ${ATDM_KNOWN_SYSTEM_NAMES_LIST[@]}; do
    ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING=$ATDM_CONFIG_BUILD_NAME
    RET=$(source ${ATDM_UTIL_SCRIPT_GET_KNOWN_SYSTEM_INFO})
    assertEquals "Hostname '$ATDM_CONFIG_BUILD_NAME' matches known ATDM host '$ATDM_CONFIG_BUILD_NAME' and system '$ATDM_CONFIG_BUILD_NAME'" "$RET"
  done

  # Set the bad ATDM_KNOWN_SYSTEM_NAMES
  ATDM_KNOWN_SYSTEM_NAMES_LIST=(
    dne-name1
    dne-name2
    dne-name3
    )

  # Check a few bad systems names for failure
  for ATDM_CONFIG_BUILD_NAME in ${ATDM_KNOWN_SYSTEM_NAMES_LIST[@]}; do
    ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING=$ATDM_CONFIG_BUILD_NAME
    source ${ATDM_UTIL_SCRIPT_GET_KNOWN_SYSTEM_INFO}
    assertEquals "" "$ATDM_SYSTEM_NAME"
    assertEquals "$ATDM_CONFIG_BUILD_NAME" "$realHostname"
  done

  # Ensure that cts1empire is the default on cts1
  ATDM_CONFIG_BUILD_NAME=default
  ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING=$ATDM_CONFIG_BUILD_NAME
  SNLSYSTEM=cts1
  RET=$(source ${ATDM_UTIL_SCRIPT_GET_KNOWN_SYSTEM_INFO})
  assertEquals "Hostname '$ATDM_CONFIG_BUILD_NAME' matches known ATDM host '$ATDM_CONFIG_BUILD_NAME' and system 'cts1empire'" "$RET"

  # Ensure that cts1 is still selected when it's in the build name
  ATDM_CONFIG_BUILD_NAME=cts1-default
  ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING=$ATDM_CONFIG_BUILD_NAME
  SNLSYSTEM=cts1
  RET=$(source ${ATDM_UTIL_SCRIPT_GET_KNOWN_SYSTEM_INFO})
  assertEquals "Hostname '$ATDM_CONFIG_BUILD_NAME' matches known ATDM host '$ATDM_CONFIG_BUILD_NAME' and system 'cts1'" "$RET"
}

#
# Run the unit tests
#
. ${SHUNIT2_DIR}/shunit2
