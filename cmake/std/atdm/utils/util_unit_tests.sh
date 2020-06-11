#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`

if [[ "$(uname)" == "Darwin" ]]; then
  ATDM_CONFIG_SCRIPT_DIR=".."
  ATDM_UTIL_SCRIPT_ATDM_CONFIG_HELPER_FUNCS="atdm_config_helper_funcs.sh"
  ATDM_UTIL_SCRIPT_GET_KNOWN_SYSTEM_INFO="get_known_system_info.sh"
  SHUNIT2_DIR="${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2"
else
  ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/..`
  ATDM_UTIL_SCRIPT_ATDM_CONFIG_HELPER_FUNCS=`readlink -f ${CURRENT_SCRIPTS_DIR}/atdm_config_helper_funcs.sh`
  ATDM_UTIL_SCRIPT_GET_KNOWN_SYSTEM_INFO=`readlink -f ${CURRENT_SCRIPTS_DIR}/get_known_system_info.sh`
  SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
fi

#
# Test atdm get_known_system_info script
#
test_atdm_get_known_system_info() {
   # TODO: D.2 needs to be tested by setting these values:
  SNLSYSTEM=
  SEMS_PLATFORM=
  ATDM_SYSTEM_NAME=
  SNLCLUSTER=

  # Get the good ATDM_KNOWN_SYSTEM_NAMES_LIST
  ATDM_CONFIG_BUILD_NAME=unit_test
  ATDM_CONFIG_DISABLE_WARNINGS=ON
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

source ${ATDM_UTIL_SCRIPT_ATDM_CONFIG_HELPER_FUNCS}
#
# Test atdm utility functions
#

test_atdm_remove_substrings_from_env_var() {
  for DELIM in ":" "," "-" "_"; do
    # Don't remove anything from path
    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="$EXPECTED_ENV_VAR"
    STRINGS="/test/path1"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    # Remove dir from end of path
    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="$EXPECTED_ENV_VAR$DELIM/test/path1"
    STRINGS="/test/path1"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    # Remove dir at beginning of path
    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="/test/path1$DELIM$EXPECTED_ENV_VAR"
    STRINGS="/test/path1"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    # Remove dir at beginning of path with similar path beside it
    EXPECTED_ENV_VAR="/paths$DELIM/paths2"
    ENV_VAR="/test/path1$DELIM$EXPECTED_ENV_VAR"
    STRINGS="/test/path1"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    # Remove dirs at beginning and end of path
    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="/test/path1$DELIM$EXPECTED_ENV_VAR$DELIM/test/path2"
    STRINGS="/test/path1 /test/path2"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    # Only one dir matches to be removed (other non-matching dir is ignored)
    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="/test/path2$DELIM/test/path1$DELIM$EXPECTED_ENV_VAR"
    STRINGS="/test/path1 /test/path2"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}
  done

  for DELIM in ";" "." "/" " "; do
    RET=$(atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS")
    assertEquals "ERROR: atdm_remove_substrings_from_env_var: \"$DELIM\" is an invalid delimiter." "$RET"
  done
}

#
# Run the unit tests
#
. ${SHUNIT2_DIR}/shunit2
