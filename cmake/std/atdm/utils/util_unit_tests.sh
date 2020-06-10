#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#ATDM_CONFIG_SCRIPT_DIR=".."
ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/..`
#ATDM_UTIL_SCRIPT="atdm_config_helper_funcs.sh"
ATDM_UTIL_SCRIPT=`readlink -f ${CURRENT_SCRIPTS_DIR}/atdm_config_helper_funcs.sh`

source ${ATDM_UTIL_SCRIPT}

#
# Test atdm utility functions
#

test_atdm_remove_substrings_from_env_var() {
  for DELIM in ":" "," "-" "_"; do
    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="$EXPECTED_ENV_VAR"
    STRINGS="/test/path1"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="$EXPECTED_ENV_VAR$DELIM/test/path1"
    STRINGS="/test/path1"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="/test/path1$DELIM$EXPECTED_ENV_VAR"
    STRINGS="/test/path1"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    EXPECTED_ENV_VAR="/paths$DELIM/paths2"
    ENV_VAR="/test/path1$DELIM$EXPECTED_ENV_VAR"
    STRINGS="/test/path1"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="/test/path1$DELIM$EXPECTED_ENV_VAR$DELIM/test/path2"
    STRINGS="/test/path1 /test/path2"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}

    EXPECTED_ENV_VAR="/paths"
    ENV_VAR="/test/path2$DELIM/test/path1$DELIM$EXPECTED_ENV_VAR"
    STRINGS="/test/path1 /test/path2"
    atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS"
    ${_ASSERT_EQUALS_} ${EXPECTED_ENV_VAR} ${ENV_VAR}
  done

  for DELIM in ";" "." "/" " "; do
    RET=$(atdm_remove_substrings_from_env_var ENV_VAR "$DELIM" "$STRINGS")
    assertEquals "$RET" "ERROR: atdm_remove_substrings_from_env_var: \"$DELIM\" is an invalid delimiter."
  done
}

testAll() {
  test_atdm_remove_substrings_from_env_var
}


#
# Run the unit tests
#
SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
#SHUNIT2_DIR="${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2"
. ${SHUNIT2_DIR}/shunit2
