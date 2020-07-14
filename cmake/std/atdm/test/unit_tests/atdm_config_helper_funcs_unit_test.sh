#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`

if [[ "$(uname)" == "Darwin" ]]; then
  ATDM_CONFIG_SCRIPT_DIR="../.."
  ATDM_UTIL_SCRIPT_ATDM_CONFIG_HELPER_FUNCS="${ATDM_CONFIG_SCRIPT_DIR}/utils/atdm_config_helper_funcs.sh"
  SHUNIT2_DIR="${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2"
else
  ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/../..`
  ATDM_UTIL_SCRIPT_ATDM_CONFIG_HELPER_FUNCS=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/utils/atdm_config_helper_funcs.sh`
  SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
fi

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
