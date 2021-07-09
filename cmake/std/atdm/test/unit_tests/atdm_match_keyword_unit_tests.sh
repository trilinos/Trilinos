#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/../..`

. ${ATDM_CONFIG_SCRIPT_DIR}/utils/define_atdm_match_keyword.sh

#
# Unit tests
#

testBasicMatches() {
  
  atdm_match_keyword "boo" "foo"
  ${_ASSERT_EQUALS_} $? 1
  
  atdm_match_keyword "foo" "foo"
  ${_ASSERT_EQUALS_} $? 0
  
  atdm_match_keyword "foo" "FOO"
  ${_ASSERT_EQUALS_} $? 0
  
  atdm_match_keyword "FOO" "FOO"
  ${_ASSERT_EQUALS_} $? 0
  
  atdm_match_keyword "defaultfoo" "foo"
  ${_ASSERT_EQUALS_} $? 1
  
  atdm_match_keyword "default-foo" "foo"
  ${_ASSERT_EQUALS_} $? 0
  
  atdm_match_keyword "default_foo" "foo"
  ${_ASSERT_EQUALS_} $? 0
  
  atdm_match_keyword "foodefault" "foo"
  ${_ASSERT_EQUALS_} $? 1

  atdm_match_keyword "foo-default" "foo"
  ${_ASSERT_EQUALS_} $? 0

  atdm_match_keyword "foo_default" "foo"
  ${_ASSERT_EQUALS_} $? 0

  atdm_match_keyword "first-foo-last" "foo"
  ${_ASSERT_EQUALS_} $? 0

  atdm_match_keyword "first_foo-last" "foo"
  ${_ASSERT_EQUALS_} $? 0

  atdm_match_keyword "first-foo_last" "foo"
  ${_ASSERT_EQUALS_} $? 0

  atdm_match_keyword "first_foo_last" "foo"
  ${_ASSERT_EQUALS_} $? 0

  atdm_match_keyword "first_foo_last" "FOO"
  ${_ASSERT_EQUALS_} $? 0

  atdm_match_keyword "FIRST_FOO_LAST" "foo"
  ${_ASSERT_EQUALS_} $? 0

}

#
# Run the unit tests
#

SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
. ${SHUNIT2_DIR}/shunit2
