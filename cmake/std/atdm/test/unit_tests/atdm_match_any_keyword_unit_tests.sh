#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/../..`

. ${ATDM_CONFIG_SCRIPT_DIR}/utils/define_atdm_match_keyword.sh

#
# Unit tests
#

testMatch_1() {
  
  atdm_match_any_keyword "boo" "foo"
  ${_ASSERT_EQUALS_} $? 1
  
  atdm_match_any_keyword "foo" "foo"
  ${_ASSERT_EQUALS_} $? 0
  
  atdm_match_any_keyword "foo" "food"
  ${_ASSERT_EQUALS_} $? 1
  
  atdm_match_any_keyword "first-foo-last" "food"
  ${_ASSERT_EQUALS_} $? 1
  
  atdm_match_any_keyword "first-food-last" "foo"
  ${_ASSERT_EQUALS_} $? 1
  
  atdm_match_any_keyword "foo" "cow" "food" "foo" "dogs"
  ${_ASSERT_EQUALS_} $? 0
  
  atdm_match_any_keyword "foo" "cow" "food" "dogs"
  ${_ASSERT_EQUALS_} $? 1

}

#
# Run the unit tests
#

. ${ATDM_CONFIG_SCRIPT_DIR}/test/shunit2/shunit2
