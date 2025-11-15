#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`

source ${CURRENT_SCRIPTS_DIR}/get-changed-trilinos-packages-helpers.sh

#
# Unit tests
#

test_comma_list_to_list() {

  list=$(comma_list_to_list "a")
  assertEquals "${list}" "a"
  
  list=$(comma_list_to_list "aaa,bbb,ccc")
  assertEquals "${list}" "aaa bbb ccc"

}


test_list_to_comma_list() {

  comma_list=$(list_to_comma_list "a")
  assertEquals "${comma_list}" "a"

  comma_list=$(list_to_comma_list "aaa bbb ccc")
  assertEquals "${comma_list}" "aaa,bbb,ccc"

}


test_list_contains_ele() {

  list_contains_ele "aaa" "aaa bbb ccc"
  ${_ASSERT_EQUALS_} $? 0

  list_contains_ele "aa" "aaa bbb ccc"
  ${_ASSERT_EQUALS_} $? 1

  list_contains_ele "bbb" "aaa bbb ccc"
  ${_ASSERT_EQUALS_} $? 0

  list_contains_ele "ccc" "aaa bbb ccc"
  ${_ASSERT_EQUALS_} $? 0

  list_contains_ele "aa" ""
  ${_ASSERT_EQUALS_} $? 1

  list_contains_ele "aa" ""
  ${_ASSERT_EQUALS_} $? 1

  list_contains_ele "" ""
  ${_ASSERT_EQUALS_} $? 1

  list_contains_ele "aaa" "bbb"
  ${_ASSERT_EQUALS_} $? 1

  list_contains_ele "bbb" "aa"
  ${_ASSERT_EQUALS_} $? 1

}


test_comma_list_contains_ele() {

  comma_list_contains_ele "aaa" "aaa,bbb,ccc"
  ${_ASSERT_EQUALS_} $? 0

  comma_list_contains_ele "aa" "aaa,bbb,ccc"
  ${_ASSERT_EQUALS_} $? 1

  comma_list_contains_ele "bbb" "aaa,bbb,ccc"
  ${_ASSERT_EQUALS_} $? 0

  comma_list_contains_ele "ccc" "aaa,bbb,ccc"
  ${_ASSERT_EQUALS_} $? 0

  comma_list_contains_ele "aa" ""
  ${_ASSERT_EQUALS_} $? 1

  comma_list_contains_ele "aa" ""
  ${_ASSERT_EQUALS_} $? 1

  comma_list_contains_ele "" ""
  ${_ASSERT_EQUALS_} $? 1

  comma_list_contains_ele "aaa" "bbb"
  ${_ASSERT_EQUALS_} $? 1

  comma_list_contains_ele "aa" "bbb"
  ${_ASSERT_EQUALS_} $? 1

}


test_trilinos_filter_packages_to_test() {

  generate_trilinos_package_dependencies_xml_file

  TRILINOS_EXCLUDE_PACKAGES_FROM_PR_TESTING=

  filtered_packages=$(trilinos_filter_packages_to_test "")
  assertEquals "${filtered_packages}" ""

  filtered_packages=$(trilinos_filter_packages_to_test "Teuchos,Tpetra")
  ${_ASSERT_EQUALS_} "${filtered_packages}" "Teuchos,Tpetra"

  filtered_packages=$(trilinos_filter_packages_to_test "Teuchos,Tpetra,PyTrilinos,Panzer")
  ${_ASSERT_EQUALS_} "${filtered_packages}" "Teuchos,Tpetra,PyTrilinos,Panzer"

  filtered_packages=$(trilinos_filter_packages_to_test "Teuchos,Tpetra,PyTrilinos,Panzer")
  ${_ASSERT_EQUALS_} "${filtered_packages}" "Teuchos,Tpetra,PyTrilinos,Panzer"

  TRILINOS_EXCLUDE_PACKAGES_FROM_PR_TESTING=(PyTrilinos)

  filtered_packages=$(trilinos_filter_packages_to_test "")
  assertEquals "${filtered_packages}" ""

  filtered_packages=$(trilinos_filter_packages_to_test "Teuchos,Tpetra")
  ${_ASSERT_EQUALS_} "${filtered_packages}" "Teuchos,Tpetra"

  filtered_packages=$(trilinos_filter_packages_to_test "Teuchos,Tpetra,PyTrilinos,Panzer")
  ${_ASSERT_EQUALS_} "${filtered_packages}" "Teuchos,Tpetra,Panzer"

  filtered_packages=$(trilinos_filter_packages_to_test "Teuchos,Tpetra,PyTrilinos,Panzer")
  ${_ASSERT_EQUALS_} "${filtered_packages}" "Teuchos,Tpetra,Panzer"

}


test_trilinos_filter_packages_to_test() {

  generate_trilinos_package_dependencies_xml_file

  all_toplevel_packages=$(trilinos_get_all_toplevel_packages)
  #echo "all_toplevel_packages='${all_toplevel_packages}'"
  assertContains "${all_toplevel_packages}" "TrilinosFrameworkTests,"
  assertContains "${all_toplevel_packages}" ",TrilinosATDMConfigTests,"
  assertContains "${all_toplevel_packages}" ",Teuchos,"
  assertContains "${all_toplevel_packages}" ",Tpetra,"
  assertContains "${all_toplevel_packages}" ",PyTrilinos2,"
  assertContains "${all_toplevel_packages}" ",Panzer,"

}


#
# Run the unit tests
#

SHUNIT2_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/../test/shunit2`
. ${SHUNIT2_DIR}/shunit2
