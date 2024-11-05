# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include("${CMAKE_CURRENT_LIST_DIR}/../utils/FindListElement.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../utils/MessageWrapper.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../utils/Join.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../utils/TribitsDeprecatedHelpers.cmake")


# Define the valid categories that will be recognized in the CATEGORIES keyword
set(${PROJECT_NAME}_VALID_CATEGORIES  BASIC  CONTINUOUS  NIGHTLY  HEAVY  WEEKLY  PERFORMANCE)

# TODO: ABove, We may want only the final project to define these categories
# and not just be general categories for all projects based on ProjectArch.
# Given the logic below, only the categories BASIC and NIGHTLY are specially
# recognized.

# This is a string used in help and error messages
join(${PROJECT_NAME}_VALID_CATEGORIES_STR ", "  FALSE ${${PROJECT_NAME}_VALID_CATEGORIES})

#
# Check for invalid categories called as:
#
#  tribits_get_invalid_categories(INVLAID_CATEGORIES CATEGORIES_LIST)
#
#  The list of categories to check comes in through the ARGN list.
#
function(tribits_get_invalid_categories  INVALID_CATEGORIES_OUT)
  #message("TRIBITS_GET_INVALID_CATEGORIES: ${INVALID_CATEGORIES_OUT} ${ARGN}")
  set(INVALID_CATEGORIES "")
  foreach(CATEGORY_IN ${ARGN})
    #print_var(CATEGORY_IN)
    set(FOUND_CATEGORY FALSE)
    find_list_element(${PROJECT_NAME}_VALID_CATEGORIES ${CATEGORY_IN} FOUND_CATEGORY)
    if (NOT FOUND_CATEGORY)
      #message(STATUS "Not found in list of valid categories!")
      set(INVALID_CATEGORIES ${INVALID_CATEGORIES} ${CATEGORY_IN})
    endif()
    #print_var(INVALID_CATEGORIES)
  endforeach()
  set(${INVALID_CATEGORIES_OUT} ${INVALID_CATEGORIES} PARENT_SCOPE)
  #print_var(${INVALID_CATEGORIES_OUT})
endfunction()


#
# Assert there are no invalid categories called as:
#
#  tribits_assert_valid_categories(CATEGORIES_LIST)
#
#  The list of categories to check comes in through the ARGN list.
#
function(tribits_filter_and_assert_categories  CATEGORIES_VAR_INOUT)
  #message("TRIBITS_ASSERT_VALID_CATEGORIES: ${ARGN}")
  set(INVALID_CATEGORIES "DUMMYCAT")
  tribits_get_invalid_categories(INVALID_CATEGORIES ${${CATEGORIES_VAR_INOUT}})
  #print_var(INVALID_CATEGORIES)
  if (INVALID_CATEGORIES)
    message_wrapper(SEND_ERROR "Error: The categories '${INVALID_CATEGORIES}' are not"
      " in the list of valid categories '${${PROJECT_NAME}_VALID_CATEGORIES_STR}'!")
  endif()
  set(CATEGORIES_OUT)
  foreach(CATEGORY  ${${CATEGORIES_VAR_INOUT}})
    if (CATEGORY STREQUAL "WEEKLY")
      tribits_deprecated("The test category 'WEEKLY' is deprecated"
        " and is replaced with 'HEAVY'.  Please change to use 'HEAVY' instead.")
      list(APPEND  CATEGORIES_OUT  "HEAVY")
    else()
      list(APPEND  CATEGORIES_OUT  ${CATEGORY})
    endif()
  endforeach()
  set(${CATEGORIES_VAR_INOUT} ${CATEGORIES_OUT} PARENT_SCOPE)
endfunction()
