# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

include(FindListElement)
include(MessageWrapper)
include(Join)
include(TribitsDeprecatedHelpers)


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
