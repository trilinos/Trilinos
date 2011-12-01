# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

INCLUDE(FindListElement)
INCLUDE(MessageWrapper)
INCLUDE(Join)


# Define the valid categories that will be recognized in the CATEGORIES keyword
SET(${PROJECT_NAME}_VALID_CATEGORIES  BASIC  CONTINUOUS NIGHTLY  PERFORMANCE)

# TODO: ABove, We may want only the final project to define these categories
# and not just be general categories for all projects based on ProjectArch.
# Given the logic below, only the categories BASIC and NIGHTLY are specially
# recognized.

# This is a string used in help and error messages
JOIN(${PROJECT_NAME}_VALID_CATEGORIES_STR ", "  FALSE ${${PROJECT_NAME}_VALID_CATEGORIES})

#
# Check for invalid categories called as:
#
#  TRIBITS_GET_INVALID_CATEGORIES(INVLAID_CATEGORIES CATEGORIES_LIST)
#
#  The list of categories to check comes in through the ARGN list.
#
FUNCTION(TRIBITS_GET_INVALID_CATEGORIES  INVALID_CATEGORIES_OUT)
  #MESSAGE("TRIBITS_GET_INVALID_CATEGORIES: ${INVALID_CATEGORIES_OUT} ${ARGN}")
  SET(INVALID_CATEGORIES "")
  FOREACH(CATEGORY_IN ${ARGN})
    #PRINT_VAR(CATEGORY_IN)
    SET(FOUND_CATEGORY FALSE)
    FIND_LIST_ELEMENT(${PROJECT_NAME}_VALID_CATEGORIES ${CATEGORY_IN} FOUND_CATEGORY)
    IF (NOT FOUND_CATEGORY)
      #MESSAGE(STATUS "Not found in list of valid categories!")
      SET(INVALID_CATEGORIES ${INVALID_CATEGORIES} ${CATEGORY_IN})
    ENDIF()
    #PRINT_VAR(INVALID_CATEGORIES)
  ENDFOREACH()
  SET(${INVALID_CATEGORIES_OUT} ${INVALID_CATEGORIES} PARENT_SCOPE)
  #PRINT_VAR(${INVALID_CATEGORIES_OUT})
ENDFUNCTION()


#
# Assert there are no invalid categories called as:
#
#  TRIBITS_ASSERT_VALID_CATEGORIES(CATEGORIES_LIST)
#
#  The list of categories to check comes in through the ARGN list.
#
FUNCTION(TRIBITS_ASSERT_VALID_CATEGORIES)
  #MESSAGE("TRIBITS_ASSERT_VALID_CATEGORIES: ${ARGN}")
  SET(INVALID_CATEGORIES "DUMMYCAT")
  TRIBITS_GET_INVALID_CATEGORIES(INVALID_CATEGORIES ${ARGN})
  #PRINT_VAR(INVALID_CATEGORIES)
  IF (INVALID_CATEGORIES)
    MESSAGE_WRAPPER(SEND_ERROR "Error: The categories '${INVALID_CATEGORIES}' are not"
      " in the list of valid categories '${${PROJECT_NAME}_VALID_CATEGORIES_STR}'!")
  ENDIF()
ENDFUNCTION()
