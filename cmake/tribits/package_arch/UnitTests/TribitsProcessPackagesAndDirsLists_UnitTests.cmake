# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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

MESSAGE("PROJECT_NAME = ${PROJECT_NAME}")
MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = ${${PROJECT_NAME}_TRIBITS_DIR}")

SET( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/package_arch"
  )

INCLUDE(TribitsProcessPackagesAndDirsLists)
INCLUDE(UnitTestHelpers)
INCLUDE(GlobalSet)


#####################################################################
#
# Unit tests for code in TribitsProcessPackagesAndDirsLists.cmake
#
#####################################################################


FUNCTION(UNITTEST_BASIC_PACKAGE_LIST_READ)

  MESSAGE("\n***")
  MESSAGE("*** Testing the basic reading of packages list")
  MESSAGE("***\n")

  SET( ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
    Package0     packages/package0  PS
    Package1     packages/package1  SS
    Package2     packages/package2  EX
    )

  SET(PACKAGE_ABS_DIR "DummyBase")
  SET(${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK TRUE)

  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${PROJECT_NAME} ".")

  UNITTEST_COMPARE_CONST( ${PROJECT_NAME}_PACKAGES
    "Package0;Package1;Package2")
  UNITTEST_COMPARE_CONST( ${PROJECT_NAME}_PACKAGE_DIRS
    "packages/package0;packages/package1;packages/package2")
  UNITTEST_COMPARE_CONST( ${PROJECT_NAME}_NUM_PACKAGES 3 )
  UNITTEST_COMPARE_CONST( ${PROJECT_NAME}_LAST_PACKAGE_IDX 2 )
  UNITTEST_COMPARE_CONST( ${PROJECT_NAME}_REVERSE_PACKAGES
    "Package2;Package1;Package0")

ENDFUNCTION()



#####################################################################
#
# Execute the unit tests
#
#####################################################################

# Assume that all unit tests will pass by default
GLOBAL_SET(UNITTEST_OVERALL_PASS TRUE)
GLOBAL_SET(UNITTEST_OVERALL_NUMPASSED 0)
GLOBAL_SET(UNITTEST_OVERALL_NUMRUN 0)

# Set common/base options
SET(PROJECT_NAME "DummyProject")

UNITTEST_BASIC_PACKAGE_LIST_READ()

# Pass in the number of expected tests that must pass!
UNITTEST_FINAL_RESULT(5)
