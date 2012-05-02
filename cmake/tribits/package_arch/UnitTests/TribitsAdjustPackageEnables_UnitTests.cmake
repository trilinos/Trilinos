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

MESSAGE("PROJECT_NAME = ${PROJECT_NAME}")
MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = ${${PROJECT_NAME}_TRIBITS_DIR}")

SET( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/package_arch"
  )

INCLUDE(TribitsAdjustPackageEnables)
INCLUDE(UnitTestHelpers)
INCLUDE(GlobalSet)


#####################################################################
#
# Unit tests for code in TribitsAdjustPackageEnables.cmake
#
#####################################################################


FUNCTION(UNITTEST_READ_PACKAGES_LIST_WITH_EXTRA_REPO)

  MESSAGE("\n***")
  MESSAGE("*** Testing the reading of packages list with extra repo")
  MESSAGE("***\n")

  # Debugging
  #SET(TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS ON)
  #SET(${PROJECT_NAME}_VERBOSE_CONFIGURE ON)

  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${PROJECT_NAME} ".")
  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO_NAME} ${EXTRA_REPO_DIR})

  UNITTEST_COMPARE_CONST( ${PROJECT_NAME}_PACKAGES
    "Teuchos;RTOp;Ex2Package1;Ex2Package2")
  UNITTEST_COMPARE_CONST( ${PROJECT_NAME}_PACKAGE_DIRS
    "packages/teuchos;packages/rtop;extraRepoTwoPackages/package1;extraRepoTwoPackages/package2")
  UNITTEST_COMPARE_CONST( ${PROJECT_NAME}_NUM_PACKAGES 4 )

ENDFUNCTION()


FUNCTION(UNITTEST_STANDARD_PROJECT_DEFAULT_EMAIL_ADDRESS_BASE)

  MESSAGE("\n***")
  MESSAGE("*** Testing the case where the TriBITS project has a default email address base and uses standard package regression email list names")
  MESSAGE("***\n")

  SET(${PROJECT_NAME}_VERBOSE_CONFIGURE ON)

  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${PROJECT_NAME} ".")
  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO_NAME} ${EXTRA_REPO_DIR})
  TRIBITS_READ_ALL_PACKAGE_DEPENDENCIES()

  UNITTEST_COMPARE_CONST(Teuchos_REGRESSION_EMAIL_LIST teuchos-regression@repo.site.gov)
  UNITTEST_COMPARE_CONST(RTOp_REGRESSION_EMAIL_LIST thyra-regression@software.sandia.gov)
  UNITTEST_COMPARE_CONST(Ex2Package1_REGRESSION_EMAIL_LIST ex2-package1-override@some.ornl.gov)
  UNITTEST_COMPARE_CONST(Ex2Package2_REGRESSION_EMAIL_LIST ex2package2-regression@project.site.gov)

ENDFUNCTION()


FUNCTION(UNITTEST_SINGLE_REPOSITORY_EMAIL_LIST)

  MESSAGE("\n***")
  MESSAGE("*** Test setting a single regression email address for all the packages in the first repo but defer to hard-coded package email addresses")
  MESSAGE("***\n")

  # Debugging
  #SET(TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS ON)
  #SET(${PROJECT_NAME}_VERBOSE_CONFIGURE ON)

  SET(${PROJECT_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESSS "my-repo@some.url.com")
  SET(${PROJECT_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE OFF) # Will cause to be ignored!

  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${PROJECT_NAME} ".")
  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO_NAME} ${EXTRA_REPO_DIR})
  TRIBITS_READ_ALL_PACKAGE_DEPENDENCIES()

  UNITTEST_COMPARE_CONST(Teuchos_REGRESSION_EMAIL_LIST "my-repo@some.url.com")
  UNITTEST_COMPARE_CONST(RTOp_REGRESSION_EMAIL_LIST thyra-regression@software.sandia.gov)
  UNITTEST_COMPARE_CONST(Ex2Package1_REGRESSION_EMAIL_LIST ex2-package1-override@some.ornl.gov)
  UNITTEST_COMPARE_CONST(Ex2Package2_REGRESSION_EMAIL_LIST ex2package2-regression@project.site.gov)

ENDFUNCTION()


FUNCTION(UNITTEST_SINGLE_REPOSITORY_EMAIL_LIST_OVERRIDE_0)

  MESSAGE("\n***")
  MESSAGE("*** Test setting a single regression email address for all the packages in the first repo with override")
  MESSAGE("***\n")

  # Debugging
  #SET(TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS ON)
  #SET(${PROJECT_NAME}_VERBOSE_CONFIGURE ON)

  SET(${PROJECT_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESSS "my-repo@some.url.com")
  SET(${PROJECT_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST ON)
  SET(${PROJECT_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE OFF) # Will cause to be ignored!

  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${PROJECT_NAME} ".")
  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO_NAME} ${EXTRA_REPO_DIR})
  TRIBITS_READ_ALL_PACKAGE_DEPENDENCIES()

  UNITTEST_COMPARE_CONST(Teuchos_REGRESSION_EMAIL_LIST "my-repo@some.url.com")
  UNITTEST_COMPARE_CONST(RTOp_REGRESSION_EMAIL_LIST "my-repo@some.url.com")
  UNITTEST_COMPARE_CONST(Ex2Package1_REGRESSION_EMAIL_LIST ex2-package1-override@some.ornl.gov)
  UNITTEST_COMPARE_CONST(Ex2Package2_REGRESSION_EMAIL_LIST ex2package2-regression@project.site.gov)

ENDFUNCTION()


FUNCTION(UNITTEST_SINGLE_REPOSITORY_EMAIL_LIST_OVERRIDE_1)

  MESSAGE("\n***")
  MESSAGE("*** Test setting a single regression email address for all the packages in the second repo with override")
  MESSAGE("***\n")

  # Debugging
  #SET(TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS ON)
  #SET(${PROJECT_NAME}_VERBOSE_CONFIGURE ON)

  SET(${EXTRA_REPO_NAME}_REPOSITORY_MASTER_EMAIL_ADDRESSS "extra-repo@some.url.com")
  SET(${EXTRA_REPO_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST ON)

  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${PROJECT_NAME} ".")
  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO_NAME} ${EXTRA_REPO_DIR})
  TRIBITS_READ_ALL_PACKAGE_DEPENDENCIES()

  UNITTEST_COMPARE_CONST(Teuchos_REGRESSION_EMAIL_LIST teuchos-regression@repo.site.gov)
  UNITTEST_COMPARE_CONST(RTOp_REGRESSION_EMAIL_LIST thyra-regression@software.sandia.gov)
  UNITTEST_COMPARE_CONST(Ex2Package1_REGRESSION_EMAIL_LIST extra-repo@some.url.com)
  UNITTEST_COMPARE_CONST(Ex2Package2_REGRESSION_EMAIL_LIST extra-repo@some.url.com)

ENDFUNCTION()


FUNCTION(UNITTEST_SINGLE_PROJECT_EMAIL_LIST)

  MESSAGE("\n***")
  MESSAGE("*** Test setting a single regression email address for all the packages in a TriBITS Project but defer to hard-coded package email addresses")
  MESSAGE("***\n")

  # Debugging
  #SET(${PROJECT_NAME}_VERBOSE_CONFIGURE ON)

  SET(${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESSS "my-project@some.url.com")
  SET(${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESSS_BASE OFF)
  SET(${PROJECT_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE OFF)

  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${PROJECT_NAME} ".")
  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO_NAME} ${EXTRA_REPO_DIR})
  TRIBITS_READ_ALL_PACKAGE_DEPENDENCIES()

  UNITTEST_COMPARE_CONST(Teuchos_REGRESSION_EMAIL_LIST "my-project@some.url.com")
  UNITTEST_COMPARE_CONST(RTOp_REGRESSION_EMAIL_LIST thyra-regression@software.sandia.gov)
  UNITTEST_COMPARE_CONST(Ex2Package1_REGRESSION_EMAIL_LIST ex2-package1-override@some.ornl.gov)
  UNITTEST_COMPARE_CONST(Ex2Package2_REGRESSION_EMAIL_LIST my-project@some.url.com)

ENDFUNCTION()


FUNCTION(UNITTEST_SINGLE_PROJECT_EMAIL_LIST_OVERRIDE)

  MESSAGE("\n***")
  MESSAGE("*** Test setting a single regression email address for all the packages in a TriBITS Project and overriding hard-coded package email addresses")
  MESSAGE("***\n")


  # Debugging
  #SET(${PROJECT_NAME}_VERBOSE_CONFIGURE ON)

  SET(${PROJECT_NAME}_PROJECT_MASTER_EMAIL_ADDRESSS "my-project@some.url.com")
  SET(${PROJECT_NAME}_PROJECT_EMAIL_URL_ADDRESSS_BASE OFF)
  SET(${PROJECT_NAME}_REPOSITORY_EMAIL_URL_ADDRESSS_BASE OFF)
  SET(${PROJECT_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST ON)
  SET(${EXTRA_REPO_NAME}_REPOSITORY_OVERRIDE_PACKAGE_EMAIL_LIST ON)

  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${PROJECT_NAME} ".")
  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO_NAME} ${EXTRA_REPO_DIR})
  TRIBITS_READ_ALL_PACKAGE_DEPENDENCIES()

  UNITTEST_COMPARE_CONST(Teuchos_REGRESSION_EMAIL_LIST "my-project@some.url.com")
  UNITTEST_COMPARE_CONST(RTOp_REGRESSION_EMAIL_LIST my-project@some.url.com)
  UNITTEST_COMPARE_CONST(Ex2Package1_REGRESSION_EMAIL_LIST my-project@some.url.com)
  UNITTEST_COMPARE_CONST(Ex2Package2_REGRESSION_EMAIL_LIST my-project@some.url.com)

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

SET(PROJECT_SOURCE_DIR "${${PROJECT_NAME}_TRIBITS_DIR}/package_arch/UnitTests/MockTrilinos")
PRINT_VAR(PROJECT_SOURCE_DIR)
SET(REPOSITORY_DIR ".")
PRINT_VAR(REPOSITORY_DIR)

# Set the mock project name last to override the outer project
SET(PROJECT_NAME "Trilinos")

SET( Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  Teuchos             packages/teuchos                PS
  RTOp                packages/rtop                   PS
  )

SET(EXTRA_REPO_NAME extraRepoTwoPackages)
SET(EXTRA_REPO_DIR extraRepoTwoPackages)

INCLUDE(${PROJECT_SOURCE_DIR}/${EXTRA_REPO_NAME}/PackagesList.cmake)

SET(${PROJECT_NAME}_ALL_REPOSITORIES "." "${EXTRA_REPO_NAME}")

SET( ${PROJECT_NAME}_ASSERT_MISSING_PACKAGES ON )


#
# Run the unit tests
#

UNITTEST_READ_PACKAGES_LIST_WITH_EXTRA_REPO()
UNITTEST_STANDARD_PROJECT_DEFAULT_EMAIL_ADDRESS_BASE()
UNITTEST_SINGLE_REPOSITORY_EMAIL_LIST()
UNITTEST_SINGLE_REPOSITORY_EMAIL_LIST_OVERRIDE_0()
UNITTEST_SINGLE_REPOSITORY_EMAIL_LIST_OVERRIDE_1()
UNITTEST_SINGLE_PROJECT_EMAIL_LIST()
UNITTEST_SINGLE_PROJECT_EMAIL_LIST_OVERRIDE()
#UNITTEST_SINGLE_PROJECT_REPOSITORY_EMAIL_LISTS()
#UNITTEST_SINGLE_PROJECT_REPOSITORY_EMAIL_LISTS_OVERRIDE()

# Pass in the number of expected tests that must pass!
UNITTEST_FINAL_RESULT(27)
