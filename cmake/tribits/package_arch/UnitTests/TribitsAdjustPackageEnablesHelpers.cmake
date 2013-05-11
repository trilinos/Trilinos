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

INCLUDE(TribitsAdjustPackageEnables)
INCLUDE(TribitsProcessTplsLists)
INCLUDE(UnitTestHelpers)
INCLUDE(GlobalSet)


#####################################################################
#
# Helper macros for unit tests
#
#####################################################################


MACRO(UNITTEST_HELPER_READ_AND_PROESS_PACKAGES)

  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${PROJECT_NAME} ".")
  TRIBITS_PROCESS_TPLS_LISTS(${PROJECT_NAME} ".")
  TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO_NAME} ${EXTRA_REPO_DIR})
  TRIBITS_PROCESS_TPLS_LISTS(${EXTRA_REPO_NAME} ${EXTRA_REPO_DIR})
  TRIBITS_READ_ALL_PACKAGE_DEPENDENCIES()
  SET_DEFAULT(${PROJECT_NAME}_ENABLE_ALL_PACKAGES OFF)
  SET_DEFAULT(${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE OFF)
  SET(DO_PROCESS_MPI_ENABLES ON) # Should not be needed but CMake is not working!
  TRIBITS_ADJUST_PACKAGE_ENABLES(TRUE)
 
ENDMACRO()


#####################################################################
#
# Set common/base options
#
#####################################################################

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

SET(Trilinos_TPLS_FINDMODS_CLASSIFICATIONS
  MPI            cmake/TPLs/    PS
  BLAS           cmake/TPLs/    PS
  LAPACK         cmake/TPLs/    PS
  Boost          cmake/TPLs/    SS
  )

SET(EXTRA_REPO_NAME extraRepoTwoPackages)
SET(EXTRA_REPO_DIR extraRepoTwoPackages)

INCLUDE(${PROJECT_SOURCE_DIR}/${EXTRA_REPO_NAME}/PackagesList.cmake)

SET(${EXTRA_REPO_NAME}_TPLS_FINDMODS_CLASSIFICATIONS)

SET(${PROJECT_NAME}_ALL_REPOSITORIES "." "${EXTRA_REPO_NAME}")

SET( ${PROJECT_NAME}_ASSERT_MISSING_PACKAGES ON )
