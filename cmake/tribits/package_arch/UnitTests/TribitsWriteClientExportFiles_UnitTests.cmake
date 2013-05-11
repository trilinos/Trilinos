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

MESSAGE("CURRENT_TEST_DIRECTORY = ${CURRENT_TEST_DIRECTORY}")

INCLUDE(${CMAKE_CURRENT_LIST_DIR}/TribitsAdjustPackageEnablesHelpers.cmake)
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsWriteClientExportFiles)


#####################################################################
#
# Unit tests for code in TribitsWriteClientExportFiles.cmake
#
#####################################################################


MACRO(SETUP_WRITE_SPECIALIZED_PACKAGE_EXPORT_MAKEFILE_TEST_STUFF)

  # These would be set the TriBITS env probing code or by CMake
  SET(${PROJECT_NAME}_ENABLE_C ON)
  SET(${PROJECT_NAME}_ENABLE_CXX ON)
  SET(${PROJECT_NAME}_ENABLE_Fortran ON)

  # These would be set automatically by CMake if we were not in script mode!
  SET(CMAKE_LINK_LIBRARY_FLAG -l)
  SET(CMAKE_LIBRARY_PATH_FLAG -L)

ENDMACRO()


#
# A) Test basic package processing and reading dependencies 
#


FUNCTION(UNITTEST_WRITE_SPECIALIZED_PACKAGE_EXPORT_MAKEFILE_RTOP_BEFORE_LIBS)

  MESSAGE("\n***")
  MESSAGE("*** Testing the generation of a specialized export makefile for RTOp *before* libs")
  MESSAGE("***\n")

  SETUP_WRITE_SPECIALIZED_PACKAGE_EXPORT_MAKEFILE_TEST_STUFF()

  # Debugging
  SET(${PROJECT_NAME}_VERBOSE_CONFIGURE ON)
  SET(TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP ON)

  SET(${PROJECT_NAME}_ENABLE_RTOp ON)
  SET(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES ON)

  UNITTEST_HELPER_READ_AND_PROESS_PACKAGES()

  # These vars would be set up by the FindTPL<TPLNAME>.cmake modules if they
  # were called
  SET(TPL_BLAS_LIBRARIES "blaspath/lib/libblas.a")
  SET(TPL_BLAS_LIBRARY_DIRS "blashpath/lib")
  SET(TPL_BLAS_INCLUDE_DIRS "blaspath/include")
  SET(TPL_LAPACK_LIBRARIES "lapackpath/lib/liblapack.a")
  SET(TPL_LAPACK_LIBRARY_DIRS "lapackhpath/lib")
  SET(TPL_LAPACK_INCLUDE_DIRS "lapackhpath/include")

  # These vars should be generated automatically by TRIBITS_PACKAGE() that
  # begins with the upstreams packages.
  SET(Teuchos_LIBRARY_DIRS "teuchos/core/src;teuchos/numeric/src")
  SET(Teuchos_INCLUDE_DIRS "teuchos/core/include;teuchos/numeric/include")
  SET(Teuchos_LIBRARIES "teuchoscore;teuchosnumeric")
  SET(Teuchos_HAS_NATIVE_LIBRARIES TRUE)

  SET(GENERATED_EXPORT_MAKEFILE
    "${CURRENT_TEST_DIRECTORY}/Makefile.export.RTOp.before")

  TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES(
    PACKAGE_NAME RTOp
    EXPORT_FILE_VAR_PREFIX RTOp1
    WRITE_EXPORT_MAKLEFILE "${GENERATED_EXPORT_MAKEFILE}"
    )

  UNITTEST_FILE_REGEX("${GENERATED_EXPORT_MAKEFILE}"
    REGEX_STRINGS
       "RTOp1_INCLUDE_DIRS= -Iteuchos/core/include -Iteuchos/numeric/include"
       "RTOp1_LIBRARY_DIRS= -Lteuchos/core/src -Lteuchos/numeric/src"
       "RTOp1_LIBRARIES= -lteuchoscore -lteuchosnumeric"
       "RTOp1_TPL_INCLUDE_DIRS= -Ilapackhpath/include -Iblaspath/include"
       "RTOp1_TPL_LIBRARIES= -llapackpath/lib/liblapack.a -lblaspath/lib/libblas.a"
       "RTOp1_PACKAGE_LIST= Teuchos"
       "RTOp1_TPL_LIST= LAPACK BLAS"
    )

ENDFUNCTION()


FUNCTION(UNITTEST_WRITE_SPECIALIZED_PACKAGE_EXPORT_MAKEFILE_RTOP_AFTER_LIBS)

  MESSAGE("\n***")
  MESSAGE("*** Testing the generation of a specialized export makefile for RTOp *after* libs")
  MESSAGE("***\n")

  SETUP_WRITE_SPECIALIZED_PACKAGE_EXPORT_MAKEFILE_TEST_STUFF()

  # Debugging
  SET(${PROJECT_NAME}_VERBOSE_CONFIGURE ON)
  SET(TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP ON)

  SET(${PROJECT_NAME}_ENABLE_RTOp ON)
  SET(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES ON)

  UNITTEST_HELPER_READ_AND_PROESS_PACKAGES()

  # These vars would be set up by the FindTPL<TPLNAME>.cmake modules if they
  # were called
  SET(TPL_BLAS_LIBRARIES "blaspath/lib/libblas.a")
  SET(TPL_BLAS_LIBRARY_DIRS "blashpath/lib")
  SET(TPL_BLAS_INCLUDE_DIRS "blaspath/include")
  SET(TPL_LAPACK_LIBRARIES "lapackpath/lib/liblapack.a")
  SET(TPL_LAPACK_LIBRARY_DIRS "lapackhpath/lib")
  SET(TPL_LAPACK_INCLUDE_DIRS "lapackhpath/include")

  # These vars should be generated automatically by TRIBITS_PACKAGE() that
  # begins with the upstreams packages.
  SET(Teuchos_LIBRARY_DIRS "teuchos/core/src;teuchos/numeric/src")
  SET(Teuchos_INCLUDE_DIRS "teuchos/core/include;teuchos/numeric/include")
  SET(Teuchos_LIBRARIES "teuchoscore;teuchosnumeric")
  SET(Teuchos_HAS_NATIVE_LIBRARIES TRUE)
  SET(RTOp_LIBRARY_DIRS "rtop/src;teuchos/core/src;teuchos/numeric/src")
  SET(RTOp_INCLUDE_DIRS "rtop/include;teuchos/core/include;teuchos/numeric/include")
  SET(RTOp_LIBRARIES "rtop")
  SET(RTOp_HAS_NATIVE_LIBRARIES TRUE)

  SET(GENERATED_EXPORT_MAKEFILE
    "${CURRENT_TEST_DIRECTORY}/Makefile.export.RTOp.after")

  TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES(
    PACKAGE_NAME RTOp
    EXPORT_FILE_VAR_PREFIX RTOp2
    WRITE_EXPORT_MAKLEFILE "${GENERATED_EXPORT_MAKEFILE}"
    )

  UNITTEST_FILE_REGEX("${GENERATED_EXPORT_MAKEFILE}"
    REGEX_STRINGS
       "RTOp2_INCLUDE_DIRS= -Irtop/include -Iteuchos/core/include -Iteuchos/numeric/include"
       "RTOp2_LIBRARY_DIRS= -Lrtop/src -Lteuchos/core/src -Lteuchos/numeric/src"
       "RTOp2_LIBRARIES= -lrtop -lteuchoscore -lteuchosnumeric"
       "RTOp2_TPL_INCLUDE_DIRS= -Ilapackhpath/include -Iblaspath/include"
       "RTOp2_TPL_LIBRARIES= -llapackpath/lib/liblapack.a -lblaspath/lib/libblas.a"
       "RTOp2_PACKAGE_LIST= RTOp Teuchos"
       "RTOp2_TPL_LIST= LAPACK BLAS"
    )

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

#
# Run the unit tests
#

UNITTEST_WRITE_SPECIALIZED_PACKAGE_EXPORT_MAKEFILE_RTOP_BEFORE_LIBS()
UNITTEST_WRITE_SPECIALIZED_PACKAGE_EXPORT_MAKEFILE_RTOP_AFTER_LIBS()

# Pass in the number of expected tests that must pass!
UNITTEST_FINAL_RESULT(14)
