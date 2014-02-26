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


TRIBITS_ADD_ADVANCED_TEST( RawHelloWorld
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  CATEGORIES CONTINUOUS
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/RawHelloWorld
    PASS_REGULAR_EXPRESSION_ALL
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*tribits/doc/examples/UnitTests/TriBITS_RawHelloWorld"
  TEST_1 CMND make
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target hello_world_lib"
      "Built target hello_world"
      "Built target unit_tests"
  TEST_2 CMND ${CMAKE_CTEST_COMMAND}
    PASS_REGULAR_EXPRESSION_ALL
      ": test .*  Passed"
      ": unit_test .*  Passed"
      "100% tests passed, 0 tests failed out of 2"
  )


TRIBITS_ADD_ADVANCED_TEST( TribitsHelloWorld
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  CATEGORIES CONTINUOUS
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsHelloWorld
    PASS_REGULAR_EXPRESSION_ALL
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*tribits/doc/examples/UnitTests/TriBITS_TribitsHelloWorld"
  TEST_1 CMND make
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target hello_world_lib"
      "Built target hello_world"
      "Built target HelloWorld_unit_tests"
  TEST_2 CMND ${CMAKE_CTEST_COMMAND}
    PASS_REGULAR_EXPRESSION_ALL
      ": HelloWorld_hello_world .*  Passed"
      ": HelloWorld_unit_tests .*  Passed"
      "100% tests passed, 0 tests failed out of 2"
  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_PT
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  CATEGORIES CONTINUOUS
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Explicitly enabled packages on input .by user.:  0"
      "Explicitly disabled packages on input .by user or by default.:  WrapExternal 1"
      "Enabling all packages that are not currently disabled because of TribitsExProj_ENABLE_ALL_PACKAGES=ON "
      "Setting TribitsExProj_ENABLE_SimpleCxx=ON"
      "Setting TribitsExProj_ENABLE_MixedLanguage=ON"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackages=ON"
      "Final set of enabled packages:  SimpleCxx MixedLanguage PackageWithSubpackages 3"
      "Final set of enabled SE packages:  SimpleCxx MixedLanguage PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages 6"
      "Final set of non-enabled packages:  WrapExternal 1"
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*tribits/doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_PT"
  TEST_1 CMND make
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target simplecxx"
      "Built target mixedlanguage"
      "Built target pws_a"
      "Built target pws_b"
      "Built target pws_c"
  TEST_2 CMND ${CMAKE_CTEST_COMMAND}
    PASS_REGULAR_EXPRESSION_ALL
      "SimpleCxx_HelloWorldTests .* Passed"
      "MixedLanguage_RayTracerTests .* Passed"
      "PackageWithSubpackagesSubpackageA_test_of_a .* Passed"
      "PackageWithSubpackagesSubpackageB_test_of_b .* Passed"
      "PackageWithSubpackagesSubpackageC_test_of_c .* Passed"
      "100% tests passed, 0 tests failed out of 5"
  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_WrapExternal
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  CATEGORIES CONTINUOUS
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=OFF
      -DTribitsExProj_ENABLE_INSTALL_CMAKE_CONFIG_FILES=OFF
      -DTribitsExProj_ENABLE_WrapExternal=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Explicitly enabled packages on input .by user.:  WrapExternal 1"
      "Explicitly disabled packages on input .by user or by default.:  0"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackages=ON because WrapExternal has an optional dependence on PackageWithSubpackages"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackagesSubpackageA=ON because PackageWithSubpackages has an optional dependence on PackageWithSubpackagesSubpackageA"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackagesSubpackageB=ON because PackageWithSubpackages has an optional dependence on PackageWithSubpackagesSubpackageB"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackagesSubpackageC=ON because PackageWithSubpackages has an optional dependence on PackageWithSubpackagesSubpackageC"
      "Setting TribitsExProj_ENABLE_SimpleCxx=ON because PackageWithSubpackagesSubpackageA has an optional dependence on SimpleCxx"
      "Final set of enabled packages:  SimpleCxx PackageWithSubpackages WrapExternal 3"
      "Final set of enabled SE packages:  SimpleCxx PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages WrapExternal 6"
      "Final set of non-enabled packages:  MixedLanguage 1"
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*tribits/doc/examples/UnitTests/TriBITS_TribitsExampleProject_WrapExternal"
  TEST_1 CMND make
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target WrapExternal_run_external_func"
  TEST_2 CMND ${CMAKE_CTEST_COMMAND}
    PASS_REGULAR_EXPRESSION_ALL
      "WrapExternal_run_external_func .* Passed"
      "100% tests passed, 0 tests failed out of 1"
  )
