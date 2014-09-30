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


SET(COMMON_ENV_ARGS_PASSTHROUGH
  -DTPL_ENABLE_MPI=${TPL_ENABLE_MPI}
  )

IF (TPL_ENABLE_MPI)
  APPEND_SET(COMMON_ENV_ARGS_PASSTHROUGH
    -DMPI_C_COMPILER=${MPI_C_COMPILER}
    -DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}
    -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER}
    -DMPI_EXEC=${MPI_EXEC}
    -DMPI_EXEC_DEFAULT_NUMPROCS=${MPI_EXEC_DEFAULT_NUMPROCS}
    -DMPI_EXEC_MAX_NUMPROCS=${MPI_EXEC_MAX_NUMPROCS}
    -DMPI_EXEC_NUMPROCS_FLAG=${MPI_EXEC_NUMPROCS_FLAG}
    -DMPI_EXEC_PRE_NUMPROCS_FLAGS=${MPI_EXEC_PRE_NUMPROCS_FLAGS}
    -DMPI_EXEC_POST_NUMPROCS_FLAGS=${MPI_EXEC_POST_NUMPROCS_FLAGS}
    )
  SET(TEST_MPI_1_SUFFIX "_MPI_1")
ELSE()
  APPEND_SET(COMMON_ENV_ARGS_PASSTHROUGH
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    )
  SET(TEST_MPI_1_SUFFIX "")
ENDIF()


TRIBITS_ADD_ADVANCED_TEST( RawHelloWorld
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  CATEGORIES CONTINUOUS
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/RawHelloWorld
    PASS_REGULAR_EXPRESSION_ALL
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*doc/examples/UnitTests/TriBITS_RawHelloWorld"
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
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsHelloWorld
    PASS_REGULAR_EXPRESSION_ALL
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsHelloWorld"
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



TRIBITS_ADD_ADVANCED_TEST( TribitsHelloWorld_ScaleTimeout
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  CATEGORIES CONTINUOUS

  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DDART_TESTING_TIMEOUT=200.0
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsHelloWorld
    PASS_REGULAR_EXPRESSION_ALL "Generating done"
    MESSAGE "Configure with default 1.0 scaling"
  TEST_1 CMND grep ARGS "^TimeOut: " DartConfiguration.tcl
    PASS_REGULAR_EXPRESSION "TimeOut: 200.0"
  TEST_2 CMND grep ARGS "^DART_TESTING_TIMEOUT:" CMakeCache.txt
    MESSAGE "DART_TESTING_TIMEOUT in cache does not change!"
    PASS_REGULAR_EXPRESSION "DART_TESTING_TIMEOUT:STRING=200.0"

  TEST_3 CMND ${CMAKE_COMMAND}
    ARGS
      -DTribitsHelloWorld_SCALE_TEST_TIMEOUT=1.5
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsHelloWorld
    MESSAGE "Configure with 1.5 scaling"
    PASS_REGULAR_EXPRESSION_ALL
      "DART_TESTING_TIMEOUT=200.0 being scaled by TribitsHelloWorld_SCALE_TEST_TIMEOUT=1.5 to 300"
      "Generating done"

  TEST_4 CMND grep ARGS "^TimeOut: " DartConfiguration.tcl
    PASS_REGULAR_EXPRESSION "TimeOut: 300"
  TEST_5 CMND grep ARGS "^DART_TESTING_TIMEOUT:" CMakeCache.txt
    MESSAGE "DART_TESTING_TIMEOUT in cache does not change!"
    PASS_REGULAR_EXPRESSION "DART_TESTING_TIMEOUT:STRING=200.0"

  TEST_6 CMND ${CMAKE_COMMAND}
    ARGS
      -DTribitsHelloWorld_SCALE_TEST_TIMEOUT=2.0
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsHelloWorld
    MESSAGE "Configure with 2.0 scaling"
    PASS_REGULAR_EXPRESSION_ALL
      "DART_TESTING_TIMEOUT=200.0 being scaled by TribitsHelloWorld_SCALE_TEST_TIMEOUT=2.0 to 400"
      "Generating done"
  TEST_7 CMND grep ARGS "^TimeOut: " DartConfiguration.tcl
    PASS_REGULAR_EXPRESSION "TimeOut: 400"
  TEST_8 CMND grep ARGS "^DART_TESTING_TIMEOUT:" CMakeCache.txt
    MESSAGE "DART_TESTING_TIMEOUT in cache does not change!"
    PASS_REGULAR_EXPRESSION "DART_TESTING_TIMEOUT:STRING=200.0"

  )


ASSERT_DEFINED(TPL_ENABLE_MPI)
IF (TPL_ENABLE_MPI)
  SET(TPL_MPI_FILE_TRACE
    "-- File Trace: TPL        INCLUDE    .*/tpls/FindTPLMPI.cmake")
ELSE()
  SET(TPL_MPI_FILE_TRACE "")
ENDIF()


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_PT_NoFortran
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  CATEGORIES CONTINUOUS
  TEST_0 CMND ${CMAKE_COMMAND}
    MESSAGE "Do the initial configure"
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=OFF
      -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      -DTribitsExProj_TRACE_FILE_PROCESSING=ON
      -DTribitsExProj_ENABLE_CPACK_PACKAGING=ON
      -DTribitsExProj_DUMP_CPACK_SOURCE_IGNORE_FILES=ON
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Warning: Setting TribitsExProj_ENABLE_WrapExternal=OFF because "
      "Warning: Setting TribitsExProj_ENABLE_MixedLanguage=OFF because "
      "Explicitly enabled packages on input .by user.:  0"
      "Explicitly disabled packages on input .by user or by default.:  MixedLanguage WrapExternal 2"
      "Enabling all SE packages that are not currently disabled because of TribitsExProj_ENABLE_ALL_PACKAGES=ON "
      "Setting TribitsExProj_ENABLE_SimpleCxx=ON"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackages=ON"
      "Final set of enabled packages:  SimpleCxx PackageWithSubpackages 2"
      "Final set of enabled SE packages:  SimpleCxx PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages 5"
      "Final set of non-enabled packages:  MixedLanguage WrapExternal 2"
      "Performing Test HAVE_SIMPLECXX___INT64"
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_PT_NoFortran"
      "-- File Trace: PROJECT    INCLUDE    .*/TribitsExampleProject/Version.cmake"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/cmake/CallbackSetupExtraOptions.cmake"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/PackagesList.cmake"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/TPLsList.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/simple_cxx/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/mixed_language/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/package_with_subpackages/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/package_with_subpackages/A/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/package_with_subpackages/B/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/package_with_subpackages/C/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/wrap_external/cmake/Dependencies.cmake"
      "-- File Trace: PROJECT    CONFIGURE  .*/TribitsExampleProject/cmake/ctest/CTestCustom.cmake.in"
      "-- File Trace: REPOSITORY READ       .*/TribitsExampleProject/Copyright.txt"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/Version.cmake"
      "${TPL_MPI_FILE_TRACE}"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/simple_cxx/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/simple_cxx/test/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/package_with_subpackages/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/package_with_subpackages/A/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/package_with_subpackages/A/tests/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/package_with_subpackages/B/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/package_with_subpackages/B/tests/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/package_with_subpackages/C/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/package_with_subpackages/C/tests/CMakeLists.txt"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/cmake/CallbackDefineRepositoryPackaging.cmake"
      "-- File Trace: PROJECT    INCLUDE    .*/TribitsExampleProject/cmake/CallbackDefineProjectPackaging.cmake"
  TEST_1 CMND make
    MESSAGE "Build the default 'all' target using raw 'make'"
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target simplecxx"
      "Built target pws_a"
      "Built target pws_b"
      "Built target pws_c"
  TEST_2 CMND ${CMAKE_CTEST_COMMAND}
    MESSAGE "Run all the tests with raw 'ctest'"
    PASS_REGULAR_EXPRESSION_ALL
      "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
      "PackageWithSubpackagesSubpackageA_test_of_a .* Passed"
      "PackageWithSubpackagesSubpackageB_test_of_b .* Passed"
      "PackageWithSubpackagesSubpackageC_test_of_c .* Passed"
      "100% tests passed, 0 tests failed out of 4"
  TEST_3 CMND make ARGS package_source
    MESSAGE "Create the tarball"
    PASS_REGULAR_EXPRESSION_ALL
      "Run CPack packaging tool for source..."
      "CPack: Create package using TGZ"
      "CPack: Install projects"
      "CPack: - Install directory: .*/doc/examples/TribitsExampleProject"
      "CPack: Create package"
      "CPack: - package: .*/doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_PT_NoFortran/tribitsexproj-1.1-Source.tar.gz generated."
      "CPack: Create package using TBZ2"
      "CPack: Install projects"
      "CPack: - Install directory: .*/doc/examples/TribitsExampleProject"
      "CPack: Create package"
      "CPack: - package: .*/doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_PT_NoFortran/tribitsexproj-1.1-Source.tar.bz2 generated."
  TEST_4 CMND tar ARGS -xzf tribitsexproj-1.1-Source.tar.gz
    MESSAGE "Untar the tarball"
  TEST_5 CMND diff ARGS -qr ${CMAKE_CURRENT_SOURCE_DIR}/../TribitsExampleProject tribitsexproj-1.1-Source
    MESSAGE "Make sure right directoires are excluced"
    PASS_REGULAR_EXPRESSION_ALL
      "Only in .*/TribitsExampleProject/cmake: ctest"
      "Only in .*/TribitsExampleProject: .gitignore"
      "Only in .*/TribitsExampleProject/packages: mixed_language"
      "Only in .*/TribitsExampleProject/packages: wrap_external"
  )


ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Fortran)
IF (${PROJECT_NAME}_ENABLE_Fortran)
  TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_PT
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1
    CATEGORIES CONTINUOUS
    TEST_0 CMND ${CMAKE_COMMAND}
      MESSAGE "Do the initial configure"
      ARGS
        ${COMMON_ENV_ARGS_PASSTHROUGH}
        -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
        -DTribitsExProj_ENABLE_TESTS=ON
        ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
      PASS_REGULAR_EXPRESSION_ALL
        "Warning: Setting TribitsExProj_ENABLE_WrapExternal=OFF because "
        "Explicitly enabled packages on input .by user.:  0"
        "Explicitly disabled packages on input .by user or by default.:  WrapExternal 1"
        "Enabling all SE packages that are not currently disabled because of TribitsExProj_ENABLE_ALL_PACKAGES=ON "
        "Setting TribitsExProj_ENABLE_SimpleCxx=ON"
        "Setting TribitsExProj_ENABLE_MixedLanguage=ON"
        "Setting TribitsExProj_ENABLE_PackageWithSubpackages=ON"
        "Final set of enabled packages:  SimpleCxx MixedLanguage PackageWithSubpackages 3"
        "Final set of enabled SE packages:  SimpleCxx MixedLanguage PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages 6"
        "Final set of non-enabled packages:  WrapExternal 1"
        "Configuring done"
        "Generating done"
        "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_PT"
    TEST_1 CMND make
      MESSAGE "Build the default 'all' target using raw 'make'"
      ARGS ${CTEST_BUILD_FLAGS}
      PASS_REGULAR_EXPRESSION_ALL
        "Built target simplecxx"
        "Built target mixedlanguage"
        "Built target pws_a"
        "Built target pws_b"
        "Built target pws_c"
    TEST_2 CMND ${CMAKE_CTEST_COMMAND}
      MESSAGE "Run all the tests with raw 'ctest'"
      PASS_REGULAR_EXPRESSION_ALL
        "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
        "MixedLanguage_RayTracerTests${TEST_MPI_1_SUFFIX} .* Passed"
        "PackageWithSubpackagesSubpackageA_test_of_a .* Passed"
        "PackageWithSubpackagesSubpackageB_test_of_b .* Passed"
        "PackageWithSubpackagesSubpackageC_test_of_c .* Passed"
        "100% tests passed, 0 tests failed out of 5"
    )
ENDIF()


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_SimpleCxx_DEBUG_int64
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  CATEGORIES CONTINUOUS
  TEST_0 CMND ${CMAKE_COMMAND}
    MESSAGE "Do the initial configure"
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_DEBUG=ON
      -DTribitsExProj_ENABLE_SimpleCxx=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      -DHAVE_SIMPLECXX___INT64=ON
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Final set of enabled packages:  SimpleCxx 1"
      "Configuring done"
      "Generating done"
  TEST_1 CMND make
    MESSAGE "Build the default 'all' target using raw 'make'"
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target simplecxx"
  TEST_2 CMND ${CMAKE_CTEST_COMMAND}
    MESSAGE "Run all the tests with raw 'ctest'"
    PASS_REGULAR_EXPRESSION_ALL
      "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
      "100% tests passed, 0 tests failed out of 1"
  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_WrapExternal
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  CATEGORIES CONTINUOUS

  TEST_0 CMND cp
    ARGS -r ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
      .

  TEST_1 CMND ${CMAKE_COMMAND}
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_DEBUG=OFF
      -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=OFF
      -DTribitsExProj_ENABLE_INSTALL_CMAKE_CONFIG_FILES=OFF
      -DTribitsExProj_ENABLE_WrapExternal=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Explicitly enabled packages on input .by user.:  WrapExternal 1"
      "Explicitly disabled packages on input .by user or by default.:  0"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackages=ON because WrapExternal has a required dependence on PackageWithSubpackages"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackagesSubpackageA=ON because PackageWithSubpackages has a required dependence on PackageWithSubpackagesSubpackageA"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackagesSubpackageB=ON because PackageWithSubpackages has a required dependence on PackageWithSubpackagesSubpackageB"
      "Setting TribitsExProj_ENABLE_PackageWithSubpackagesSubpackageC=ON because PackageWithSubpackages has a required dependence on PackageWithSubpackagesSubpackageC"
      "Setting TribitsExProj_ENABLE_SimpleCxx=ON because PackageWithSubpackagesSubpackageA has a required dependence on SimpleCxx"
      "Final set of enabled packages:  SimpleCxx PackageWithSubpackages WrapExternal 3"
      "Final set of enabled SE packages:  SimpleCxx PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages WrapExternal 6"
      "Final set of non-enabled packages:  MixedLanguage 1"
      "This package has no unfiltered binary files so consider out of date"
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_WrapExternal"
  TEST_2 CMND make ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target WrapExternal_run_external_func"
  TEST_3 CMND ${CMAKE_CTEST_COMMAND}
    PASS_REGULAR_EXPRESSION_ALL
      "WrapExternal_run_external_func${TEST_MPI_1_SUFFIX} .* Passed"
      "100% tests passed, 0 tests failed out of 1"

  TEST_4 CMND sleep ARGS 1s 
     MESSAGE "Sleep for 1 sec for systems were time stamps are only accurate to 1 sec"
  TEST_5 CMND touch
     ARGS TribitsExampleProject/packages/simple_cxx/src/SimpleCxx_HelloWorld.cpp
     MESSAGE "Test that changing upstream source will trigger rebuild"
  TEST_6 CMND ${CMAKE_COMMAND} ARGS TribitsExampleProject
     MESSAGE "Recofigure with changed upstream source"
    PASS_REGULAR_EXPRESSION_ALL
      "The upstream SE package source file ./src/SimpleCxx_HelloWorld.cpp is more recent than this package's binary file ./WrapExternal_run_external_func.exe"
      "Blowing away WrapExternal build dir external_func/ so it will build from scratch"
  TEST_7 CMND make ARGS ${CTEST_BUILD_FLAGS}
    MESSAGE "Rebuild only exteranl_func"
    PASS_REGULAR_EXPRESSION_ALL
      "Built target simplecxx"
      "Built target pws_a"
      "Built target pws_b"
      "Built target pws_c"
      "Generating external_func/libexternal_func.a"
      "Linking CXX executable WrapExternal_run_external_func.exe"

  TEST_8 CMND ${CMAKE_COMMAND}
     ARGS  -DSimpleCxx_ENABLE_DEBUG=ON  TribitsExampleProject
     MESSAGE "Recofigure changing the debug mode to trigger rebuild"
    PASS_REGULAR_EXPRESSION_ALL
      "The upstream SE package binary file ./src/SimpleCxx_config.h is more recent than this package's binary file ./WrapExternal_run_external_func.exe"
      "Blowing away WrapExternal build dir external_func/ so it will build from scratch"
  TEST_9 CMND make ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Generating external_func/libexternal_func.a"
      "Linking CXX executable WrapExternal_run_external_func.exe"

  TEST_10 CMND sleep ARGS 1s 
     MESSAGE "Sleep for 1 sec for systems were time stamps are only accurate to 1 sec"
  TEST_11 CMND touch
     ARGS TribitsExampleProject/packages/wrap_external/external_func/configure.py
     MESSAGE "Test that changing the external file will trigger rebuild"
  TEST_12 CMND ${CMAKE_COMMAND} ARGS TribitsExampleProject
     MESSAGE "Recofigure with changes external file"
    PASS_REGULAR_EXPRESSION_ALL
      "The this package's source file ./external_func/configure.py is more recent than this package's binary file ./WrapExternal_run_external_func.exe"
      "Blowing away WrapExternal build dir external_func/ so it will build from scratch"
  TEST_13 CMND make ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Generating external_func/libexternal_func.a"
      "Linking CXX executable WrapExternal_run_external_func.exe"

  TEST_14 CMND ${CMAKE_COMMAND} ARGS TribitsExampleProject
     MESSAGE "Recofigure with no changes that will not do anything"
    PASS_REGULAR_EXPRESSION_ALL
      "This package's most recent binary file ./WrapExternal_run_external_func.exe is more recent than its upstream SE package source or binary files or this package's source files"
  TEST_15 CMND make ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target build_external_func"
      "Built target WrapExternal_run_external_func"

  )
