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


#
# Common arguments
#


SET(COMMON_ENV_ARGS_PASSTHROUGH
  -DTPL_ENABLE_MPI=${TPL_ENABLE_MPI}
  -DHeaderOnlyTpl_INCLUDE_DIRS=${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/tpls/HeaderOnlyTpl
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


#
# RawHelloWorld
#


TRIBITS_ADD_ADVANCED_TEST( RawHelloWorld
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
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
  TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
    PASS_REGULAR_EXPRESSION_ALL
      ": test .*  Passed"
      ": unit_test .*  Passed"
      "100% tests passed, 0 tests failed out of 2"
  )


#
# TribitsHelloWorld
#


TRIBITS_ADD_ADVANCED_TEST( TribitsHelloWorld
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
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
  TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
    PASS_REGULAR_EXPRESSION_ALL
      ": HelloWorld_hello_world .*  Passed"
      ": HelloWorld_unit_tests .*  Passed"
      "100% tests passed, 0 tests failed out of 2"
  )



TRIBITS_ADD_ADVANCED_TEST( TribitsHelloWorld_ScaleTimeout
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

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


#
# TribitsExampleProject
#


ASSERT_DEFINED(TPL_ENABLE_MPI)
IF (TPL_ENABLE_MPI)
  SET(TPL_MPI_FILE_TRACE
    "-- File Trace: TPL        INCLUDE    .*/tpls/FindTPLMPI.cmake")
  SET(FINAL_ENABLED_TPLS "MPI HeaderOnlyTpl 2")
ELSE()
  SET(TPL_MPI_FILE_TRACE "")
  SET(FINAL_ENABLED_TPLS "HeaderOnlyTpl 1")
ENDIF()


IF (EXISTS "${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject/.gitignore")
  SET(REGEX_FOR_GITIGNORE "Only in .*/TribitsExampleProject: .gitignore")
ELSE()
  SET(REGEX_FOR_GITIGNORE)
ENDIF()


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_ST_NoFortran
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0 CMND ${CMAKE_COMMAND}
    MESSAGE "Do the initial configure (and test a lot of things at once)"
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=OFF
      -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON
      -DTribitsExProj_TRACE_FILE_PROCESSING=ON
      -DTribitsExProj_ENABLE_CPACK_PACKAGING=ON
      -DTribitsExProj_DUMP_CPACK_SOURCE_IGNORE_FILES=ON
      -DTribitsExProj_DUMP_PACKAGE_DEPENDENCIES=ON
      -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=ON
      -DTribitsExProj_ENABLE_INSTALL_CMAKE_CONFIG_FILES=ON
      -DCMAKE_INSTALL_PREFIX=install
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Warning: Setting TribitsExProj_ENABLE_WrapExternal=OFF because "
      "Warning: Setting TribitsExProj_ENABLE_MixedLang=OFF because "
      "Printing package dependencies ..."
      "-- TribitsExProj_PACKAGES: SimpleCxx MixedLang WithSubpackages WrapExternal"
      "-- TribitsExProj_SE_PACKAGES: SimpleCxx MixedLang WithSubpackagesA WithSubpackagesB WithSubpackagesC WithSubpackages WrapExternal"

      "-- SimpleCxx_LIB_REQUIRED_DEP_TPLS: HeaderOnlyTpl"
      "-- MixedLang: No dependencies!"
      "-- WithSubpackagesA_LIB_REQUIRED_DEP_PACKAGES: SimpleCxx"
      "-- WithSubpackagesB_LIB_REQUIRED_DEP_PACKAGES: SimpleCxx"
      "-- WithSubpackagesB_LIB_OPTIONAL_DEP_PACKAGES: WithSubpackagesA"
      "-- WithSubpackagesB_TEST_OPTIONAL_DEP_PACKAGES: MixedLang"
      "-- WithSubpackagesC_LIB_REQUIRED_DEP_PACKAGES: WithSubpackagesA WithSubpackagesB"
      "-- WithSubpackages_LIB_REQUIRED_DEP_PACKAGES: WithSubpackagesA"
      "-- WithSubpackages_LIB_OPTIONAL_DEP_PACKAGES: WithSubpackagesB WithSubpackagesC"
      "-- WrapExternal_LIB_REQUIRED_DEP_PACKAGES: WithSubpackagesA"
      "-- WrapExternal_LIB_OPTIONAL_DEP_PACKAGES: MixedLang"
      "-- SimpleCxx: No library dependencies!"
      "-- WithSubpackagesA_FULL_ENABLED_DEP_PACKAGES: SimpleCxx"
      "-- WithSubpackagesB_FULL_ENABLED_DEP_PACKAGES: WithSubpackagesA SimpleCxx"
      "-- WithSubpackagesC_FULL_ENABLED_DEP_PACKAGES: WithSubpackagesB WithSubpackagesA SimpleCxx"
      "-- WithSubpackages_FULL_ENABLED_DEP_PACKAGES: WithSubpackagesC WithSubpackagesB WithSubpackagesA SimpleCxx"
      "Explicitly enabled packages on input .by user.:  0"
      "Explicitly disabled packages on input .by user or by default.:  MixedLang WrapExternal 2"
      "Enabling all SE packages that are not currently disabled because of TribitsExProj_ENABLE_ALL_PACKAGES=ON "
      "Setting TribitsExProj_ENABLE_SimpleCxx=ON"
      "Setting TribitsExProj_ENABLE_WithSubpackages=ON"
      "Setting TPL_ENABLE_HeaderOnlyTpl=ON because it is required by the enabled package SimpleCxx" 
      "Final set of enabled packages:  SimpleCxx WithSubpackages 2"
      "Final set of enabled SE packages:  SimpleCxx WithSubpackagesA WithSubpackagesB WithSubpackagesC WithSubpackages 5"
      "Final set of enabled TPLs:  ${FINAL_ENABLED_TPLS}"
      "Final set of non-enabled packages:  MixedLang WrapExternal 2"
      "Processing enabled TPL: HeaderOnlyTpl"
      "-- File Trace: TPL        INCLUDE    .+/TribitsExampleProject/cmake/tpls/FindTPLHeaderOnlyTpl.cmake"
      "-- TPL_HeaderOnlyTpl_INCLUDE_DIRS='.+/doc/examples/tpls/HeaderOnlyTpl'"
      "Performing Test HAVE_SIMPLECXX___INT64"
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_ST_NoFortran"
      "-- File Trace: PROJECT    INCLUDE    .*/TribitsExampleProject/Version.cmake"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/cmake/CallbackSetupExtraOptions.cmake"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/PackagesList.cmake"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/TPLsList.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/simple_cxx/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/mixed_lang/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/with_subpackages/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/with_subpackages/a/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/with_subpackages/b/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/with_subpackages/c/cmake/Dependencies.cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TribitsExampleProject/packages/wrap_external/cmake/Dependencies.cmake"
      "-- File Trace: PROJECT    CONFIGURE  .*/TribitsExampleProject/cmake/ctest/CTestCustom.cmake.in"
      "-- File Trace: REPOSITORY READ       .*/TribitsExampleProject/Copyright.txt"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/Version.cmake"
      "${TPL_MPI_FILE_TRACE}"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/simple_cxx/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/simple_cxx/test/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/with_subpackages/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/with_subpackages/a/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/with_subpackages/a/tests/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/with_subpackages/b/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/with_subpackages/b/tests/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/with_subpackages/c/CMakeLists.txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TribitsExampleProject/packages/with_subpackages/c/tests/CMakeLists.txt"
      "-- File Trace: REPOSITORY INCLUDE    .*/TribitsExampleProject/cmake/CallbackDefineRepositoryPackaging.cmake"
      "-- File Trace: PROJECT    INCLUDE    .*/TribitsExampleProject/cmake/CallbackDefineProjectPackaging.cmake"
  TEST_1 CMND make ARGS ${CTEST_BUILD_FLAGS}
    MESSAGE "Build the default 'all' target using raw 'make'"
    PASS_REGULAR_EXPRESSION_ALL
      "Built target simplecxx"
      "Built target pws_a"
      "Built target pws_b"
      "Built target pws_c"
  TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
    MESSAGE "Run all the tests with raw 'ctest'"
    PASS_REGULAR_EXPRESSION_ALL
      "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
      "WithSubpackagesA_test_of_a .* Passed"
      "WithSubpackagesB_test_of_b .* Passed"
      "WithSubpackagesC_test_of_c .* Passed"
      "WithSubpackagesC_test_of_c_util.* Passed"
      "100% tests passed, 0 tests failed out of 5"

  TEST_3 CMND make ARGS install ${CTEST_BUILD_FLAGS}
    MESSAGE "Build 'install' target using raw 'make'"
    PASS_REGULAR_EXPRESSION_ALL
      "Installing: .+/install/include/TribitsExProj_version.h"
      "Installing: .+/install/lib/cmake/TribitsExProj/TribitsExProjConfig.cmake"
      "Installing: .+/install/include/Makefile.export.TribitsExProj"
      "Installing: .+/install/lib/cmake/TribitsExProj/TribitsExProjConfigVersion.cmake"
      "Installing: .+/install/include/TribitsExProjConfig.cmake"
      "Installing: .+/install/lib/cmake/SimpleCxx/SimpleCxxConfig.cmake"
      "Installing: .+/install/lib/cmake/SimpleCxx/SimpleCxxTargets.cmake"
      "Installing: .+/install/lib/cmake/SimpleCxx/SimpleCxxTargets-release.cmake"
      "Installing: .+/install/include/Makefile.export.SimpleCxx"
      "Installing: .+/install/lib/libsimplecxx.a"
      "Installing: .+/install/include/SimpleCxx_HelloWorld.hpp"
      "Installing: .+/install/lib/cmake/WithSubpackages/WithSubpackagesConfig.cmake"
      "Installing: .+/install/include/Makefile.export.WithSubpackages"
      "Installing: .+/install/lib/libpws_a.a"
      "Installing: .+/install/include/A.hpp"
      "Installing: .+/install/lib/cmake/WithSubpackagesA/WithSubpackagesAConfig.cmake"
      "Installing: .+/install/lib/cmake/WithSubpackagesA/WithSubpackagesATargets.cmake"
      "Installing: .+/install/lib/cmake/WithSubpackagesA/WithSubpackagesATargets-release.cmake"
      "Installing: .+/install/include/Makefile.export.WithSubpackagesA"
      "Installing: .+/install/lib/libpws_b.a"
      "Installing: .+/install/include/B.hpp"
      "Installing: .+/install/lib/cmake/WithSubpackagesB/WithSubpackagesBConfig.cmake"
      "Installing: .+/install/lib/cmake/WithSubpackagesB/WithSubpackagesBTargets.cmake"
      "Installing: .+/install/lib/cmake/WithSubpackagesB/WithSubpackagesBTargets-release.cmake"
      "Installing: .+/install/include/Makefile.export.WithSubpackagesB"
      "Installing: .+/install/lib/libpws_c.a"
      "Installing: .+/install/include/C.hpp"
      "Installing: .+/install/lib/cmake/WithSubpackagesC/WithSubpackagesCConfig.cmake"
      "Installing: .+/install/lib/cmake/WithSubpackagesC/WithSubpackagesCTargets.cmake"
      "Installing: .+/install/lib/cmake/WithSubpackagesC/WithSubpackagesCTargets-release.cmake"
      "Installing: .+/install/include/Makefile.export.WithSubpackagesC"

  TEST_4 CMND ${CMAKE_COMMAND}
    ARGS
      -DDUMMY_PROJECT_NAME=DummyProject
      -DDUMMY_PROJECT_DIR=dummy_client_of_WithSubpackages
      -DEXPORT_VAR_PREFIX=WithSubpackages
      -DEXPORT_CONFIG_FILE=../install/lib/cmake/WithSubpackages/WithSubpackagesConfig.cmake
      -DCMAKE_COMMAND=${CMAKE_COMMAND}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/RunDummyPackageClientBulid.cmake
    MESSAGE "Create and configure a dummy project that includes WithSubpackagesConfig.cmake"
    PASS_REGULAR_EXPRESSION_ALL
      "WithSubpackages_CXX_FLAGS = '.*'"
      "WithSubpackages_C_FLAGS = '.*'"
      "WithSubpackages_FORTRAN_FLAGS = '.*'"
      "WithSubpackages_EXTRA_LD_FLAGS = '.*'"
      "WithSubpackages_SHARED_LIB_RPATH_COMMAND = '.*'"
      "WithSubpackages_BUILD_SHARED_LIBS = '.*'"
      "WithSubpackages_LINKER = '.+'"
      "WithSubpackages_AR = '.+'"
      "WithSubpackages_INCLUDE_DIRS = '.+/install/lib/cmake/WithSubpackages/../../../include'"
      "WithSubpackages_LIBRARY_DIRS = '.+/install/lib/cmake/WithSubpackages/../../../lib'"
      "WithSubpackages_LIBRARIES = 'pws_c.pws_b.pws_a.simplecxx'"
      "WithSubpackages_TPL_INCLUDE_DIRS = '.+/doc/examples/tpls/HeaderOnlyTpl'"
      "WithSubpackages_TPL_LIBRARY_DIRS = ''"
      "WithSubpackages_TPL_LIBRARIES = ''"
      "WithSubpackages_MPI_LIBRARIES = ''"
      "WithSubpackages_MPI_LIBRARY_DIRS = ''"
      "WithSubpackages_MPI_INCLUDE_DIRS = ''"
      "WithSubpackages_MPI_EXEC = '${MPI_EXEC}'"
      "WithSubpackages_MPI_EXEC_MAX_NUMPROCS = '${MPI_EXEC_MAX_NUMPROCS}'"
      "WithSubpackages_MPI_EXEC_NUMPROCS_FLAG = '${MPI_EXEC_NUMPROCS_FLAG}'"
      "WithSubpackages_PACKAGE_LIST = 'WithSubpackagesC.WithSubpackagesB.WithSubpackagesA.SimpleCxx'"
      "WithSubpackages_TPL_LIST = 'HeaderOnlyTpl'"

  TEST_5 CMND make ARGS package_source
    MESSAGE "Create the tarball"
    PASS_REGULAR_EXPRESSION_ALL
      "Run CPack packaging tool for source..."
      "CPack: Create package using TGZ"
      "CPack: Install projects"
      "CPack: - Install directory: .*/doc/examples/TribitsExampleProject"
      "CPack: Create package"
      "CPack: - package: .*/doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_ST_NoFortran/tribitsexproj-1.1-Source.tar.gz generated."
      "CPack: Create package using TBZ2"
      "CPack: Install projects"
      "CPack: - Install directory: .*/doc/examples/TribitsExampleProject"
      "CPack: Create package"
      "CPack: - package: .*/doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_ST_NoFortran/tribitsexproj-1.1-Source.tar.bz2 generated."
  TEST_6 CMND tar ARGS -xzf tribitsexproj-1.1-Source.tar.gz
    MESSAGE "Untar the tarball"
  TEST_7 CMND diff ARGS -qr ${CMAKE_CURRENT_SOURCE_DIR}/../TribitsExampleProject tribitsexproj-1.1-Source
    MESSAGE "Make sure right directoires are excluced"
    PASS_REGULAR_EXPRESSION_ALL
      "Only in .*/TribitsExampleProject/cmake: ctest"
      ${REGEX_FOR_GITIGNORE}
      "Only in .*/TribitsExampleProject/packages: mixed_lang"
      "Only in .*/TribitsExampleProject/packages: wrap_external"

  )


IF (NOT ${PROJECT_NAME}_HOSTTYPE STREQUAL "Windows")

  TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_PT_NoFortran_ConfigTiming
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1
  
    TEST_0 CMND ${CMAKE_COMMAND}
      MESSAGE "Do the initial configure with basic configure timing"
      ARGS
        ${COMMON_ENV_ARGS_PASSTHROUGH}
        -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -DTribitsExProj_ENABLE_Fortran=OFF
        -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
        -DTribitsExProj_ENABLE_TESTS=ON
        -DTribitsExProj_ENABLE_CPACK_PACKAGING=ON
        -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=ON
        -DTribitsExProj_ENABLE_CONFIGURE_TIMING=ON
        ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
      PASS_REGULAR_EXPRESSION_ALL
        "Total time to read in and process all package dependencies: "
        "Total time to adjust package and TPL enables: "
        "Total time to probe and setup the environment: "
        "Total time to configure enabled TPLs: "
        "Total time to configure enabled packages: "
        "Total time to set up for CPack packaging: "
        "Total time to configure TribitsExProj: "
        "Final set of enabled packages:  SimpleCxx WithSubpackages 2"
        "Final set of enabled SE packages:  SimpleCxx WithSubpackagesA WithSubpackages 3"
  
    TEST_1 CMND ${CMAKE_COMMAND}
      MESSAGE "Reconfigure to test out timing of all packages"
      ARGS
        -DTribitsExProj_ENABLE_PACKAGE_CONFIGURE_TIMING=ON
        .
      PASS_REGULAR_EXPRESSION_ALL
        "Total time to read in and process all package dependencies: "
        "Total time to adjust package and TPL enables: "
        "Total time to probe and setup the environment: "
        "Total time to configure enabled TPLs: "
        "-- Total time to configure package SimpleCxx: "
        "-- Total time to configure package WithSubpackages: "
        "Total time to configure enabled packages: "
        "Total time to set up for CPack packaging: "
        "Total time to configure TribitsExProj: "
  
    TEST_2 CMND ${CMAKE_COMMAND}
      MESSAGE "Reconfigure to test out timing of just one package"
      ARGS
        -DTribitsExProj_ENABLE_PACKAGE_CONFIGURE_TIMING=OFF
        -DSimpleCxx_PACKAGE_CONFIGURE_TIMING=ON
        .
      PASS_REGULAR_EXPRESSION_ALL
        "Total time to read in and process all package dependencies: "
        "Total time to adjust package and TPL enables: "
        "Total time to probe and setup the environment: "
        "Total time to configure enabled TPLs: "
        "-- Total time to configure package SimpleCxx: "
        "Total time to configure enabled packages: "
        "Total time to set up for CPack packaging: "
        "Total time to configure TribitsExProj: "
  
    )

ENDIF()


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_PT_NoFortran
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0 CMND ${CMAKE_COMMAND}
    MESSAGE "Do the initial configure (and test a lot of things at once)"
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=OFF
      -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      -DTribitsExProj_DUMP_PACKAGE_DEPENDENCIES=ON
      -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=ON
      -DTribitsExProj_ENABLE_INSTALL_CMAKE_CONFIG_FILES=ON
      -DCMAKE_INSTALL_PREFIX=install
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "-- SimpleCxx: No library dependencies!"
      "-- WithSubpackagesA_FULL_ENABLED_DEP_PACKAGES: SimpleCxx"
      "-- WithSubpackages_FULL_ENABLED_DEP_PACKAGES: WithSubpackagesA SimpleCxx"
      "Explicitly enabled packages on input .by user.:  0"
      "Explicitly disabled packages on input .by user or by default.:  MixedLang WrapExternal 2"
      "Final set of enabled packages:  SimpleCxx WithSubpackages 2"
      "Final set of enabled SE packages:  SimpleCxx WithSubpackagesA WithSubpackages 3"
      "Final set of non-enabled packages:  MixedLang WrapExternal 2"
      "Final set of non-enabled SE packages:  MixedLang WithSubpackagesB WithSubpackagesC WrapExternal 4"

  TEST_1 CMND make ARGS ${CTEST_BUILD_FLAGS}
    MESSAGE "Build the default 'all' target using raw 'make'"
    PASS_REGULAR_EXPRESSION_ALL
      "Built target simplecxx"
      "Built target pws_a"

  TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
    MESSAGE "Run all the tests with raw 'ctest'"
    PASS_REGULAR_EXPRESSION_ALL
      "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
      "WithSubpackagesA_test_of_a .* Passed"
      "100% tests passed, 0 tests failed out of 2"

  TEST_3 CMND make ARGS install ${CTEST_BUILD_FLAGS}
    MESSAGE "Build 'install' target using raw 'make'"
    PASS_REGULAR_EXPRESSION_ALL
      "Installing: .+/install/lib/cmake/WithSubpackages/WithSubpackagesConfig.cmake"
      "Installing: .+/install/include/Makefile.export.WithSubpackages"

  TEST_4 CMND ${CMAKE_COMMAND}
    ARGS
      -DDUMMY_PROJECT_NAME=DummyProject
      -DDUMMY_PROJECT_DIR=dummy_client_of_WithSubpackages
      -DEXPORT_VAR_PREFIX=WithSubpackages
      -DEXPORT_CONFIG_FILE=../install/lib/cmake/WithSubpackages/WithSubpackagesConfig.cmake
      -DCMAKE_COMMAND=${CMAKE_COMMAND}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/RunDummyPackageClientBulid.cmake
    MESSAGE "Create and configure a dummy project that includes WithSubpackagesConfig.cmake"
    PASS_REGULAR_EXPRESSION_ALL
      "WithSubpackages_INCLUDE_DIRS = '.+/install/lib/cmake/WithSubpackages/../../../include'"
      "WithSubpackages_LIBRARY_DIRS = '.+/install/lib/cmake/WithSubpackages/../../../lib'"
      "WithSubpackages_LIBRARIES = 'pws_a.simplecxx'"
      "WithSubpackages_PACKAGE_LIST = 'WithSubpackagesA.SimpleCxx'"
  )


ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Fortran)
IF (${PROJECT_NAME}_ENABLE_Fortran)


  TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_ST
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1
    TEST_0 CMND ${CMAKE_COMMAND}
      MESSAGE "Do the initial configure"
      ARGS
        ${COMMON_ENV_ARGS_PASSTHROUGH}
        -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -DTribitsExProj_ENABLE_Fortran=ON
        -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
        -DTribitsExProj_ENABLE_TESTS=ON
        -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON
        -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=ON
        -DTribitsExProj_ENABLE_INSTALL_CMAKE_CONFIG_FILES=ON
        ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
      PASS_REGULAR_EXPRESSION_ALL
        "Warning: Setting TribitsExProj_ENABLE_WrapExternal=OFF because "
        "Explicitly enabled packages on input .by user.:  0"
        "Explicitly disabled packages on input .by user or by default.:  WrapExternal 1"
        "Enabling all SE packages that are not currently disabled because of TribitsExProj_ENABLE_ALL_PACKAGES=ON "
        "Setting TribitsExProj_ENABLE_SimpleCxx=ON"
        "Setting TribitsExProj_ENABLE_MixedLang=ON"
        "Setting TribitsExProj_ENABLE_WithSubpackages=ON"
        "Final set of enabled packages:  SimpleCxx MixedLang WithSubpackages 3"
        "Final set of enabled SE packages:  SimpleCxx MixedLang WithSubpackagesA WithSubpackagesB WithSubpackagesC WithSubpackages 6"
        "Final set of non-enabled packages:  WrapExternal 1"
        "Configuring done"
        "Generating done"
        "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_ST"
    TEST_1 CMND make
      MESSAGE "Build the default 'all' target using raw 'make'"
      ARGS ${CTEST_BUILD_FLAGS}
      PASS_REGULAR_EXPRESSION_ALL
        "Built target simplecxx"
        "Built target mixedlang"
        "Built target pws_a"
        "Built target pws_b"
        "Built target pws_c"
    TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
      MESSAGE "Run all the tests with raw 'ctest'"
      PASS_REGULAR_EXPRESSION_ALL
        "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
        "MixedLang_RayTracerTests${TEST_MPI_1_SUFFIX} .* Passed"
        "WithSubpackagesA_test_of_a .* Passed"
        "WithSubpackagesB_test_of_b .* Passed"
        "WithSubpackagesB_test_of_b_mixed_lang.* Passed"
        "WithSubpackagesC_test_of_c_util.* Passed" 
        "WithSubpackagesC_test_of_c .* Passed"
        "WithSubpackagesC_test_of_c_b_mixed_lang.* Passed"
        "100% tests passed, 0 tests failed out of 8"
    )


  TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_ST_LibPrefix
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1
    TEST_0 CMND ${CMAKE_COMMAND}
      MESSAGE "Do the initial configure"
      ARGS
        ${COMMON_ENV_ARGS_PASSTHROUGH}
        -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -DTribitsExProj_ENABLE_Fortran=ON
        -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
        -DTribitsExProj_ENABLE_TESTS=ON
        -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON
        -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=ON
        -DTribitsExProj_ENABLE_INSTALL_CMAKE_CONFIG_FILES=ON
        -DTribitsExProj_LIBRARY_NAME_PREFIX=tep_
        ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
      PASS_REGULAR_EXPRESSION_ALL
        "Warning: Setting TribitsExProj_ENABLE_WrapExternal=OFF because "
        "Explicitly enabled packages on input .by user.:  0"
        "Explicitly disabled packages on input .by user or by default.:  WrapExternal 1"
        "Enabling all SE packages that are not currently disabled because of TribitsExProj_ENABLE_ALL_PACKAGES=ON "
        "Setting TribitsExProj_ENABLE_SimpleCxx=ON"
        "Setting TribitsExProj_ENABLE_MixedLang=ON"
        "Setting TribitsExProj_ENABLE_WithSubpackages=ON"
        "Final set of enabled packages:  SimpleCxx MixedLang WithSubpackages 3"
        "Final set of enabled SE packages:  SimpleCxx MixedLang WithSubpackagesA WithSubpackagesB WithSubpackagesC WithSubpackages 6"
        "Final set of non-enabled packages:  WrapExternal 1"
        "Configuring done"
        "Generating done"
        "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_ST_LibPrefix"
    TEST_1 CMND make
      MESSAGE "Build the default 'all' target using raw 'make'"
      ARGS ${CTEST_BUILD_FLAGS}
      PASS_REGULAR_EXPRESSION_ALL
        "Built target tep_simplecxx"
        "Built target tep_mixedlang"
        "Built target tep_pws_a"
        "Built target tep_pws_b"
        "Built target tep_pws_c"
    TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
      MESSAGE "Run all the tests with raw 'ctest'"
      PASS_REGULAR_EXPRESSION_ALL
        "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
        "MixedLang_RayTracerTests${TEST_MPI_1_SUFFIX} .* Passed"
        "WithSubpackagesA_test_of_a .* Passed"
        "WithSubpackagesB_test_of_b .* Passed"
        "WithSubpackagesB_test_of_b_mixed_lang.* Passed"
        "WithSubpackagesC_test_of_c_util.* Passed" 
        "WithSubpackagesC_test_of_c .* Passed"
        "WithSubpackagesC_test_of_c_b_mixed_lang.* Passed"
        "100% tests passed, 0 tests failed out of 8"
    )


  TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_ST_LibUsage
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1

    TEST_0 
      MESSAGE "Do the initial configure to get the package enables"
        " and the  env probe in place."
      CMND ${CMAKE_COMMAND}
      ARGS
        ${COMMON_ENV_ARGS_PASSTHROUGH}
        -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -DTribitsExProj_ENABLE_Fortran=ON
        -DTribitsExProj_ENABLE_WrapExternal=OFF
        -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
        -DTribitsExProj_ENABLE_TESTS=ON
        -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON
        ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
      PASS_REGULAR_EXPRESSION_ALL
        "Final set of enabled packages:  SimpleCxx MixedLang WithSubpackages 3"
        "Final set of enabled SE packages:  SimpleCxx MixedLang WithSubpackagesA WithSubpackagesB WithSubpackagesC WithSubpackages 6"
        "Configuring done"
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    #
    # Testing pasing libs to TRIBITS_ADD_LIBRARY()
    #

    TEST_1
      MESSAGE "Show deprecated warning when trying to link lib from upstream SE package."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_UPSTREAM_DEPLIBS_ERROR=ON
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: 'simplecxx' in DEPSLIBS is not a lib in this SE package"
        "packages/with_subpackages/b/src/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    TEST_2
      MESSAGE "Show deprecated warning when passing a lib from this SE package through IMPORTEDLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_SE_PKG_LIB_IMPORTEDLIBS_ERROR=ON
        -DSPKB_SHOW_UPSTREAM_DEPLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Lib 'pws_b' in IMPORTEDLIBS is in this SE package "
        "packages/with_subpackages/b/tests/testlib/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    TEST_3
      MESSAGE "Show deprecated warning when passing a lib from upstream"
        " SE package through IMPORTEDLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_UPSTREAM_SE_PKG_LIB_IMPORTEDLIBS_ERROR=ON
        -DSPKB_SHOW_SE_PKG_LIB_IMPORTEDLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Lib 'simplecxx' being passed through IMPORTEDLIBS"
        "TribitsExampleProject/packages/with_subpackages/b/cmake/Dependencies.cmake"
        "packages/with_subpackages/b/tests/testlib/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    TEST_4
      MESSAGE "Show deprecated warning when passing a TESTONLY lib through DEPLBIS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKC_SHOW_TESTONLY_DEPLBIS_ERROR=ON
        -DSPKB_SHOW_UPSTREAM_SE_PKG_LIB_IMPORTEDLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: 'b_mixed_lang' in DEPLIBS is a TESTONLY lib "
        "TribitsExampleProject/packages/with_subpackages/c/cmake/Dependencies.cmake"
        "packages/with_subpackages/c/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    TEST_5
      MESSAGE "Show deprecated warning when passing a TESTONLY lib through IMPORTEDLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKC_SHOW_TESTONLY_IMPORTEDLIBS_ERROR=ON
        -DSPKC_SHOW_TESTONLY_DEPLBIS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: 'b_mixed_lang' in IMPORTEDLIBS is a TESTONLY lib"
        "packages/with_subpackages/c/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    TEST_6
      MESSAGE "Show deprecated warning when passing m through DEPLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_M_DEPLIBS_ERROR=ON
        -DSPKC_SHOW_TESTONLY_IMPORTEDLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: 'm' in DEPSLIBS is not a lib defined in the current cmake "
        "packages/with_subpackages/b/src/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    #
    # Testing pasing libs to TRIBITS_ADD_EXECUTABLE()
    #

    TEST_7
      MESSAGE "Show error when trying to link an INSTALLABLE exec against a TESTONLY lib."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_TESTONLY_INSTALLABLE_ERROR=ON
        -DSPKB_SHOW_M_DEPLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "ERROR: TESTONLY lib 'b_test_utils' not allowed with INSTALLABLE executable"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_8
      MESSAGE "Show error when trying to link against non-TESTONLY lib."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_NON_TESTONLY_LIB_ERROR=ON
        -DSPKB_SHOW_TESTONLY_INSTALLABLE_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        " ERROR: 'simplecxx' in TESTONLYLIBS not a TESTONLY lib"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_9
      MESSAGE "Show error when trying to link SE package lib using TESTONLYLIBS"
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_IMPORTED_LIBS_THIS_PKG_ERROR=ON
        -DSPKB_SHOW_NON_TESTONLY_LIB_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "ERROR: Lib 'pws_b' in IMPORTEDLIBS is in this SE package"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_10
      MESSAGE "Show deprecated warning when trying to link TESTONLY lib using DEBLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_TESTONLY_DEBLIBS_WARNING=ON
        -DSPKB_SHOW_IMPORTED_LIBS_THIS_PKG_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Passing TESTONLY lib 'b_mixed_lang' through DEPLIBS is deprecated"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    TEST_11
      MESSAGE "Show deprecated warning when trying to link non-TESTONLY lib using DEBLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_NONTESTONLY_DEBLIBS_WARNING=ON
        -DSPKB_SHOW_TESTONLY_DEBLIBS_WARNING=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Passing non-TESTONLY lib 'pws_b' through DEPLIBS is deprecated"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    TEST_12
      MESSAGE "Show deprecated warning when trying to link external lib using DEBLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_EXTERNAL_DEBLIBS_WARNING=ON
        -DSPKB_SHOW_NONTESTONLY_DEBLIBS_WARNING=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Passing external lib 'm' through DEPLIBS is deprecated"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage"

    )


  TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_ST_LibUsage_LibPrefix
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1

    TEST_0 
      MESSAGE "Do the initial configure to get the package enables"
        " and the  env probe in place."
      CMND ${CMAKE_COMMAND}
      ARGS
        ${COMMON_ENV_ARGS_PASSTHROUGH}
        -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -DTribitsExProj_ENABLE_Fortran=ON
        -DTribitsExProj_ENABLE_WrapExternal=OFF
        -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
        -DTribitsExProj_ENABLE_TESTS=ON
        -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON
        -DTribitsExProj_LIBRARY_NAME_PREFIX=tep_
        ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
      PASS_REGULAR_EXPRESSION_ALL
        "Final set of enabled packages:  SimpleCxx MixedLang WithSubpackages 3"
        "Final set of enabled SE packages:  SimpleCxx MixedLang WithSubpackagesA WithSubpackagesB WithSubpackagesC WithSubpackages 6"
        "Configuring done"
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage_LibPrefix"

    #
    # Testing pasing libs to TRIBITS_ADD_LIBRARY()
    #

    TEST_1
      MESSAGE "Show deprecated warning when trying to link lib from upstream SE package."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_UPSTREAM_DEPLIBS_ERROR=ON
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: 'simplecxx' in DEPSLIBS is not a lib in this SE package"
        "packages/with_subpackages/b/src/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage_LibPrefix"

    TEST_2
      MESSAGE "Show deprecated warning when passing a lib from this SE package through IMPORTEDLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_SE_PKG_LIB_IMPORTEDLIBS_ERROR=ON
        -DSPKB_SHOW_UPSTREAM_DEPLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Lib 'pws_b' in IMPORTEDLIBS is in this SE package "
        "packages/with_subpackages/b/tests/testlib/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage_LibPrefix"

    TEST_3
      MESSAGE "Show deprecated warning when passing a lib from upstream"
        " SE package through IMPORTEDLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_UPSTREAM_SE_PKG_LIB_IMPORTEDLIBS_ERROR=ON
        -DSPKB_SHOW_SE_PKG_LIB_IMPORTEDLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Lib 'simplecxx' being passed through IMPORTEDLIBS"
        "TribitsExampleProject/packages/with_subpackages/b/cmake/Dependencies.cmake"
        "packages/with_subpackages/b/tests/testlib/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage_LibPrefix"

    TEST_4
      MESSAGE "Show deprecated warning when passing a TESTONLY lib through DEPLBIS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKC_SHOW_TESTONLY_DEPLBIS_ERROR=ON
        -DSPKB_SHOW_UPSTREAM_SE_PKG_LIB_IMPORTEDLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: 'b_mixed_lang' in DEPLIBS is a TESTONLY lib "
        "TribitsExampleProject/packages/with_subpackages/c/cmake/Dependencies.cmake"
        "packages/with_subpackages/c/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage_LibPrefix"

    TEST_5
      MESSAGE "Show deprecated warning when passing a TESTONLY lib through IMPORTEDLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKC_SHOW_TESTONLY_IMPORTEDLIBS_ERROR=ON
        -DSPKC_SHOW_TESTONLY_DEPLBIS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: 'b_mixed_lang' in IMPORTEDLIBS is a TESTONLY lib"
        "packages/with_subpackages/c/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage_LibPrefix"

    TEST_6
      MESSAGE "Show deprecated warning when passing m through DEPLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_M_DEPLIBS_ERROR=ON
        -DSPKC_SHOW_TESTONLY_IMPORTEDLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: 'm' in DEPSLIBS is not a lib defined in the current cmake "
        "packages/with_subpackages/b/src/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage_LibPrefix"

    #
    # Testing pasing libs to TRIBITS_ADD_EXECUTABLE()
    #

    TEST_7
      MESSAGE "Show error when trying to link an INSTALLABLE exec against a TESTONLY lib."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_TESTONLY_INSTALLABLE_ERROR=ON
        -DSPKB_SHOW_M_DEPLIBS_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "ERROR: TESTONLY lib 'b_test_utils' not allowed with INSTALLABLE executable"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_8
      MESSAGE "Show error when trying to link against non-TESTONLY lib."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_NON_TESTONLY_LIB_ERROR=ON
        -DSPKB_SHOW_TESTONLY_INSTALLABLE_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        " ERROR: 'simplecxx' in TESTONLYLIBS not a TESTONLY lib"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_9
      MESSAGE "Show error when trying to link SE package lib using TESTONLYLIBS"
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_IMPORTED_LIBS_THIS_PKG_ERROR=ON
        -DSPKB_SHOW_NON_TESTONLY_LIB_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "ERROR: Lib 'pws_b' in IMPORTEDLIBS is in this SE package"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_10
      MESSAGE "Show deprecated warning when trying to link TESTONLY lib using DEBLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_TESTONLY_DEBLIBS_WARNING=ON
        -DSPKB_SHOW_IMPORTED_LIBS_THIS_PKG_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Passing TESTONLY lib 'b_mixed_lang' through DEPLIBS is deprecated"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage_LibPrefix"

    TEST_11
      MESSAGE "Show deprecated warning when trying to link non-TESTONLY lib using DEBLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_NONTESTONLY_DEBLIBS_WARNING=ON
        -DSPKB_SHOW_TESTONLY_DEBLIBS_WARNING=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Passing non-TESTONLY lib 'pws_b' through DEPLIBS is deprecated"
        "packages/with_subpackages/b/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Generating done"
        "Build files have been written to: .*/TriBITS_TribitsExampleProject_ALL_ST_LibUsage_LibPrefix"

    )


ENDIF()


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_SimpleCxx_DEBUG_int64
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
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
  TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
    MESSAGE "Run all the tests with raw 'ctest'"
    PASS_REGULAR_EXPRESSION_ALL
      "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
      "100% tests passed, 0 tests failed out of 1"
  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_WrapExternal
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0
    MESSAGE "Copy TribitsExampleProject so that we can modify it."
    CMND cp
    ARGS -r ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject .

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
      "Setting TribitsExProj_ENABLE_WithSubpackagesA=ON because WrapExternal has a required dependence on WithSubpackagesA"
      "Setting TribitsExProj_ENABLE_MixedLang=ON because WrapExternal has an optional dependence on MixedLang"
      "Setting TribitsExProj_ENABLE_SimpleCxx=ON because WithSubpackagesA has a required dependence on SimpleCxx"
      "Final set of enabled packages:  SimpleCxx MixedLang WithSubpackages WrapExternal 4"
      "Final set of enabled SE packages:  SimpleCxx MixedLang WithSubpackagesA WithSubpackages WrapExternal 5"
      "Final set of non-enabled packages:  0"
      "Final set of non-enabled SE packages:  WithSubpackagesB WithSubpackagesC 2"
      "This package has no unfiltered binary files so consider out of date"
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_WrapExternal"
  TEST_2 CMND make ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target WrapExternal_run_external_func"
  TEST_3 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
    PASS_REGULAR_EXPRESSION_ALL
      "WrapExternal_run_external_func${TEST_MPI_1_SUFFIX} .* Passed"
      "100% tests passed, 0 tests failed out of 1"

  TEST_4 CMND sleep ARGS 1s
     MESSAGE "Sleep for 1 sec for systems were time stamps are only accurate to 1 sec"
  TEST_5 CMND touch
     ARGS TribitsExampleProject/packages/with_subpackages/a/A.cpp
     MESSAGE "Test that changing upstream source will trigger rebuild"
  TEST_6 CMND ${CMAKE_COMMAND} ARGS TribitsExampleProject
      -DWrapExternal_SHOW_MOST_RECENT_FILES=TRUE
     MESSAGE "Recofigure with changed upstream source"
    PASS_REGULAR_EXPRESSION_ALL
      "Most recent file in ./packages/with_subpackages/ is ./a/A.cpp"
      "Overall most recent modified file is in ./packages/with_subpackages/ and is ./a/A.cpp"
      "The upstream SE package source file ./a/A.cpp is more recent than this package's binary file ./WrapExternal_run_external_func.exe"
      "Blowing away WrapExternal build dir external_func/ so it will build from scratch"
  TEST_7 CMND make ARGS ${CTEST_BUILD_FLAGS}
    MESSAGE "Rebuild only exteranl_func"
    PASS_REGULAR_EXPRESSION_ALL
      "Built target simplecxx"
      "Built target mixedlang"
      "Built target pws_a"
      "Generating external_func/libexternal_func.a"
      "Linking CXX executable WrapExternal_run_external_func.exe"

  TEST_8 CMND sleep ARGS 1s
     MESSAGE "Sleep for 1 sec for systems were time stamps are only accurate to 1 sec"
  TEST_9 CMND ${CMAKE_COMMAND}
     ARGS  -DSimpleCxx_ENABLE_DEBUG=ON  -DWrapExternal_SHOW_MOST_RECENT_FILES=FALSE
       TribitsExampleProject
     MESSAGE "Recofigure changing the debug mode to trigger rebuild"
    PASS_REGULAR_EXPRESSION_ALL
      "The upstream SE package binary file ./src/SimpleCxx_config.h is more recent than this package's binary file ./WrapExternal_run_external_func.exe"
      "Blowing away WrapExternal build dir external_func/ so it will build from scratch"
  TEST_10 CMND make ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Generating external_func/libexternal_func.a"
      "Linking CXX executable WrapExternal_run_external_func.exe"

  TEST_11 CMND sleep ARGS 1s
     MESSAGE "Sleep for 1 sec for systems were time stamps are only accurate to 1 sec"
  TEST_12 CMND touch
     ARGS TribitsExampleProject/packages/wrap_external/external_func/configure.py
     MESSAGE "Test that changing the external file will trigger rebuild"
  TEST_13 CMND ${CMAKE_COMMAND} ARGS TribitsExampleProject
     MESSAGE "Recofigure with changes external file"
    PASS_REGULAR_EXPRESSION_ALL
      "The this package's source file ./external_func/configure.py is more recent than this package's binary file ./WrapExternal_run_external_func.exe"
      "Blowing away WrapExternal build dir external_func/ so it will build from scratch"
  TEST_14 CMND make ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Generating external_func/libexternal_func.a"
      "Linking CXX executable WrapExternal_run_external_func.exe"

  TEST_15 CMND ${CMAKE_COMMAND} ARGS TribitsExampleProject
     MESSAGE "Recofigure with no changes that will not do anything"
    PASS_REGULAR_EXPRESSION_ALL
      "This package's most recent binary file ./WrapExternal_run_external_func.exe is more recent than its upstream SE package source or binary files or this package's source files"
  TEST_16 CMND make ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target build_external_func"
      "Built target WrapExternal_run_external_func"

  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0
    MESSSAGE "Do a verbose configure and check that the packages coupling variables are correct."
    CMND ${CMAKE_COMMAND}
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=OFF
      -DTribitsExProj_ENABLE_DEBUG=OFF
      -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=OFF
      -DTribitsExProj_ENABLE_INSTALL_CMAKE_CONFIG_FILES=OFF
      -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
      -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      -DTribitsExProj_VERBOSE_CONFIGURE=ON
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Final set of enabled packages:  SimpleCxx WithSubpackages WrapExternal 3"
      "Final set of enabled SE packages:  SimpleCxx WithSubpackagesA WithSubpackagesB WithSubpackagesC WithSubpackages WrapExternal 6"

      "HeaderOnlyTpl_INCLUDE_DIRS='.+/doc/examples/tpls/HeaderOnlyTpl'"
      "-- TPL_HeaderOnlyTpl_INCLUDE_DIRS='.+/doc/examples/tpls/HeaderOnlyTpl'"

      "SimpleCxx_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src[;].+/TribitsExampleProject/packages/simple_cxx/src[;].+/doc/examples/tpls/HeaderOnlyTpl'"
      "SimpleCxx_LIBRARY_DIRS='.+/packages/simple_cxx/src'"
      "SimpleCxx_LIBRARIES='simplecxx'"

      "WithSubpackagesA_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure[;].+/TribitsExampleProject/packages/with_subpackages/a[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src[;].+/TribitsExampleProject/packages/simple_cxx/src[;].+/doc/examples/tpls/HeaderOnlyTpl'"
      "WithSubpackagesA_LIBRARY_DIRS='.+/packages/with_subpackages/a[;].+/packages/simple_cxx/src'"
      "WithSubpackagesA_LIBRARIES='pws_a'"

      "b_test_utils_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/b/src[;].+/TribitsExampleProject/packages/with_subpackages/b/src[;].+/TribitsExampleProject/packages/with_subpackages/a[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src[;].+/TribitsExampleProject/packages/simple_cxx/src[;].+/doc/examples/tpls/HeaderOnlyTpl[;].+/TribitsExampleProject/packages/with_subpackages/b/tests/testlib'"

      "WithSubpackagesB_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/b/src[;].+/TribitsExampleProject/packages/with_subpackages/b/src[;].+/TribitsExampleProject/packages/with_subpackages/a[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src[;].+/TribitsExampleProject/packages/simple_cxx/src[;].+/doc/examples/tpls/HeaderOnlyTpl'"
      "WithSubpackagesB_LIBRARY_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/b/src;.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/a;.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src'"
      "WithSubpackagesB_LIBRARIES='pws_b'"

      "WithSubpackagesC_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/b/src[;].+/TribitsExampleProject/packages/with_subpackages/b/src[;].+/TribitsExampleProject/packages/with_subpackages/a[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src[;].+/TribitsExampleProject/packages/simple_cxx/src[;].+/doc/examples/tpls/HeaderOnlyTpl[;].+/TribitsExampleProject/packages/with_subpackages/c'"
      "WithSubpackagesC_LIBRARY_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/b/src[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/a[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/c'"
      "WithSubpackagesC_LIBRARIES='pws_c'"

      "WithSubpackages_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/b/src[;].+/TribitsExampleProject/packages/with_subpackages/b/src[;].+/TribitsExampleProject/packages/with_subpackages/a[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src[;].+/TribitsExampleProject/packages/simple_cxx/src[;].+/doc/examples/tpls/HeaderOnlyTpl[;].+/TribitsExampleProject/packages/with_subpackages/c'"
      "WithSubpackages_LIBRARY_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/b/src[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/a[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/with_subpackages/c'"
      "WithSubpackages_LIBRARIES='pws_c.pws_b.pws_a'"

      "WrapExternal_INCLUDE_DIRS='.+/packages/wrap_external/external_func'"
      "WrapExternal_LIBRARY_DIRS=''"
      "WrapExternal_LIBRARIES='external_func.pws_a'"

      "pws_b_TARGET_NAME='pws_b'"
      "b_test_TARGET_NAME='WithSubpackagesB_b_test'"
      "test_of_b_TEST_NAME='WithSubpackagesB_test_of_b'"
      "c_util_TEST_NAME='WithSubpackagesC_test_of_c_util${TEST_MPI_1_SUFFIX}'"

      "-- simplecxx:LINK_LIBS=''"
      "-- HelloWorldTests:LINK_LIBS='simplecxx'"
      "-- pws_a:LINK_LIBS='simplecxx'"
      "-- a_test:LINK_LIBS='pws_a'"
      "-- pws_b:LINK_LIBS='pws_a[;]simplecxx'"
      "-- b_test:LINK_LIBS='pws_b'"
      "-- c_util:LINK_LIBS='pws_b[;]pws_a'"
      "-- pws_c:LINK_LIBS='pws_b[;]pws_a'"
      "-- c_test:LINK_LIBS='pws_c'"
      "-- run_external_func:LINK_LIBS='external_func[;]pws_a'"

  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_WrapExternal_VerboseConfigure
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0
    MESSSAGE "Do a verbose configure and check that the packages coupling variables are correct."
    CMND ${CMAKE_COMMAND}
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=ON
      -DTribitsExProj_ENABLE_DEBUG=OFF
      -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=OFF
      -DTribitsExProj_ENABLE_INSTALL_CMAKE_CONFIG_FILES=OFF
      -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
      -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      -DTribitsExProj_VERBOSE_CONFIGURE=ON
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Final set of enabled packages:  SimpleCxx MixedLang WithSubpackages WrapExternal 4"
      "Final set of enabled SE packages:  SimpleCxx MixedLang WithSubpackagesA WithSubpackagesB WithSubpackagesC WithSubpackages WrapExternal 7"

      # ToDo Fill in more checks for libs, variables, etc. 

      "-- pws_b_TARGET_NAME='pws_b'"
      "-- b_test_TARGET_NAME='WithSubpackagesB_b_test'"
      "-- test_of_b_TEST_NAME='WithSubpackagesB_test_of_b'"
      "-- c_util_TEST_NAME='WithSubpackagesC_test_of_c_util${TEST_MPI_1_SUFFIX}'"
      "-- test_of_c_b_mixed_lang_TARGET_NAME='WithSubpackagesC_test_of_c_b_mixed_lang'"
      "-- test_of_c_b_mixed_lang_TEST_NAME='WithSubpackagesC_test_of_c_b_mixed_lang${TEST_MPI_1_SUFFIX}'"

      "-- simplecxx:LINK_LIBS=''"
      "-- HelloWorldTests:LINK_LIBS='simplecxx'"
      "-- mixedlang:LINK_LIBS=''"
      "-- RayTracerTests:LINK_LIBS='mixedlang'"
      "-- pws_a:LINK_LIBS='simplecxx'"
      "-- a_test:LINK_LIBS='pws_a'"
      "-- pws_b:LINK_LIBS='pws_a[;]simplecxx'"
      "-- b_test:LINK_LIBS='pws_b[;]mixedlang'"
      "-- b_mixed_lang:LINK_LIBS='pws_b[;]mixedlang'"
      "-- test_of_b_mixed_lang:LINK_LIBS='b_test_utils[;]b_mixed_lang[;]pws_b[;]mixedlang'"
      "-- c_util:LINK_LIBS='pws_b[;]pws_a'"
      "-- pws_c:LINK_LIBS='pws_b[;]pws_a'"
      "-- c_test:LINK_LIBS='pws_c'"
      "-- c_b_mixed_lang:LINK_LIBS='pws_c[;]b_mixed_lang'"
      "-- test_of_c_b_mixed_lang:LINK_LIBS='c_b_mixed_lang[;]pws_c'"
      "-- run_external_func:LINK_LIBS='external_func[;]pws_a'"

  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0
    MESSAGE "Copy TribitsExampleProject so that we can modify it."
    CMND cp
    ARGS -r ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject .

  TEST_1
    MESSAGE "Copy the with_subpackages package dir to base dir."
    CMND cp
    ARGS -r TribitsExampleProject/packages/with_subpackages/
      TribitsExampleProject/.

  TEST_2
    MESSAGE "Remove the packages/with_subpackages package dir."
    CMND rm
    ARGS -r TribitsExampleProject/packages/with_subpackages/

  TEST_3
    MESSSAGE "Override with_packages/ source dir and configure"
    CMND ${CMAKE_COMMAND}
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=OFF
      -DTribitsExProj_ENABLE_DEBUG=OFF
      -DTribitsExProj_ENABLE_EXPORT_MAKEFILES=OFF
      -DTribitsExProj_ENABLE_INSTALL_CMAKE_CONFIG_FILES=OFF
      -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
      -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      -DTribitsExProj_VERBOSE_CONFIGURE=ON
      -DWithSubpackages_SOURCE_DIR_OVERRIDE=with_subpackages
      TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Final set of enabled packages:  SimpleCxx WithSubpackages WrapExternal 3"
      "Final set of enabled SE packages:  SimpleCxx WithSubpackagesA WithSubpackagesB WithSubpackagesC WithSubpackages WrapExternal 6"
      "-- NOTE: WithSubpackages_SOURCE_DIR_OVERRIDE='with_subpackages' is overriding default path 'packages/with_subpackages'"
      "-- TribitsExProj_PACKAGE_DIRS='packages/simple_cxx[;]packages/mixed_lang[;]with_subpackages[;]packages/wrap_external'"
      "-- File Trace: PACKAGE    INCLUDE    .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/cmake/Dependencies[.]cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/a/cmake/Dependencies[.]cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/b/cmake/Dependencies[.]cmake"
      "-- File Trace: PACKAGE    INCLUDE    .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/c/cmake/Dependencies[.]cmake"
      "-- WithSubpackages_BINARY_DIR='.*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages'"
      "-- WithSubpackages_SOURCE_DIR='.*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages'"
      "-- WithSubpackages_BINARY_DIR='.*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages'"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/CMakeLists[.]txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/a/CMakeLists[.]txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/a/tests/CMakeLists[.]txt"

      "-- WithSubpackagesA_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/a[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/packages/simple_cxx/src[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/packages/simple_cxx/src[;].+/doc/examples/tpls/HeaderOnlyTpl'"

      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/b/CMakeLists[.]txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/b/tests/CMakeLists[.]txt"

      "-- WithSubpackagesB_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/b/src[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/b/src[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/a[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/packages/simple_cxx/src[;].+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/packages/simple_cxx/src[;].+/doc/examples/tpls/HeaderOnlyTpl'"
      "-- WithSubpackagesB_LIBRARY_DIRS='.*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/b/src[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/a[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/packages/simple_cxx/src'"

      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/c/CMakeLists[.]txt"
      "-- File Trace: PACKAGE    ADD_SUBDIR .*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/c/tests/CMakeLists[.]txt"

      "-- WithSubpackagesC_INCLUDE_DIRS='.*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/b/src[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/b/src[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/a[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/packages/simple_cxx/src[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/packages/simple_cxx/src[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/c'"
      "-- WithSubpackagesC_LIBRARY_DIRS='.*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/b/src[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/a[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/packages/simple_cxx/src[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/c'"

     "-- WithSubpackages_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir;.+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/b/src;.+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/b/src;.+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/a;.+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/packages/simple_cxx/src;.+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/packages/simple_cxx/src;.+/doc/examples/tpls/HeaderOnlyTpl;.+/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/TribitsExampleProject/with_subpackages/c'"
      "-- WithSubpackages_LIBRARY_DIRS='.*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/b/src[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/a[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/packages/simple_cxx/src[;].*/TriBITS_TribitsExampleProject_ALL_NoFortran_OverridePackageSourceDir/with_subpackages/c'"

  TEST_4 CMND make
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target WithSubpackagesC_c_test"

  TEST_5 CMND ${CMAKE_CTEST_COMMAND}
    PASS_REGULAR_EXPRESSION_ALL
      "100% tests passed, 0 tests failed out of 6"

  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_HeaderOnlyTpl_Fail
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0
    MESSAGE "Do the initial configure with a bad TPL find path"
    CMND ${CMAKE_COMMAND}
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=OFF
      -DTribitsExProj_ENABLE_SimpleCxx=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      #-DTribitsExProj_VERBOSE_CONFIGURE=ON
      -DHeaderOnlyTpl_INCLUDE_DIRS=/path_does_not_exist
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Processing enabled TPL: HeaderOnlyTpl .enabled by SimpleCxx, disable with -DTPL_ENABLE_HeaderOnlyTpl=OFF."
      "-- Searching for headers in HeaderOnlyTpl_INCLUDE_DIRS='/path_does_not_exist'"
      "-- TIP: If the TPL 'HeaderOnlyTpl' is on your system then you can set:"
      "-DHeaderOnlyTpl_INCLUDE_DIRS='<dir0>[;]<dir1>[;]...'"
      "-DTPL_HeaderOnlyTpl_INCLUDE_DIRS='<dir0>[;]<dir1>[;]...'"
      "-- ERROR: Failed finding all of the parts of TPL 'HeaderOnlyTpl' .see above., Aborting!"
      "TIP: One way to get past the configure failure for the"
      "TPL 'HeaderOnlyTpl' is to simply disable it with:"
      "  -DTPL_ENABLE_HeaderOnlyTpl=OFF"
      "which will disable it and will recursively disable all of the"
      "downstream packages that have required dependencies on it, including"
      "the package 'SimpleCxx' which triggered its enable."
      "When you reconfigure, just grep the cmake stdout for 'HeaderOnlyTpl'"
      "and then follow the disables that occur as a result to see what impact"
      "this TPL disable has on the configuration of TribitsExProj."
      "CMake Error at .+/TribitsProcessEnabledTpl[.]cmake:[0-9]+ .MESSAGE.:"
      "  ERROR: TPL_HeaderOnlyTpl_NOT_FOUND=TRUE, aborting!"
      "Call Stack .most recent call first.:"
      "-- Configuring incomplete, errors occurred!"
  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_HeaderOnlyTpl_HardEnable_Fail
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0
    MESSAGE "Do the initial configure with a bad TPL find path"
    CMND ${CMAKE_COMMAND}
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=OFF
      -DTribitsExProj_ENABLE_SimpleCxx=ON
      -DTribitsExProj_ENABLE_TESTS=ON
      #-DTribitsExProj_VERBOSE_CONFIGURE=ON
      -DTPL_ENABLE_HeaderOnlyTpl=ON
      -DHeaderOnlyTpl_INCLUDE_DIRS=/path_does_not_exist
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Processing enabled TPL: HeaderOnlyTpl .enabled explicitly, disable with -DTPL_ENABLE_HeaderOnlyTpl=OFF."
      "TIP: Even though the TPL 'HeaderOnlyTpl' was explicitly enabled in input,"
      "it can be disabled with:"
      "  -DTPL_ENABLE_HeaderOnlyTpl=OFF"
      "which will disable it and will recursively disable all of the"
      "downstream packages that have required dependencies on it."
      "When you reconfigure, just grep the cmake stdout for 'HeaderOnlyTpl'"
      "and then follow the disables that occur as a result to see what impact" 
      "this TPL disable has on the configuration of TribitsExProj."
      "-- ERROR: Failed finding all of the parts of TPL 'HeaderOnlyTpl' .see above., Aborting!"
      "CMake Error at .+/TribitsProcessEnabledTpl[.]cmake:[0-9]+ .MESSAGE.:"
      "  ERROR: TPL_HeaderOnlyTpl_NOT_FOUND=TRUE, aborting!"
      "Call Stack .most recent call first.:"
      "-- Configuring incomplete, errors occurred!"
  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ExternalPkg
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0
    MESSAGE "Copy TribitsExampleProject so that we can copy in ExteranlPkg."
    CMND cp
    ARGS -r ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject .

  TEST_1
    MESSAGE "Configure to get the compiler options"
    CMND ${CMAKE_COMMAND}
    ARGS
      ${COMMON_ENV_ARGS_PASSTHROUGH}
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=OFF
      -DTribitsExProj_ENABLE_DEBUG=OFF
      TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Configuring done"
      "Generating done"

  TEST_2
    MESSAGE "Configure asserting existence of missing ExternalPkg"
    CMND ${CMAKE_COMMAND}
    ARGS
      -DExternalPkg_ALLOW_MISSING_EXTERNAL_PACKAGE=FALSE
      .
    PASS_REGULAR_EXPRESSION_ALL
      "Error, the package ExternalPkg directory .+/TribitsExampleProject/ExternalPkg does not exist!"
      "CMake Error at .+/TribitsProcessPackagesAndDirsLists.cmake:[0-9]+ .MESSAGE.:"
      "Configuring incomplete, errors occurred!"

  TEST_3
    MESSAGE "Copy ExteranlPkg to base dir."
    CMND cp
    ARGS -r ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/ExternalPkg
      TribitsExampleProject/.


  TEST_4
    MESSAGE "Configure all packages"
    CMND ${CMAKE_COMMAND}
    ARGS
      -DTribitsExProj_ENABLE_ALL_PACKAGES=ON 
      -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON
      -DTribitsExProj_ENABLE_TESTS=ON
     .
    PASS_REGULAR_EXPRESSION_ALL
      "Setting TribitsExProj_ENABLE_ExternalPkg=ON"
      "Final set of enabled packages:  SimpleCxx ExternalPkg .+"
      "Final set of enabled SE packages:  SimpleCxx ExternalPkg .+"
      "Configuring done"
      "Generating done"

  TEST_5 CMND make
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target externalpkg"
      "Linking CXX executable ExternalPkg_test.exe"

  TEST_6 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
    PASS_REGULAR_EXPRESSION_ALL
      "ExternalPkg_test${TEST_MPI_1_SUFFIX} [.]+   Passed"
      "100% tests passed, 0 tests failed out of"

  # ToDo: Add usage of ExternalPkg in downstream package.

  )


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_DisableWithSubpackagesB_EnableWithSubpackagesB
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0
    MESSAGE "Just do dependency analysis to test enabling of parent package"
     " with eanbled subpackages even if is disabled"
    CMND ${CMAKE_COMMAND}
    ARGS
      -DTribitsExProj_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
      -DTribitsExProj_ENABLE_Fortran=OFF
      -DTribitsExProj_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING=ON
      -DTribitsExProj_ENABLE_WithSubpackagesA=OFF
      -DTribitsExProj_ENABLE_WithSubpackagesB=ON
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject
    PASS_REGULAR_EXPRESSION_ALL
      "Explicitly enabled SE packages on input .by user.:  WithSubpackagesB 1"
      "-- Setting TribitsExProj_ENABLE_WrapExternal=OFF because WrapExternal has a required library dependence on disabled package WithSubpackagesA"
      "Enabling all parent packages that have at least one subpackage enabled [.][.][.]"
      "-- Setting TribitsExProj_ENABLE_WithSubpackages=ON because TribitsExProj_ENABLE_WithSubpackagesB=ON"
      "Final set of enabled packages:  SimpleCxx WithSubpackages 2"
      "Final set of enabled SE packages:  SimpleCxx WithSubpackagesB WithSubpackages 3"

  )
# NOTE: The above test is the *only* test that we have that checks that a
# parent package is re-enabled if any of its subpackages are enabled!