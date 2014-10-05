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
ELSE()
  SET(TPL_MPI_FILE_TRACE "")
ENDIF()


TRIBITS_ADD_ADVANCED_TEST( TribitsExampleProject_ALL_ST_NoFortran
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1

  TEST_0 CMND ${CMAKE_COMMAND}
    MESSAGE "Do the initial configure (and test a lot of things at once"
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
      "Warning: Setting TribitsExProj_ENABLE_MixedLanguage=OFF because "
      "Printing package dependencies ..."
      "-- TribitsExProj_PACKAGES: SimpleCxx MixedLanguage PackageWithSubpackages WrapExternal"
      "-- TribitsExProj_SE_PACKAGES: SimpleCxx MixedLanguage PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages WrapExternal"

      "-- SimpleCxx: No dependencies!"
      "-- MixedLanguage: No dependencies!"
      "-- PackageWithSubpackagesSubpackageA_LIB_REQUIRED_DEP_PACKAGES: SimpleCxx"
      "-- PackageWithSubpackagesSubpackageB_LIB_REQUIRED_DEP_PACKAGES: SimpleCxx"
      "-- PackageWithSubpackagesSubpackageB_LIB_OPTIONAL_DEP_PACKAGES: PackageWithSubpackagesSubpackageA"
      "-- PackageWithSubpackagesSubpackageB_TEST_OPTIONAL_DEP_PACKAGES: MixedLanguage"
      "-- PackageWithSubpackagesSubpackageC_LIB_REQUIRED_DEP_PACKAGES: PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB"
      "-- PackageWithSubpackages_LIB_REQUIRED_DEP_PACKAGES: PackageWithSubpackagesSubpackageA"
      "-- PackageWithSubpackages_LIB_OPTIONAL_DEP_PACKAGES: PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC"
      "-- WrapExternal_LIB_REQUIRED_DEP_PACKAGES: PackageWithSubpackagesSubpackageA"
      "-- WrapExternal_LIB_OPTIONAL_DEP_PACKAGES: MixedLanguage"
      "-- SimpleCxx: No library dependencies!"
      "-- PackageWithSubpackagesSubpackageA_FULL_ENABLED_DEP_PACKAGES: SimpleCxx"
      "-- PackageWithSubpackagesSubpackageB_FULL_ENABLED_DEP_PACKAGES: PackageWithSubpackagesSubpackageA SimpleCxx"
      "-- PackageWithSubpackagesSubpackageC_FULL_ENABLED_DEP_PACKAGES: PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageA SimpleCxx"
      "-- PackageWithSubpackages_FULL_ENABLED_DEP_PACKAGES: PackageWithSubpackagesSubpackageC PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageA SimpleCxx"
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
      "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_ST_NoFortran"
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
      "PackageWithSubpackagesSubpackageA_test_of_a .* Passed"
      "PackageWithSubpackagesSubpackageB_test_of_b .* Passed"
      "PackageWithSubpackagesSubpackageC_test_of_c .* Passed"
      "PackageWithSubpackagesSubpackageC_test_of_c_util.* Passed" 
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
      "Installing: .+/install/lib/cmake/PackageWithSubpackages/PackageWithSubpackagesConfig.cmake"
      "Installing: .+/install/include/Makefile.export.PackageWithSubpackages"
      "Installing: .+/install/lib/libpws_a.a"
      "Installing: .+/install/include/A.hpp"
      "Installing: .+/install/lib/cmake/PackageWithSubpackagesSubpackageA/PackageWithSubpackagesSubpackageAConfig.cmake"
      "Installing: .+/install/lib/cmake/PackageWithSubpackagesSubpackageA/PackageWithSubpackagesSubpackageATargets.cmake"
      "Installing: .+/install/lib/cmake/PackageWithSubpackagesSubpackageA/PackageWithSubpackagesSubpackageATargets-release.cmake"
      "Installing: .+/install/include/Makefile.export.PackageWithSubpackagesSubpackageA"
      "Installing: .+/install/lib/libpws_b.a"
      "Installing: .+/install/include/B.hpp"
      "Installing: .+/install/lib/cmake/PackageWithSubpackagesSubpackageB/PackageWithSubpackagesSubpackageBConfig.cmake"
      "Installing: .+/install/lib/cmake/PackageWithSubpackagesSubpackageB/PackageWithSubpackagesSubpackageBTargets.cmake"
      "Installing: .+/install/lib/cmake/PackageWithSubpackagesSubpackageB/PackageWithSubpackagesSubpackageBTargets-release.cmake"
      "Installing: .+/install/include/Makefile.export.PackageWithSubpackagesSubpackageB"
      "Installing: .+/install/lib/libpws_c.a"
      "Installing: .+/install/include/C.hpp"
      "Installing: .+/install/lib/cmake/PackageWithSubpackagesSubpackageC/PackageWithSubpackagesSubpackageCConfig.cmake"
      "Installing: .+/install/lib/cmake/PackageWithSubpackagesSubpackageC/PackageWithSubpackagesSubpackageCTargets.cmake"
      "Installing: .+/install/lib/cmake/PackageWithSubpackagesSubpackageC/PackageWithSubpackagesSubpackageCTargets-release.cmake"
      "Installing: .+/install/include/Makefile.export.PackageWithSubpackagesSubpackageC"

  TEST_4 CMND ${CMAKE_COMMAND}
    ARGS
      -DDUMMY_PROJECT_NAME=DummyProject
      -DDUMMY_PROJECT_DIR=dummy_client_of_PackageWithSubpackages
      -DEXPORT_VAR_PREFIX=PackageWithSubpackages
      -DEXPORT_CONFIG_FILE=../install/lib/cmake/PackageWithSubpackages/PackageWithSubpackagesConfig.cmake
      -DCMAKE_COMMAND=${CMAKE_COMMAND}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/RunDummyPackageClientBulid.cmake
    MESSAGE "Create and configure a dummy project that includes PackageWithSubpackagesConfig.cmake"
    PASS_REGULAR_EXPRESSION_ALL
      "PackageWithSubpackages_CXX_FLAGS = '.*'"
      "PackageWithSubpackages_C_FLAGS = '.*'"
      "PackageWithSubpackages_FORTRAN_FLAGS = '.*'"
      "PackageWithSubpackages_EXTRA_LD_FLAGS = '.*'"
      "PackageWithSubpackages_SHARED_LIB_RPATH_COMMAND = '.*'"
      "PackageWithSubpackages_BUILD_SHARED_LIBS = '.*'"
      "PackageWithSubpackages_LINKER = '.+'"
      "PackageWithSubpackages_AR = '.+'"
      "PackageWithSubpackages_INCLUDE_DIRS = '.+/install/lib/cmake/PackageWithSubpackages/../../../include'"
      "PackageWithSubpackages_LIBRARY_DIRS = '.+/install/lib/cmake/PackageWithSubpackages/../../../lib'"
      "PackageWithSubpackages_LIBRARIES = 'pws_c.pws_b.pws_a.simplecxx'"
      "PackageWithSubpackages_TPL_INCLUDE_DIRS = ''"
      "PackageWithSubpackages_TPL_LIBRARY_DIRS = ''"
      "PackageWithSubpackages_TPL_LIBRARIES = ''"
      "PackageWithSubpackages_MPI_LIBRARIES = ''"
      "PackageWithSubpackages_MPI_LIBRARY_DIRS = ''"
      "PackageWithSubpackages_MPI_INCLUDE_DIRS = ''"
      "PackageWithSubpackages_MPI_EXEC = '${MPI_EXEC}'"
      "PackageWithSubpackages_MPI_EXEC_MAX_NUMPROCS = '${MPI_EXEC_MAX_NUMPROCS}'"
      "PackageWithSubpackages_MPI_EXEC_NUMPROCS_FLAG = '${MPI_EXEC_NUMPROCS_FLAG}'"
      "PackageWithSubpackages_PACKAGE_LIST = 'PackageWithSubpackagesSubpackageC.PackageWithSubpackagesSubpackageB.PackageWithSubpackagesSubpackageA.SimpleCxx'"
      "PackageWithSubpackages_TPL_LIST = ''"

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
      "Only in .*/TribitsExampleProject: .gitignore"
      "Only in .*/TribitsExampleProject/packages: mixed_language"
      "Only in .*/TribitsExampleProject/packages: wrap_external"
  )


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
      "-- PackageWithSubpackagesSubpackageA_FULL_ENABLED_DEP_PACKAGES: SimpleCxx"
      "-- PackageWithSubpackages_FULL_ENABLED_DEP_PACKAGES: PackageWithSubpackagesSubpackageA SimpleCxx"
      "Explicitly enabled packages on input .by user.:  0"
      "Explicitly disabled packages on input .by user or by default.:  MixedLanguage WrapExternal 2"
      "Final set of enabled packages:  SimpleCxx PackageWithSubpackages 2"
      "Final set of enabled SE packages:  SimpleCxx PackageWithSubpackagesSubpackageA PackageWithSubpackages 3"
      "Final set of non-enabled packages:  MixedLanguage WrapExternal 2"
      "Final set of non-enabled SE packages:  MixedLanguage PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC WrapExternal 4"

  TEST_1 CMND make ARGS ${CTEST_BUILD_FLAGS}
    MESSAGE "Build the default 'all' target using raw 'make'"
    PASS_REGULAR_EXPRESSION_ALL
      "Built target simplecxx"
      "Built target pws_a"

  TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
    MESSAGE "Run all the tests with raw 'ctest'"
    PASS_REGULAR_EXPRESSION_ALL
      "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
      "PackageWithSubpackagesSubpackageA_test_of_a .* Passed"
      "100% tests passed, 0 tests failed out of 2"

  TEST_3 CMND make ARGS install ${CTEST_BUILD_FLAGS}
    MESSAGE "Build 'install' target using raw 'make'"
    PASS_REGULAR_EXPRESSION_ALL
      "Installing: .+/install/lib/cmake/PackageWithSubpackages/PackageWithSubpackagesConfig.cmake"
      "Installing: .+/install/include/Makefile.export.PackageWithSubpackages"

  TEST_4 CMND ${CMAKE_COMMAND}
    ARGS
      -DDUMMY_PROJECT_NAME=DummyProject
      -DDUMMY_PROJECT_DIR=dummy_client_of_PackageWithSubpackages
      -DEXPORT_VAR_PREFIX=PackageWithSubpackages
      -DEXPORT_CONFIG_FILE=../install/lib/cmake/PackageWithSubpackages/PackageWithSubpackagesConfig.cmake
      -DCMAKE_COMMAND=${CMAKE_COMMAND}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/RunDummyPackageClientBulid.cmake
    MESSAGE "Create and configure a dummy project that includes PackageWithSubpackagesConfig.cmake"
    PASS_REGULAR_EXPRESSION_ALL
      "PackageWithSubpackages_INCLUDE_DIRS = '.+/install/lib/cmake/PackageWithSubpackages/../../../include'"
      "PackageWithSubpackages_LIBRARY_DIRS = '.+/install/lib/cmake/PackageWithSubpackages/../../../lib'"
      "PackageWithSubpackages_LIBRARIES = 'pws_a.simplecxx'"
      "PackageWithSubpackages_PACKAGE_LIST = 'PackageWithSubpackagesSubpackageA.SimpleCxx'"
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
        "Setting TribitsExProj_ENABLE_MixedLanguage=ON"
        "Setting TribitsExProj_ENABLE_PackageWithSubpackages=ON"
        "Final set of enabled packages:  SimpleCxx MixedLanguage PackageWithSubpackages 3"
        "Final set of enabled SE packages:  SimpleCxx MixedLanguage PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages 6"
        "Final set of non-enabled packages:  WrapExternal 1"
        "Configuring done"
        "Generating done"
        "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_ST"
    TEST_1 CMND make
      MESSAGE "Build the default 'all' target using raw 'make'"
      ARGS ${CTEST_BUILD_FLAGS}
      PASS_REGULAR_EXPRESSION_ALL
        "Built target simplecxx"
        "Built target mixedlanguage"
        "Built target pws_a"
        "Built target pws_b"
        "Built target pws_c"
    TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
      MESSAGE "Run all the tests with raw 'ctest'"
      PASS_REGULAR_EXPRESSION_ALL
        "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
        "MixedLanguage_RayTracerTests${TEST_MPI_1_SUFFIX} .* Passed"
        "PackageWithSubpackagesSubpackageA_test_of_a .* Passed"
        "PackageWithSubpackagesSubpackageB_test_of_b .* Passed"
        "PackageWithSubpackagesSubpackageB_test_of_b_mixed_lang.* Passed"
        "PackageWithSubpackagesSubpackageC_test_of_c_util.* Passed" 
        "PackageWithSubpackagesSubpackageC_test_of_c .* Passed"
        "PackageWithSubpackagesSubpackageC_test_of_c_b_mixed_lang.* Passed"
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
        "Setting TribitsExProj_ENABLE_MixedLanguage=ON"
        "Setting TribitsExProj_ENABLE_PackageWithSubpackages=ON"
        "Final set of enabled packages:  SimpleCxx MixedLanguage PackageWithSubpackages 3"
        "Final set of enabled SE packages:  SimpleCxx MixedLanguage PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages 6"
        "Final set of non-enabled packages:  WrapExternal 1"
        "Configuring done"
        "Generating done"
        "Build files have been written to: .*doc/examples/UnitTests/TriBITS_TribitsExampleProject_ALL_ST_LibPrefix"
    TEST_1 CMND make
      MESSAGE "Build the default 'all' target using raw 'make'"
      ARGS ${CTEST_BUILD_FLAGS}
      PASS_REGULAR_EXPRESSION_ALL
        "Built target tep_simplecxx"
        "Built target tep_mixedlanguage"
        "Built target tep_pws_a"
        "Built target tep_pws_b"
        "Built target tep_pws_c"
    TEST_2 CMND ${CMAKE_CTEST_COMMAND} ARGS -VV
      MESSAGE "Run all the tests with raw 'ctest'"
      PASS_REGULAR_EXPRESSION_ALL
        "SimpleCxx_HelloWorldTests${TEST_MPI_1_SUFFIX} .* Passed"
        "MixedLanguage_RayTracerTests${TEST_MPI_1_SUFFIX} .* Passed"
        "PackageWithSubpackagesSubpackageA_test_of_a .* Passed"
        "PackageWithSubpackagesSubpackageB_test_of_b .* Passed"
        "PackageWithSubpackagesSubpackageB_test_of_b_mixed_lang.* Passed"
        "PackageWithSubpackagesSubpackageC_test_of_c_util.* Passed" 
        "PackageWithSubpackagesSubpackageC_test_of_c .* Passed"
        "PackageWithSubpackagesSubpackageC_test_of_c_b_mixed_lang.* Passed"
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
        "Final set of enabled packages:  SimpleCxx MixedLanguage PackageWithSubpackages 3"
        "Final set of enabled SE packages:  SimpleCxx MixedLanguage PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages 6"
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
        "packages/package_with_subpackages/B/src/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "packages/package_with_subpackages/B/tests/testlib/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "TribitsExampleProject/packages/package_with_subpackages/B/cmake/Dependencies.cmake"
        "packages/package_with_subpackages/B/tests/testlib/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "TribitsExampleProject/packages/package_with_subpackages/C/cmake/Dependencies.cmake"
        "packages/package_with_subpackages/C/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "packages/package_with_subpackages/C/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "packages/package_with_subpackages/B/src/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "ERROR: TESTONLY lib 'b_mixed_lang' not allowed with INSTALLABLE executable"
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_8
      MESSAGE "Show error when trying to link against non-TESTONLY lib."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_NON_TESTONLY_LIB_ERROR=ON
        -DSPKB_SHOW_TESTONLY_INSTALLABLE_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        " ERROR: 'simplecxx' in TESTONLYLIBS not a TESTONLY lib"
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_9
      MESSAGE "Show error when trying to link SE package lib using TESTONLYLIBS"
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_IMPORTED_LIBS_THIS_PKG_ERROR=ON
        -DSPKB_SHOW_NON_TESTONLY_LIB_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "ERROR: Lib 'pws_b' in IMPORTEDLIBS is in this SE package"
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_10
      MESSAGE "Show deprecated warning when trying to link TESTONLY lib using DEBLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_TESTONLY_DEBLIBS_WARNING=ON
        -DSPKB_SHOW_IMPORTED_LIBS_THIS_PKG_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Passing TESTONLY lib 'b_mixed_lang' through DEPLIBS is deprecated"
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
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
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
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
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
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
        "Final set of enabled packages:  SimpleCxx MixedLanguage PackageWithSubpackages 3"
        "Final set of enabled SE packages:  SimpleCxx MixedLanguage PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages 6"
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
        "packages/package_with_subpackages/B/src/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "packages/package_with_subpackages/B/tests/testlib/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "TribitsExampleProject/packages/package_with_subpackages/B/cmake/Dependencies.cmake"
        "packages/package_with_subpackages/B/tests/testlib/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "TribitsExampleProject/packages/package_with_subpackages/C/cmake/Dependencies.cmake"
        "packages/package_with_subpackages/C/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "packages/package_with_subpackages/C/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "packages/package_with_subpackages/B/src/CMakeLists.txt:.* .TRIBITS_ADD_LIBRARY."
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
        "ERROR: TESTONLY lib 'b_mixed_lang' not allowed with INSTALLABLE executable"
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_8
      MESSAGE "Show error when trying to link against non-TESTONLY lib."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_NON_TESTONLY_LIB_ERROR=ON
        -DSPKB_SHOW_TESTONLY_INSTALLABLE_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        " ERROR: 'simplecxx' in TESTONLYLIBS not a TESTONLY lib"
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_9
      MESSAGE "Show error when trying to link SE package lib using TESTONLYLIBS"
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_IMPORTED_LIBS_THIS_PKG_ERROR=ON
        -DSPKB_SHOW_NON_TESTONLY_LIB_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "ERROR: Lib 'pws_b' in IMPORTEDLIBS is in this SE package"
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
        "Configuring incomplete, errors occurred!"

    TEST_10
      MESSAGE "Show deprecated warning when trying to link TESTONLY lib using DEBLIBS."
      CMND ${CMAKE_COMMAND}
      ARGS -DSPKB_SHOW_TESTONLY_DEBLIBS_WARNING=ON
        -DSPKB_SHOW_IMPORTED_LIBS_THIS_PKG_ERROR=
        .
      PASS_REGULAR_EXPRESSION_ALL
        "WARNING: Passing TESTONLY lib 'b_mixed_lang' through DEPLIBS is deprecated"
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
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
        "packages/package_with_subpackages/B/tests/CMakeLists.txt:.* .TRIBITS_ADD_EXECUTABLE_AND_TEST."
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

  TEST_0 CMND cp
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
      "Setting TribitsExProj_ENABLE_PackageWithSubpackagesSubpackageA=ON because WrapExternal has a required dependence on PackageWithSubpackagesSubpackageA"
      "Setting TribitsExProj_ENABLE_MixedLanguage=ON because WrapExternal has an optional dependence on MixedLanguage"
      "Setting TribitsExProj_ENABLE_SimpleCxx=ON because PackageWithSubpackagesSubpackageA has a required dependence on SimpleCxx"
      "Final set of enabled packages:  SimpleCxx MixedLanguage PackageWithSubpackages WrapExternal 4"
      "Final set of enabled SE packages:  SimpleCxx MixedLanguage PackageWithSubpackagesSubpackageA PackageWithSubpackages WrapExternal 5"
      "Final set of non-enabled packages:  0"
      "Final set of non-enabled SE packages:  PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC 2"
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
     ARGS TribitsExampleProject/packages/package_with_subpackages/A/A.cpp
     MESSAGE "Test that changing upstream source will trigger rebuild"
  TEST_6 CMND ${CMAKE_COMMAND} ARGS TribitsExampleProject
      -DWrapExternal_SHOW_MOST_RECENT_FILES=TRUE
     MESSAGE "Recofigure with changed upstream source"
    PASS_REGULAR_EXPRESSION_ALL
      "Most recent file in ./packages/package_with_subpackages/ is ./A/A.cpp"
      "Overall most recent modified file is in ./packages/package_with_subpackages/ and is ./A/A.cpp"
      "The upstream SE package source file ./A/A.cpp is more recent than this package's binary file ./WrapExternal_run_external_func.exe"
      "Blowing away WrapExternal build dir external_func/ so it will build from scratch"
  TEST_7 CMND make ARGS ${CTEST_BUILD_FLAGS}
    MESSAGE "Rebuild only exteranl_func"
    PASS_REGULAR_EXPRESSION_ALL
      "Built target simplecxx"
      "Built target mixedlanguage"
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
  CATEGORIES CONTINUOUS

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
      "Final set of enabled packages:  SimpleCxx PackageWithSubpackages WrapExternal 3"
      "Final set of enabled SE packages:  SimpleCxx PackageWithSubpackagesSubpackageA PackageWithSubpackagesSubpackageB PackageWithSubpackagesSubpackageC PackageWithSubpackages WrapExternal 6"

      "SimpleCxx_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure..+/packages/simple_cxx/src..+/TribitsExampleProject/packages/simple_cxx/src'"
      "SimpleCxx_LIBRARY_DIRS='.+/packages/simple_cxx/src'"
      "SimpleCxx_LIBRARIES='simplecxx'"

      "PackageWithSubpackagesSubpackageA_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure..+/TribitsExampleProject/packages/package_with_subpackages/A..+/packages/simple_cxx/src..+/TribitsExampleProject/packages/simple_cxx/src'"
      "PackageWithSubpackagesSubpackageA_LIBRARY_DIRS='.+/packages/package_with_subpackages/A..+/packages/simple_cxx/src'"
      "PackageWithSubpackagesSubpackageA_LIBRARIES='pws_a'"

      "PackageWithSubpackagesSubpackageB_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure..+/TribitsExampleProject/packages/package_with_subpackages/B..+/TribitsExampleProject/packages/package_with_subpackages/A..+/packages/simple_cxx/src..+/TribitsExampleProject/packages/simple_cxx/src'"
      "PackageWithSubpackagesSubpackageB_LIBRARY_DIRS='.+/packages/package_with_subpackages/B..+/packages/package_with_subpackages/A..+/packages/simple_cxx/src'"
      "PackageWithSubpackagesSubpackageB_LIBRARIES='pws_b'"

      "PackageWithSubpackagesSubpackageC_INCLUDE_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/package_with_subpackages/B/src..+/TribitsExampleProject/packages/package_with_subpackages/B/src..+/TribitsExampleProject/packages/package_with_subpackages/A..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src..+/TribitsExampleProject/packages/simple_cxx/src..+/TribitsExampleProject/packages/package_with_subpackages/C'"
      "PackageWithSubpackagesSubpackageC_LIBRARY_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/package_with_subpackages/B/src..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/package_with_subpackages/A..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/package_with_subpackages/C'"
      "PackageWithSubpackagesSubpackageC_LIBRARIES='pws_c'"

      "PackageWithSubpackages_INCLUDE_DIRS='.*/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/package_with_subpackages/B/src..+/TribitsExampleProject/packages/package_with_subpackages/B/src..+/TribitsExampleProject/packages/package_with_subpackages/A..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src..+/TribitsExampleProject/packages/simple_cxx/src..+/TribitsExampleProject/packages/package_with_subpackages/C'"
      "PackageWithSubpackages_LIBRARY_DIRS='.+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/package_with_subpackages/B/src..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/package_with_subpackages/A..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/simple_cxx/src..+/TriBITS_TribitsExampleProject_ALL_NoFortran_WrapExternal_VerboseConfigure/packages/package_with_subpackages/C'"
      "PackageWithSubpackages_LIBRARIES='pws_c.pws_b.pws_a'"

      "WrapExternal_INCLUDE_DIRS='.+/packages/wrap_external/external_func'"
      "WrapExternal_LIBRARY_DIRS=''"
      "WrapExternal_LIBRARIES='external_func.pws_a'"

  )
