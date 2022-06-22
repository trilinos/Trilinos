IF(COMMAND TRIBITS_PACKAGE_DECL)
  SET(KOKKOSKERNELS_HAS_TRILINOS ON CACHE BOOL "")
ELSE()
  SET(KOKKOSKERNELS_HAS_TRILINOS OFF CACHE BOOL "")
ENDIF()

GET_DIRECTORY_PROPERTY(KOKKOSKERNELS_HAS_PARENT PARENT_DIRECTORY)

SET(KOKKOSKERNELS_TOP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(KOKKOSKERNELS_TOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

IF(NOT KOKKOSKERNELS_HAS_TRILINOS)
  cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
  IF (Spack_WORKAROUND)
    #if we are explicitly using Spack for development,
    #nuke the Spack compiler
    SET(SPACK_CXX $ENV{SPACK_CXX})
    IF(SPACK_CXX)
      SET(CMAKE_CXX_COMPILER ${SPACK_CXX} CACHE STRING "the C++ compiler" FORCE)
      SET(ENV{CXX} ${SPACK_CXX})
    ENDIF()
  ENDIF()
  IF(NOT DEFINED ${PROJECT_NAME})
    PROJECT(KokkosKernels CXX)
  ENDIF()
  SET(KokkosKernels_VERSION_MAJOR 3)
  SET(KokkosKernels_VERSION_MINOR 6)
  SET(KokkosKernels_VERSION_PATCH 01)
  SET(KokkosKernels_VERSION "${KokkosKernels_VERSION_MAJOR}.${KokkosKernels_VERSION_MINOR}.${KokkosKernels_VERSION_PATCH}")
  MATH(EXPR KOKKOSKERNELS_VERSION "${KokkosKernels_VERSION_MAJOR} * 10000 + ${KokkosKernels_VERSION_MINOR} * 100 + ${KokkosKernels_VERSION_PATCH}")
ENDIF()

MESSAGE(STATUS "Setting policy CMP0074 to use <Package>_ROOT variables")
CMAKE_POLICY(SET CMP0074 NEW)

INCLUDE(GNUInstallDirs)
IF (KOKKOSKERNELS_HAS_TRILINOS)
 SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
 SET(KOKKOSKERNELS_HEADER_INSTALL_DIR ${TRILINOS_INCDIR})
 SET(KOKKOS_ENABLE_CUDA_UVM ${Kokkos_ENABLE_CUDA_UVM})
ELSEIF(KOKKOSKERNELS_HAS_PARENT)
 SET(KOKKOSKERNELS_HEADER_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/kokkos-kernels")
 SET(KOKKOS_ENABLE_CUDA_UVM ${Kokkos_ENABLE_CUDA_UVM})
ELSE()
 SET(KOKKOSKERNELS_HEADER_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}")
ENDIF()
INCLUDE(cmake/fake_tribits.cmake)
INCLUDE(cmake/kokkoskernels_tribits.cmake)

OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF)

KOKKOSKERNELS_PACKAGE()

IF (NOT KOKKOSKERNELS_HAS_TRILINOS)
  SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
  KOKKOSKERNELS_ADD_OPTION(
    "ENABLE_EXAMPLES"
    OFF
          BOOL
          "Whether to build examples. Default: OFF"
  )
  KOKKOSKERNELS_ADD_OPTION(
          "ENABLE_TESTS"
          OFF
          BOOL
          "Whether to build tests. Default: OFF"
  )
  KOKKOSKERNELS_ADD_OPTION(
          "ENABLE_TESTS_AND_PERFSUITE"
          OFF
          BOOL
          "Whether to build tests including Perfsuite. Default: OFF"
  )
  IF(KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
    set(BLT_CODE_CHECK_TARGET_NAME "fix-for-blt" CACHE STRING "Docstring")
    set(INFRASTRUCTURE_ONLY ON CACHE BOOL "Only build the RAJAPerf infrastructure, no builtin kernels")
    add_definitions("-DRAJAPERF_INFRASTRUCTURE_ONLY")
    add_subdirectory(tpls/rajaperf)
    include_directories(tpls/rajaperf/src)
  ENDIF()
ENDIF ()

KOKKOSKERNELS_ADD_OPTION(
        "ENABLE_DOCS"
        OFF
        BOOL
        "Whether to build docs. Default: OFF"
)

SET(KokkosKernels_INSTALL_TESTING OFF CACHE INTERNAL
        "Whether to build tests and examples against installation")
IF (KokkosKernels_INSTALL_TESTING)
  # Force testing on if we are doing intall testing
  SET(KOKKOSKERNELS_ENABLE_TESTS ON)
  SET(KOKKOSKERNELS_ENABLE_EXAMPLES ON)
  # Don't build, load installed kernels
  FIND_PACKAGE(KokkosKernels REQUIRED)
  # Still need to figure out which backends
  INCLUDE(cmake/kokkos_backends.cmake)
  # Only build the tests
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(perf_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(unit_test)
  KOKKOSKERNELS_ADD_EXAMPLE_DIRECTORIES(example)
ELSE()
  # Regular build, not install testing
  # Do all the regular option processing
  IF (NOT KOKKOSKERNELS_HAS_TRILINOS AND NOT KOKKOSKERNELS_HAS_PARENT)
   # This is a standalone build
   FIND_PACKAGE(Kokkos REQUIRED)
   MESSAGE(STATUS "Found Kokkos at ${Kokkos_DIR}")
   KOKKOS_CHECK(OPTIONS CUDA_UVM RETURN_VALUE KOKKOS_ENABLE_CUDA_UVM)
  ENDIF()

  INCLUDE(cmake/kokkos_backends.cmake)

  #This will collect imported TPLs we need to export
  #in our config file
  GLOBAL_SET(KOKKOSKERNELS_TPL_EXPORTS)

  # If building in debug mode, define the HAVE_KOKKOSKERNELS_DEBUG macro.
  KOKKOSKERNELS_ADD_DEBUG_OPTION()

  #
  # "Optimization level" for KokkosKernels computational kernels.  The
  # higher the level, the more code variants get generated, and thus the
  # longer the compile times.  However, more code variants mean both
  # better performance overall, and more uniform performance for corner
  # cases.  Values of current interest (24 Apr 2014) are 0, 1, and 2.
  #
  KOKKOSKERNELS_ADD_OPTION_AND_DEFINE(
    LINALG_OPT_LEVEL
    KOKKOSLINALG_OPT_LEVEL
    "Optimization level for KokkosKernels computational kernels: a nonnegative integer.  Higher levels result in better performance that is more uniform for corner cases, but increase build time and library size.  The default value is 1, which should give performance within ten percent of optimal on most platforms, for most problems. Default: 1"
    "1")

  # Enable experimental features of KokkosKernels if set at configure
  # time. Default is no.
  KOKKOSKERNELS_ADD_OPTION_AND_DEFINE(
    ENABLE_EXPERIMENTAL
    HAVE_KOKKOSKERNELS_EXPERIMENTAL
    "Enable building and installation of experimental KokkosKernels features. Default: OFF"
    OFF)

  KOKKOSKERNELS_ADD_OPTION(
    ADD_DEFAULT_ETI
    ON
    BOOL
    "Whether to include a set of default ETI instantiations (otherwise only those explicitly requested will be included"
  )

  KOKKOSKERNELS_ADD_OPTION(
    ETI_ONLY
    OFF
    BOOL
    "Whether to restrict availability of kernels to ETI types only. Turning this on guarantees that kernels are never built inside of object files which simply call KokkosKernels functions. Default: OFF"
    )

  KOKKOSKERNELS_ADD_OPTION(
    TEST_ETI_ONLY
    ON
    BOOL
    "Whether to restrict testing to ETI types. Default: ON"
    )

  KOKKOSKERNELS_ADD_OPTION(
    ENABLED_COMPONENTS
    "ALL"
    STRING
    "A list of components to enable in testing and building"
    VALID_ENTRIES BATCHED BLAS GRAPH SPARSE ALL
  )

  # ==================================================================
  # Enable Device Types for ETI (exec- + mem-space)
  # ==================================================================
  INCLUDE(cmake/kokkoskernels_eti_devices.cmake)
  # ==================================================================
  # Enable Scalar Types for ETI
  # ==================================================================
  INCLUDE(cmake/kokkoskernels_eti_floats.cmake)
  # ==================================================================
  # Enable Ordinal Types for ETI
  # ==================================================================
  INCLUDE(cmake/kokkoskernels_eti_ordinals.cmake)
  # ==================================================================
  # Enable Offset Types for ETI
  # ==================================================================
  INCLUDE(cmake/kokkoskernels_eti_offsets.cmake)
  # ==================================================================
  # Enable Layout Types for ETI
  # ==================================================================
  INCLUDE(cmake/kokkoskernels_eti_layouts.cmake)
  # ==================================================================
  # Enable Third Party Libraries
  # ==================================================================
  INCLUDE(cmake/kokkoskernels_tpls.cmake)
  INCLUDE(cmake/kokkoskernels_features.cmake)
  INCLUDE(cmake/kokkos_requirements.cmake)
  # ==================================================================
  # CMake Summary
  # ==================================================================
  MESSAGE("")
  MESSAGE("=======================")
  MESSAGE("KokkosKernels ETI Types")
  MESSAGE("   Devices:  ${DEVICE_LIST}")
  MESSAGE("   Scalars:  ${SCALAR_LIST}")
  MESSAGE("   Ordinals: ${ORDINAL_LIST}")
  MESSAGE("   Offsets:  ${OFFSET_LIST}")
  MESSAGE("   Layouts:  ${LAYOUT_LIST}")
  MESSAGE("")
  MESSAGE("KokkosKernels TPLs")
  FOREACH(TPL ${KOKKOSKERNELS_TPL_LIST})
    PAD_STRING("${TPL}:" TPL_PADDED 12)
    MESSAGE("   ${TPL_PADDED} ${${TPL}_FOUND_INFO}")
  ENDFOREACH()
  MESSAGE("=======================")
  MESSAGE("")
  # Skip building Kokkos Kernels if we are doing an installation test
  ADD_SUBDIRECTORY(src)
  IF (KokkosKernels_ENABLE_INSTALL_TEST)
    ADD_SUBDIRECTORY(install_test)
    MESSAGE("The install test has been enabled, you will need to peform: make install before running the tests otherwise install_test will fail")
  ENDIF ()
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(perf_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(unit_test)
  KOKKOSKERNELS_ADD_EXAMPLE_DIRECTORIES(example)

  KOKKOSKERNELS_PACKAGE_POSTPROCESS()
  IF (KokkosKernels_ENABLE_DOCS)
    ADD_SUBDIRECTORY(docs)
  ENDIF ()
ENDIF()
