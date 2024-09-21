IF(COMMAND TRIBITS_PACKAGE_DECL)
  SET(KOKKOSKERNELS_HAS_TRILINOS ON CACHE BOOL "")
ELSE()
  SET(KOKKOSKERNELS_HAS_TRILINOS OFF CACHE BOOL "")
ENDIF()

GET_DIRECTORY_PROPERTY(KOKKOSKERNELS_HAS_PARENT PARENT_DIRECTORY)

SET(KOKKOSKERNELS_TOP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(KOKKOSKERNELS_TOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

SET(KokkosKernels_VERSION_MAJOR 4)
SET(KokkosKernels_VERSION_MINOR 4)
SET(KokkosKernels_VERSION_PATCH 1)
SET(KokkosKernels_VERSION "${KokkosKernels_VERSION_MAJOR}.${KokkosKernels_VERSION_MINOR}.${KokkosKernels_VERSION_PATCH}")

#Set variables for config file
MATH(EXPR KOKKOSKERNELS_VERSION "${KokkosKernels_VERSION_MAJOR} * 10000 + ${KokkosKernels_VERSION_MINOR} * 100 + ${KokkosKernels_VERSION_PATCH}")
MATH(EXPR KOKKOSKERNELS_VERSION_MAJOR "${KOKKOSKERNELS_VERSION} / 10000")
MATH(EXPR KOKKOSKERNELS_VERSION_MINOR "${KOKKOSKERNELS_VERSION} / 100 % 100")
MATH(EXPR KOKKOSKERNELS_VERSION_PATCH "${KOKKOSKERNELS_VERSION} % 100")

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
ENDIF()

INCLUDE(GNUInstallDirs)
IF (KOKKOSKERNELS_HAS_TRILINOS)
 SET(TRILINOS_INCDIR ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
 SET(KOKKOSKERNELS_HEADER_INSTALL_DIR ${TRILINOS_INCDIR})
ELSEIF(KOKKOSKERNELS_HAS_PARENT)
 SET(KOKKOSKERNELS_HEADER_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/kokkos-kernels")
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
          "ENABLE_PERFTESTS"
          OFF
          BOOL
          "Whether to build performance tests. Default: OFF"
  )
  KOKKOSKERNELS_ADD_OPTION(
          "ENABLE_TESTS_AND_PERFSUITE"
          OFF
          BOOL
          "Whether to build performance tests and suite. Default: OFF"
  )
  IF(KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
    set(BLT_CODE_CHECK_TARGET_NAME "fix-for-blt" CACHE STRING "Docstring")
    set(INFRASTRUCTURE_ONLY ON CACHE BOOL "Only build the RAJAPerf infrastructure, no builtin kernels")
    add_definitions("-DRAJAPERF_INFRASTRUCTURE_ONLY")
    add_subdirectory(tpls/rajaperf)
    include_directories(tpls/rajaperf/src)
    set(KokkosKernels_ENABLE_PERFTESTS ON CACHE BOOL "Whether to build tests including Perfsuite. Default: OFF" FORCE)
  ENDIF()
  IF(KokkosKernels_ENABLE_BENCHMARK)
    SET(KOKKOSKERNELS_ENABLE_BENCHMARK ON CACHE BOOL "Benchmark enabled")
    INCLUDE(cmake/kokkoskernels_benchmarks.cmake)
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
  # Assume that all components have been enabled
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(common/unit_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(batched/dense/unit_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(batched/sparse/unit_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(blas/unit_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(lapack/unit_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(graph/unit_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(sparse/unit_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(ode/unit_test)
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(perf_test)
  KOKKOSKERNELS_ADD_EXAMPLE_DIRECTORIES(example)
ELSE()
  # Regular build, not install testing
  # Do all the regular option processing
  IF (NOT KOKKOSKERNELS_HAS_TRILINOS AND NOT KOKKOSKERNELS_HAS_PARENT)
    # This is a standalone build
    FIND_PACKAGE(Kokkos REQUIRED)
    IF(Kokkos_COMPILE_LANGUAGE)
      ENABLE_LANGUAGE(${Kokkos_COMPILE_LANGUAGE})
      IF(Kokkos_COMPILE_LANGUAGE STREQUAL CUDA)
        SET(CMAKE_CUDA_ARCHITECTURES ${Kokkos_CUDA_ARCHITECTURES})
      ENDIF()
      IF(Kokkos_COMPILE_LANGUAGE STREQUAL HIP)
        SET(CMAKE_HIP_ARCHITECTURES ${Kokkos_HIP_ARCHITECTURES})
      ENDIF()
    ENDIF()
    IF(${Kokkos_VERSION} VERSION_GREATER_EQUAL "4.3.01")
      MESSAGE(STATUS "Found Kokkos version ${Kokkos_VERSION} at ${Kokkos_DIR}")
      IF((${Kokkos_VERSION} VERSION_GREATER "4.4.99"))
        MESSAGE(WARNING "Configuring with Kokkos ${Kokkos_VERSION} which is newer than the expected develop branch - version check may need update")
      ENDIF()
    ELSE()
      MESSAGE(FATAL_ERROR "Kokkos Kernels ${KokkosKernels_VERSION} requires Kokkos 4.3.01 or greater (found ${Kokkos_VERSION})")
    ENDIF()
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
    "DEPRECATED. Optimization level for KokkosKernels computational kernels: a nonnegative integer.  Higher levels result in better performance that is more uniform for corner cases, but increase build time and library size.  The default value is 1, which should give performance within ten percent of optimal on most platforms, for most problems. Default: 1"
    "1")

  if (KokkosKernels_LINALG_OPT_LEVEL AND NOT KokkosKernels_LINALG_OPT_LEVEL STREQUAL "1")
   message(WARNING "KokkosKernels_LINALG_OPT_LEVEL is deprecated!")
 endif()
  if(KokkosKernels_KOKKOSLINALG_OPT_LEVEL AND NOT KokkosKernels_KOKKOSLINALG_OPT_LEVEL STREQUAL "1")
    message(WARNING "KokkosKernels_KOKKOSLINALG_OPT_LEVEL is deprecated!")
  endif()

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
    VALID_ENTRIES BATCHED BLAS LAPACK GRAPH SPARSE ALL
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
  # Enable Library Components
  # ==================================================================
  INCLUDE(cmake/kokkoskernels_components.cmake)
  # ==================================================================
  # CMake Summary
  # ==================================================================
  MESSAGE("")
  MESSAGE("================================")
  MESSAGE("Kokkos Kernels version: ${KOKKOSKERNELS_VERSION_MAJOR}.${KOKKOSKERNELS_VERSION_MINOR}.${KOKKOSKERNELS_VERSION_PATCH}")
  MESSAGE("================================")
  MESSAGE("Kokkos Kernels ETI Types")
  MESSAGE("   Devices:  ${DEVICE_LIST}")
  MESSAGE("   Scalars:  ${SCALAR_LIST}")
  MESSAGE("   Ordinals: ${ORDINAL_LIST}")
  MESSAGE("   Offsets:  ${OFFSET_LIST}")
  MESSAGE("   Layouts:  ${LAYOUT_LIST}")
  MESSAGE("")
  MESSAGE("Kokkos Kernels components")
  MESSAGE("   COMMON:    ON")
  MESSAGE("   BATCHED:   ${KokkosKernels_ENABLE_COMPONENT_BATCHED}")
  MESSAGE("   BLAS:      ${KokkosKernels_ENABLE_COMPONENT_BLAS}")
  MESSAGE("   LAPACK:    ${KokkosKernels_ENABLE_COMPONENT_LAPACK}")
  MESSAGE("   GRAPH:     ${KokkosKernels_ENABLE_COMPONENT_GRAPH}")
  MESSAGE("   SPARSE:    ${KokkosKernels_ENABLE_COMPONENT_SPARSE}")
  MESSAGE("   ODE:       ${KokkosKernels_ENABLE_COMPONENT_ODE}")
  MESSAGE("")
  MESSAGE("Kokkos Kernels TPLs")
  IF(KOKKOSKERNELS_TPL_LIST)
    FOREACH(TPL ${KOKKOSKERNELS_TPL_LIST})
      PAD_STRING("${TPL}:" TPL_PADDED 12)
      MESSAGE("   ${TPL_PADDED} ${${TPL}_FOUND_INFO}")
    ENDFOREACH()
  ELSE()
     MESSAGE("   (None)")
  ENDIF()
  MESSAGE("================================")
  MESSAGE("")



  ## configure the library
  KOKKOSKERNELS_CONFIGURE_FILE(KokkosKernels_config.h)

  SET(KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR})

  #I don't want the complexity of sometimes interface/sometimes not
  #Just minimally force an empty dummy file
  SET(HEADERS)
  SET(SOURCES common/src/dummy.cpp)

  # ETI infrastructure to be used by components
  include(cmake/kokkoskernels_eti.cmake)
  SET(ETI_HEADERS)

  LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)

  # Only enable requested components
  # if nothing requested, build everything
  # to recover original behavior.
  INCLUDE(common/CMakeLists.txt)
  IF (KokkosKernels_ENABLE_COMPONENT_BATCHED)
    INCLUDE(batched/CMakeLists.txt)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_BLAS)
    INCLUDE(blas/CMakeLists.txt)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_LAPACK)
    INCLUDE(lapack/CMakeLists.txt)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_GRAPH)
    INCLUDE(graph/CMakeLists.txt)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_SPARSE)
    INCLUDE(sparse/CMakeLists.txt)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_ODE)
    INCLUDE(ode/CMakeLists.txt)
  ENDIF()

  FOREACH(DIR ${KK_INCLUDE_DIRS})
    KOKKOSKERNELS_INCLUDE_DIRECTORIES(${DIR})
    APPEND_GLOB(HEADERS ${DIR}/*.hpp)
  ENDFOREACH()
  # MESSAGE("HEADERS: ${HEADERS}")


  # For now, don't add the ETI headers to complete list of headers
  # This will break some IDEs and won't properly display headers in the file list
  # However, because of the behavior of TRIBITS_ADD_LIBRARY
  # I can't add the ETI headers with them getting double-installed
  # with a bunch of headers in the wrong location
  # This doesn't change pre-existing behavior before the ETI changes
  #LIST(APPEND HEADERS ${ETI_HEADERS})
  #-----------------------------------------------------------------------------
  KOKKOSKERNELS_ADD_LIBRARY(
    kokkoskernels
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
  )

  # It does not seem to be possible to automatically have headers installed
  # in the correct subdirectory structure using HEADERS and ADD_LIBRARY with tribits
  # Force the specializations to be installed in a subdirectory
  INSTALL(FILES ${ETI_HEADERS} DESTINATION ${KOKKOSKERNELS_HEADER_INSTALL_DIR}/generated_specializations_hpp/)

  IF (KOKKOSKERNELS_HAS_TRILINOS)
    #no linking commands required - tribits does this
  ELSE()
    ADD_LIBRARY(Kokkos::kokkoskernels ALIAS kokkoskernels)
    # Address kokkos/kokkos-kernels#1749
    ADD_LIBRARY(KokkosKernels::kokkoskernels ALIAS kokkoskernels)
    # all_libs target is required for TriBITS-compliance
    ADD_LIBRARY(KokkosKernels::all_libs ALIAS kokkoskernels)
    TARGET_LINK_LIBRARIES(kokkoskernels PUBLIC Kokkos::kokkos)
    FOREACH(DIR ${KK_INCLUDE_DIRS})
      TARGET_INCLUDE_DIRECTORIES(kokkoskernels PUBLIC $<BUILD_INTERFACE:${DIR}>)
    ENDFOREACH()
    TARGET_INCLUDE_DIRECTORIES(kokkoskernels PUBLIC
      $<INSTALL_INTERFACE:${KOKKOSKERNELS_HEADER_INSTALL_DIR}>)

    IF(KokkosKernels_ENABLE_BENCHMARK)
      INCLUDE(cmake/kokkoskernels_version_info.cmake)
      check_version_info()
    ENDIF()
  ENDIF()

  # FIXME_SYCL waiting for compiler support
  IF (KOKKOS_ENABLE_SYCL)
    SET(KOKKOSKERNELS_INTEL_ARCHS ${Kokkos_ARCH})
    LIST(FILTER KOKKOSKERNELS_INTEL_ARCHS INCLUDE REGEX ".*INTEL.*")
    LIST(LENGTH KOKKOSKERNELS_INTEL_ARCHS KOKKOSKERNELS_INTEL_ARCHS_NUM)
    # Trilinos doesn't see Kokkos_ARCH but sees Kokkos CMake configuration.
    IF(KOKKOSKERNELS_INTEL_ARCHS_NUM GREATER_EQUAL 1 OR Kokkos_ARCH_INTEL_XEHP OR Kokkos_ARCH_INTEL_PVC)
      IF (NOT BUILD_SHARED_LIBS)
        MESSAGE(SEND_ERROR
          "At the moment, we require KokkosKernels (and Kokkos) to be built as "
          "shared libraries to allow querying free and total device memory!"
        )
      ENDIF()
    TARGET_LINK_LIBRARIES(kokkoskernels PUBLIC ze_loader)
    ENDIF()
  ENDIF()

  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC LAPACK)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC BLAS)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC CBLAS)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC LAPACKE)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC SUPERLU)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC CHOLMOD)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC MKL)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC CUBLAS)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC CUSPARSE)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC CUSOLVER)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC ROCBLAS)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC ROCSPARSE)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC ROCSOLVER)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC METIS)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC ARMPL)
  KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC MAGMA)



  IF (KokkosKernels_ENABLE_INSTALL_TEST)
    ADD_SUBDIRECTORY(install_test)
    MESSAGE("The install test has been enabled, you will need to peform: make install before running the tests otherwise install_test will fail")
  ENDIF ()
  IF(KokkosKernels_ENABLE_TESTS OR KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
    SET(GTEST_SOURCE_DIR ${KOKKOSKERNELS_TOP_SOURCE_DIR}/tpls/gtest)

    KOKKOSKERNELS_ADD_TEST_LIBRARY(
      kokkoskernels_gtest
      HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
      SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
      )
    # Disables pthreads, this is a problem for serial builds in Trilinos & Sierra if it's enabled.
    TARGET_COMPILE_DEFINITIONS(kokkoskernels_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0")
    TARGET_INCLUDE_DIRECTORIES(kokkoskernels_gtest PUBLIC $<BUILD_INTERFACE:${GTEST_SOURCE_DIR}>)

    #Gtest minimally requires C++11
    TARGET_COMPILE_FEATURES(kokkoskernels_gtest PUBLIC cxx_std_11)
  ENDIF()
  KOKKOSKERNELS_ADD_TEST_DIRECTORIES(common/unit_test)
  IF (KokkosKernels_ENABLE_COMPONENT_BATCHED)
    KOKKOSKERNELS_ADD_TEST_DIRECTORIES(batched/dense/unit_test)
    KOKKOSKERNELS_ADD_TEST_DIRECTORIES(batched/sparse/unit_test)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_BLAS)
    KOKKOSKERNELS_ADD_TEST_DIRECTORIES(blas/unit_test)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_LAPACK)
    KOKKOSKERNELS_ADD_TEST_DIRECTORIES(lapack/unit_test)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_GRAPH)
    KOKKOSKERNELS_ADD_TEST_DIRECTORIES(graph/unit_test)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_SPARSE)
    KOKKOSKERNELS_ADD_TEST_DIRECTORIES(sparse/unit_test)
  ENDIF()
  IF (KokkosKernels_ENABLE_COMPONENT_ODE)
    KOKKOSKERNELS_ADD_TEST_DIRECTORIES(ode/unit_test)
  ENDIF()
  IF (KOKKOSKERNELS_ALL_COMPONENTS_ENABLED)
    IF (KokkosKernels_ENABLE_PERFTESTS)
      MESSAGE(STATUS "Enabling perf tests.")
      add_subdirectory(perf_test) # doesn't require KokkosKernels_ENABLE_TESTS=ON
    ENDIF ()
    IF (KokkosKernels_ENABLE_EXAMPLES)
      MESSAGE(STATUS "Enabling examples.")
      KOKKOSKERNELS_ADD_EXAMPLE_DIRECTORIES(example)
    ENDIF ()
  ELSE ()
    # all components were not enabled, so perftests and examples can't be enabled.
    # Warn if they were requested.
    IF (KokkosKernels_ENABLE_PERFTESTS)
      MESSAGE(WARNING "Could not enable perf tests because not all components were enabled")
    ENDIF ()
    IF (KokkosKernels_ENABLE_EXAMPLES)
      MESSAGE(WARNING "Could not enable examples because not all components were enabled")
    ENDIF ()
  ENDIF ()

  KOKKOSKERNELS_PACKAGE_POSTPROCESS()
  IF (KokkosKernels_ENABLE_DOCS)
    ADD_SUBDIRECTORY(docs)
  ENDIF ()
ENDIF()
