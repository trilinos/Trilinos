if (COMMAND TRIBITS_PACKAGE_DECL)
  set(KOKKOSKERNELS_HAS_TRILINOS ON CACHE BOOL "")
else()
  set(KOKKOSKERNELS_HAS_TRILINOS OFF CACHE BOOL "")
endif()

get_directory_property(KOKKOSKERNELS_HAS_PARENT PARENT_DIRECTORY)

set(KOKKOSKERNELS_TOP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(KOKKOSKERNELS_TOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

SET(KokkosKernels_VERSION_MAJOR 5)
SET(KokkosKernels_VERSION_MINOR 0)
SET(KokkosKernels_VERSION_PATCH 2)
SET(KokkosKernels_VERSION "${KokkosKernels_VERSION_MAJOR}.${KokkosKernels_VERSION_MINOR}.${KokkosKernels_VERSION_PATCH}")

#Set variables for config file
math(EXPR KOKKOSKERNELS_VERSION "${KokkosKernels_VERSION_MAJOR} * 10000 + ${KokkosKernels_VERSION_MINOR} * 100 + ${KokkosKernels_VERSION_PATCH}")
math(EXPR KOKKOSKERNELS_VERSION_MAJOR "${KOKKOSKERNELS_VERSION} / 10000")
math(EXPR KOKKOSKERNELS_VERSION_MINOR "${KOKKOSKERNELS_VERSION} / 100 % 100")
math(EXPR KOKKOSKERNELS_VERSION_PATCH "${KOKKOSKERNELS_VERSION} % 100")

if (NOT KOKKOSKERNELS_HAS_TRILINOS)
  cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
  if (Spack_WORKAROUND)
    #if we are explicitly using Spack for development,
    #nuke the Spack compiler
    set(SPACK_CXX $ENV{SPACK_CXX})
    if (SPACK_CXX)
      set(CMAKE_CXX_COMPILER ${SPACK_CXX} CACHE STRING "the C++ compiler" FORCE)
      set(ENV{CXX} ${SPACK_CXX})
    endif()
  endif()
  if (NOT DEFINED ${PROJECT_NAME})
    project(KokkosKernels VERSION ${KokkosKernels_VERSION} LANGUAGES CXX)
  endif()
endif()

include(GNUInstallDirs)
if (KOKKOSKERNELS_HAS_TRILINOS)
  set(TRILINOS_INCDIR ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
  set(KOKKOSKERNELS_HEADER_INSTALL_DIR ${TRILINOS_INCDIR})
elseif (KOKKOSKERNELS_HAS_PARENT)
  set(KOKKOSKERNELS_HEADER_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/kokkos-kernels")
else()
  set(KOKKOSKERNELS_HEADER_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}")
endif()
include(cmake/fake_tribits.cmake)
include(cmake/kokkoskernels_tribits.cmake)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

kokkoskernels_package()

if (NOT KOKKOSKERNELS_HAS_TRILINOS)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
  kokkoskernels_add_option(
    "ENABLE_EXAMPLES"
    OFF
    BOOL
    "Whether to build examples. Default: OFF"
  )
  kokkoskernels_add_option(
    "ENABLE_TESTS"
    OFF
    BOOL
    "Whether to build tests. Default: OFF"
  )
  kokkoskernels_add_option(
    "ENABLE_PERFTESTS"
    OFF
    BOOL
    "Whether to build performance tests. Default: OFF"
  )
  kokkoskernels_add_option(
    "ENABLE_TESTS_AND_PERFSUITE"
    OFF
    BOOL
    "Whether to build performance tests and suite. Default: OFF"
  )
  if (KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
    set(BLT_CODE_CHECK_TARGET_NAME "fix-for-blt" CACHE STRING "Docstring")
    set(INFRASTRUCTURE_ONLY ON CACHE BOOL "Only build the RAJAPerf infrastructure, no builtin kernels")
    add_definitions("-DRAJAPERF_INFRASTRUCTURE_ONLY")
    add_subdirectory(tpls/rajaperf)
    include_directories(tpls/rajaperf/src)
    set(KokkosKernels_ENABLE_PERFTESTS ON CACHE BOOL "Whether to build tests including Perfsuite. Default: OFF" FORCE)
  endif()
  kokkoskernels_add_option(
    "ENABLE_BENCHMARKS"
    OFF
    BOOL
    "Whether to build benchmarks. Default: OFF"
  )
  if (KokkosKernels_ENABLE_BENCHMARKS)
    include(cmake/kokkoskernels_benchmarks.cmake)
  endif()
  kokkoskernels_add_option(
    "RUN_BENCHMARKS"
    OFF
    BOOL
    "Whether to run benchmarks. Default: OFF"
  )
endif()

kokkoskernels_add_option(
  "ENABLE_DOCS"
  OFF
  BOOL
  "Whether to build docs. Default: OFF"
)

kokkoskernels_add_option(
  "ENABLE_COMPILER_WARNINGS"
  OFF
  BOOL
  "Whether to turn on more warnings. Default: OFF"
)

kokkoskernels_add_option(
  "ENABLE_WARNINGS_AS_ERRORS"
  OFF
  BOOL
  "Whether warnings should be treated as errors. Default: OFF"
)

include(cmake/kokkoskernels_warnings.cmake)

set(KokkosKernels_INSTALL_TESTING OFF CACHE INTERNAL
  "Whether to build tests and examples against installation")
if (KokkosKernels_INSTALL_TESTING)
  # Force testing on if we are doing intall testing
  set(KOKKOSKERNELS_ENABLE_TESTS ON)
  set(KOKKOSKERNELS_ENABLE_EXAMPLES ON)
  # Don't build, load installed kernels
  find_package(KokkosKernels REQUIRED)
  # Still need to figure out which backends
  include(cmake/kokkos_backends.cmake)
  # Only build the tests
  # Assume that all components have been enabled
  kokkoskernels_add_test_directories(common/unit_test)
  kokkoskernels_add_test_directories(batched/dense/unit_test)
  kokkoskernels_add_test_directories(batched/sparse/unit_test)
  kokkoskernels_add_test_directories(blas/unit_test)
  kokkoskernels_add_test_directories(lapack/unit_test)
  kokkoskernels_add_test_directories(graph/unit_test)
  kokkoskernels_add_test_directories(sparse/unit_test)
  kokkoskernels_add_test_directories(ode/unit_test)
  kokkoskernels_add_test_directories(perf_test)
  kokkoskernels_add_example_directories(example)
else()
  # Regular build, not install testing
  # Do all the regular option processing
  if (NOT KOKKOSKERNELS_HAS_TRILINOS AND NOT KOKKOSKERNELS_HAS_PARENT)
    # This is a standalone build
    find_package(Kokkos REQUIRED)
    if (Kokkos_COMPILE_LANGUAGE)
      enable_language(${Kokkos_COMPILE_LANGUAGE})
      if (Kokkos_COMPILE_LANGUAGE STREQUAL CUDA)
        set(CMAKE_CUDA_ARCHITECTURES ${Kokkos_CUDA_ARCHITECTURES})
      endif()
      if (Kokkos_COMPILE_LANGUAGE STREQUAL HIP)
        set(CMAKE_HIP_ARCHITECTURES ${Kokkos_HIP_ARCHITECTURES})
      endif()
    endif()
    if(${Kokkos_VERSION} VERSION_GREATER_EQUAL "4.7.02")
      message(STATUS "Found Kokkos version ${Kokkos_VERSION} at ${Kokkos_DIR}")
      if ((${Kokkos_VERSION} VERSION_GREATER "5.0.99"))
        message(WARNING "Configuring with Kokkos ${Kokkos_VERSION} which is newer than the expected develop branch - version check may need update")
      endif()
    else()
      message(FATAL_ERROR "Kokkos Kernels ${KokkosKernels_VERSION} requires Kokkos 4.7.02 or greater (found ${Kokkos_VERSION})")
    endif()

    if (KokkosKernels_ENABLE_COMPILER_WARNINGS)
      kokkoskernels_enable_warnings()
    endif()
  endif()

  include(cmake/kokkos_backends.cmake)

  #This will collect imported TPLs we need to export
  #in our config file
  global_set(KOKKOSKERNELS_TPL_EXPORTS)

  # If building in debug mode, define the HAVE_KOKKOSKERNELS_DEBUG macro.
  kokkoskernels_add_debug_option()

  #
  # "Optimization level" for KokkosKernels computational kernels.  The
  # higher the level, the more code variants get generated, and thus the
  # longer the compile times.  However, more code variants mean both
  # better performance overall, and more uniform performance for corner
  # cases.  Values of current interest (24 Apr 2014) are 0, 1, and 2.
  #
  kokkoskernels_add_option_and_define(
    LINALG_OPT_LEVEL
    KOKKOSLINALG_OPT_LEVEL
    "DEPRECATED. Optimization level for KokkosKernels computational kernels: a nonnegative integer.  Higher levels result in better performance that is more uniform for corner cases, but increase build time and library size.  The default value is 1, which should give performance within ten percent of optimal on most platforms, for most problems. Default: 1"
    "1")

  if (KokkosKernels_LINALG_OPT_LEVEL AND NOT KokkosKernels_LINALG_OPT_LEVEL STREQUAL "1")
    message(WARNING "KokkosKernels_LINALG_OPT_LEVEL is deprecated!")
  endif()
  if (KokkosKernels_KOKKOSLINALG_OPT_LEVEL AND NOT KokkosKernels_KOKKOSLINALG_OPT_LEVEL STREQUAL "1")
    message(WARNING "KokkosKernels_KOKKOSLINALG_OPT_LEVEL is deprecated!")
  endif()

  # Enable experimental features of KokkosKernels if set at configure
  # time. Default is no.
  kokkoskernels_add_option_and_define(
    ENABLE_EXPERIMENTAL
    HAVE_KOKKOSKERNELS_EXPERIMENTAL
    "Enable building and installation of experimental KokkosKernels features. Default: OFF"
    OFF)

  kokkoskernels_add_option(
    ADD_DEFAULT_ETI
    ON
    BOOL
    "Whether to include a set of default ETI instantiations (otherwise only those explicitly requested will be included"
  )

  kokkoskernels_add_option(
    ETI_ONLY
    OFF
    BOOL
    "Whether to restrict availability of kernels to ETI types only. Turning this on guarantees that kernels are never built inside of object files which simply call KokkosKernels functions. Default: OFF"
  )

  kokkoskernels_add_option(
    TEST_ETI_ONLY
    ON
    BOOL
    "Whether to restrict testing to ETI types. Default: ON"
  )

  kokkoskernels_add_option(
    ENABLED_COMPONENTS
    "ALL"
    STRING
    "A list of components to enable in testing and building"
    VALID_ENTRIES BATCHED BLAS LAPACK GRAPH SPARSE ALL
  )

  # ==================================================================
  # Enable Device Types for ETI (exec- + mem-space)
  # ==================================================================
  include(cmake/kokkoskernels_eti_devices.cmake)
  # ==================================================================
  # Enable Scalar Types for ETI
  # ==================================================================
  include(cmake/kokkoskernels_eti_floats.cmake)
  # ==================================================================
  # Enable Ordinal Types for ETI
  # ==================================================================
  include(cmake/kokkoskernels_eti_ordinals.cmake)
  # ==================================================================
  # Enable Offset Types for ETI
  # ==================================================================
  include(cmake/kokkoskernels_eti_offsets.cmake)
  # ==================================================================
  # Enable Layout Types for ETI
  # ==================================================================
  include(cmake/kokkoskernels_eti_layouts.cmake)
  # ==================================================================
  # Enable Third Party Libraries
  # ==================================================================
  include(cmake/kokkoskernels_tpls.cmake)
  include(cmake/kokkoskernels_features.cmake)
  include(cmake/kokkos_requirements.cmake)
  # ==================================================================
  # Enable Library Components
  # ==================================================================
  include(cmake/kokkoskernels_components.cmake)
  # ==================================================================
  # CMake Summary
  # ==================================================================
  message("")
  message("================================")
  message("Kokkos Kernels version: ${KOKKOSKERNELS_VERSION_MAJOR}.${KOKKOSKERNELS_VERSION_MINOR}.${KOKKOSKERNELS_VERSION_PATCH}")
  message("================================")
  message("Kokkos Kernels ETI Types")
  message("   Devices:  ${DEVICE_LIST}")
  message("   Scalars:  ${SCALAR_LIST}")
  message("   Ordinals: ${ORDINAL_LIST}")
  message("   Offsets:  ${OFFSET_LIST}")
  message("   Layouts:  ${LAYOUT_LIST}")
  message("")
  message("Kokkos Kernels components")
  message("   COMMON:    ON")
  message("   BATCHED:   ${KokkosKernels_ENABLE_COMPONENT_BATCHED}")
  message("   BLAS:      ${KokkosKernels_ENABLE_COMPONENT_BLAS}")
  message("   LAPACK:    ${KokkosKernels_ENABLE_COMPONENT_LAPACK}")
  message("   GRAPH:     ${KokkosKernels_ENABLE_COMPONENT_GRAPH}")
  message("   SPARSE:    ${KokkosKernels_ENABLE_COMPONENT_SPARSE}")
  message("   ODE:       ${KokkosKernels_ENABLE_COMPONENT_ODE}")
  message("")
  message("Kokkos Kernels TPLs")
  if (KOKKOSKERNELS_TPL_LIST)
    foreach(TPL ${KOKKOSKERNELS_TPL_LIST})
      pad_string("${TPL}:" TPL_PADDED 12)
      message("   ${TPL_PADDED} ${${TPL}_FOUND_INFO}")
    endforeach()
  else()
    message("   (None)")
  endif()
  message("================================")
  message("")

  ## configure the library
  kokkoskernels_configure_file(KokkosKernels_config.h)

  set(KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR})

  #I don't want the complexity of sometimes interface/sometimes not
  #Just minimally force an empty dummy file
  set(HEADERS)
  set(SOURCES common/src/dummy.cpp)

  # ETI infrastructure to be used by components
  include(cmake/kokkoskernels_eti.cmake)
  set(ETI_HEADERS)

  list(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)

  # Only enable requested components
  # if nothing requested, build everything
  # to recover original behavior.
  include(common/CMakeLists.txt)
  if (KokkosKernels_ENABLE_COMPONENT_BATCHED)
    include(batched/CMakeLists.txt)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_BLAS)
    include(blas/CMakeLists.txt)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_LAPACK)
    include(lapack/CMakeLists.txt)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_GRAPH)
    include(graph/CMakeLists.txt)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_SPARSE)
    include(sparse/CMakeLists.txt)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_ODE)
    include(ode/CMakeLists.txt)
  endif()

  foreach(DIR ${KK_INCLUDE_DIRS})
    kokkoskernels_include_directories(${DIR})
    append_glob(HEADERS ${DIR}/*.hpp)
  endforeach()
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
  install(FILES ${ETI_HEADERS} DESTINATION ${KOKKOSKERNELS_HEADER_INSTALL_DIR}/generated_specializations_hpp/)

  if (KOKKOSKERNELS_HAS_TRILINOS)
    #no linking commands required - tribits does this
  else()
    # Address kokkos/kokkos-kernels#1749
    add_library(KokkosKernels::kokkoskernels ALIAS kokkoskernels)
    # all_libs target is required for TriBITS-compliance
    add_library(KokkosKernels::all_libs ALIAS kokkoskernels)
    target_link_libraries(kokkoskernels PUBLIC Kokkos::kokkos)
    foreach(DIR ${KK_INCLUDE_DIRS})
      target_include_directories(kokkoskernels PUBLIC $<BUILD_INTERFACE:${DIR}>)
    endforeach()
    target_include_directories(kokkoskernels PUBLIC
      $<INSTALL_INTERFACE:${KOKKOSKERNELS_HEADER_INSTALL_DIR}>)

    if (KokkosKernels_ENABLE_BENCHMARKS)
      include(cmake/kokkoskernels_version_info.cmake)
      check_version_info()
    endif()
  endif()

  # FIXME_SYCL waiting for compiler support
  if (KOKKOS_ENABLE_SYCL)
    set(KOKKOSKERNELS_INTEL_ARCHS ${Kokkos_ARCH})
    list(FILTER KOKKOSKERNELS_INTEL_ARCHS INCLUDE REGEX ".*INTEL.*")
    list(LENGTH KOKKOSKERNELS_INTEL_ARCHS KOKKOSKERNELS_INTEL_ARCHS_NUM)
    # Trilinos doesn't see Kokkos_ARCH but sees Kokkos CMake configuration.
    if (KOKKOSKERNELS_INTEL_ARCHS_NUM GREATER_EQUAL 1 OR Kokkos_ARCH_INTEL_XEHP OR Kokkos_ARCH_INTEL_PVC)
      if (NOT BUILD_SHARED_LIBS)
        message(SEND_ERROR
          "At the moment, we require KokkosKernels (and Kokkos) to be built as "
          "shared libraries to allow querying free and total device memory!"
          )
      endif()
      target_link_libraries(kokkoskernels PUBLIC ze_loader)
    endif()
  endif()

  kokkoskernels_link_tpl(kokkoskernels PUBLIC LAPACK)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC BLAS)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC CBLAS)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC LAPACKE)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC SUPERLU)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC CHOLMOD)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC MKL)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC CUBLAS)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC CUSPARSE)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC CUSOLVER)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC ROCBLAS)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC ROCSPARSE)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC ROCSOLVER)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC METIS)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC ARMPL)
  kokkoskernels_link_tpl(kokkoskernels PUBLIC MAGMA)

  if (KokkosKernels_ENABLE_INSTALL_TEST)
    add_subdirectory(install_test)
    message("The install test has been enabled, you will need to peform: make install before running the tests otherwise install_test will fail")
  endif ()
  if (KokkosKernels_ENABLE_TESTS OR KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
    set(GTEST_SOURCE_DIR ${KOKKOSKERNELS_TOP_SOURCE_DIR}/tpls/gtest)

    kokkoskernels_add_test_library(
      kokkoskernels_gtest
      HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
      SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
    )
    # Disables pthreads, this is a problem for serial builds in Trilinos & Sierra if it's enabled.
    target_compile_definitions(kokkoskernels_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0")
    target_include_directories(kokkoskernels_gtest PUBLIC $<BUILD_INTERFACE:${GTEST_SOURCE_DIR}>)

    #Gtest minimally requires C++11
    target_compile_features(kokkoskernels_gtest PUBLIC cxx_std_11)
  endif()
  kokkoskernels_add_test_directories(common/unit_test)
  if (KokkosKernels_ENABLE_COMPONENT_BATCHED)
    kokkoskernels_add_test_directories(batched/dense/unit_test)
    kokkoskernels_add_test_directories(batched/sparse/unit_test)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_BLAS)
    kokkoskernels_add_test_directories(blas/unit_test)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_LAPACK)
    kokkoskernels_add_test_directories(lapack/unit_test)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_GRAPH)
    kokkoskernels_add_test_directories(graph/unit_test)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_SPARSE)
    kokkoskernels_add_test_directories(sparse/unit_test)
  endif()
  if (KokkosKernels_ENABLE_COMPONENT_ODE)
    kokkoskernels_add_test_directories(ode/unit_test)
  endif()
  if (KOKKOSKERNELS_ALL_COMPONENTS_ENABLED)
    if (KokkosKernels_ENABLE_PERFTESTS)
      message(STATUS "Enabling perf tests.")
      add_subdirectory(perf_test) # doesn't require KokkosKernels_ENABLE_TESTS=ON
    endif()
    if (KokkosKernels_ENABLE_EXAMPLES)
      message(STATUS "Enabling examples.")
      kokkoskernels_add_example_directories(example)
    endif()
    if (KokkosKernels_ENABLE_BENCHMARKS)
      message(STATUS "Enabling benchmarks.")
      add_subdirectory(benchmarks)
    endif()
  else()
    # all components were not enabled, so perftests and examples can't be enabled.
    # Warn if they were requested.
    if (KokkosKernels_ENABLE_PERFTESTS)
      message(WARNING "Could not enable perf tests because not all components were enabled")
    endif()
    if (KokkosKernels_ENABLE_EXAMPLES)
      message(WARNING "Could not enable examples because not all components were enabled")
    endif()
  endif()

  kokkoskernels_package_postprocess()
  if (KokkosKernels_ENABLE_DOCS)
    add_subdirectory(docs)
  endif()
endif()
