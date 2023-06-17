# Adding source directory to the build
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched)

LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/src)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/impl)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/unit_test)

LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/src)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/impl)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/unit_test)

IF (NOT KokkosKernels_ENABLE_COMPONENT_BLAS)
  LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas/impl)
  LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas/tpls)
  APPEND_GLOB(HEADERS ${PACKAGE_SOURCE_DIR}/blas/impl/KokkosBlas_util.hpp)
  LIST(APPEND SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/batched/KokkosBatched_Util.cpp)
ENDIF()

IF(KokkosKernels_ENABLE_TESTS OR KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
  # Adding unit-tests
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
    ${CMAKE_CURRENT_SOURCE_DIR}/batched)
ENDIF()
# NOTE: Above, the build directory 'batched' is not created unless unit tests
# are actually enabled (which are actually included from the base-level
# CMakeLists.txt file).  And the KokkosKernelsTargets.cmake file that gets
# generated from this CMake package in the build dir will be broken if these
# are listed in the `INTERFACE_INCLUDE_DIRECTORIES` property when the build
# `batched` is not created (see Trilinos PR #11966).
