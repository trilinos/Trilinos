KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../test_common)

# At some point, we may wish to make this into real "performance
# tests," in the sense that they can pass or fail.  At that point, use
# "CATEGORIES PERFORMANCE" to mark them as such.  For now, we just
# build the executables for manual use, but don't run the tests.  They
# build correctly with or without MPI, but only run them with a single
# MPI process.

SET(GTEST_SOURCE_DIR ${PACKAGE_SOURCE_DIR}/tpls/gtest)

KOKKOSKERNELS_ADD_TEST_LIBRARY(
  kokkoskernelsperf_gtest
  HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
  SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
  )
# Disables pthreads, this is a problem for serial builds in Trilinos & Sierra if it's enabled.
TARGET_COMPILE_DEFINITIONS(kokkoskernelsperf_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0")
TARGET_INCLUDE_DIRECTORIES(kokkoskernelsperf_gtest PUBLIC $<BUILD_INTERFACE:${GTEST_SOURCE_DIR}>)

#Gtest minimally requires C++11
TARGET_COMPILE_FEATURES(kokkoskernelsperf_gtest PUBLIC cxx_std_11)



ADD_COMPONENT_SUBDIRECTORY(batched)
ADD_COMPONENT_SUBDIRECTORY(graph)
ADD_COMPONENT_SUBDIRECTORY(sparse)
ADD_COMPONENT_SUBDIRECTORY(blas)
ADD_SUBDIRECTORY(performance)
#ADD_SUBDIRECTORY(common)

