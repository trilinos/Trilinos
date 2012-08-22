INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINRARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../TPL)

SET(SOURCES PerfTestMain.cpp PerfTestHost.cpp ../TPL/gtest/gtest-all.cc)
SET(LIBRARIES kokkosarray_devicehost kokkosarray_impl)

IF(KokkosArray_ENABLE_DeviceCuda)
  PACKAGE_ADD_CUDA_LIBRARY(
    perf_cuda
    SOURCES PerfTestCuda.cu
  )
  LIST(APPEND SOURCES
    PerfTestCuda.cpp
  )
  LIST(APPEND LIBRARIES
    kokkosarray_devicecuda perf_cuda
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ArrayPerfTest
  SOURCES ${SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "[  PASSED  ]"
  DEPLIBS ${LIBRARIES}
  )
