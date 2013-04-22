INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../TPL)

SET(SOURCES
  UnitTestMain.cpp 
  TestHost.cpp 
  TestTileHost.cpp 
  ../TPL/gtest/gtest-all.cc
  )

SET(LIBRARIES kokkosarray_devicehost kokkosarray_impl)

IF(KokkosArray_ENABLE_DeviceOpenMP)
  LIST( APPEND SOURCES
    TestOpenMP.cpp
  )
ENDIF()

IF(KokkosArray_ENABLE_DeviceCuda)

  TRIBITS_ADD_LIBRARY(
    unit_cuda
    SOURCES TestCudaFunctions.cu TestTileCuda.cu
    DEPLIBS
    CUDALIBRARY
  )
  LIST( APPEND SOURCES
    TestCuda.cpp
  )
  LIST( APPEND LIBRARIES
    kokkosarray_devicecuda unit_cuda
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ArrayUnitTest
  SOURCES ${SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "[  PASSED  ]"
  DEPLIBS ${LIBRARIES}
  )
  
