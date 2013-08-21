INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../core/TPL)

SET(SOURCES
  TestMain.cpp 
  TestHost.cpp 
  TestThreads.cpp
  ../../core/TPL/gtest/gtest-all.cc
  )

SET(LIBRARIES kokkoscore_devicehost kokkoscore_impl)

IF(Kokkos_ENABLE_OpenMP)
  LIST( APPEND SOURCES
    TestOpenMP.cpp
  )
ENDIF()

IF(Kokkos_ENABLE_Cuda)

  TRIBITS_ADD_LIBRARY(
    kokkos_containers_performance_cuda
    SOURCES TestCuda.cu 
    DEPLIBS
    CUDALIBRARY
  )
  LIST( APPEND SOURCES
    TestCuda.cpp
  )
  LIST( APPEND LIBRARIES
    kokkoscore_devicecuda kokkos_containers_performance_cuda
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PerformanceTest
  SOURCES ${SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "[  PASSED  ]"
  DEPLIBS ${LIBRARIES}
  )
  
