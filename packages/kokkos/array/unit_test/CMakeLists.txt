INCLUDE(PackageAddExecutableAndTest)   
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../TPL)

SET( SOURCES UnitTestMain.cpp UnitTestHost.cpp ../TPL/gtest/gtest.cpp )
SET( LIBRARIES kokkos kokkosdevicehost kokkosimpl )

IF(Kokkos_ENABLE_DeviceTPI)
  LIST( APPEND SOURCES
    UnitTestTPI.cpp
  )
  LIST( APPEND LIBRARIES
    kokkosdevicetpi
  )
ENDIF()

IF(Kokkos_ENABLE_DeviceTBB)
  LIST( APPEND SOURCES
    UnitTestTBB.cpp
  )
  LIST( APPEND LIBRARIES
    kokkosdevicetbb tbb
  )
ENDIF()

IF(Kokkos_ENABLE_DeviceCuda)
  CUDA_ADD_LIBRARY(
    unit_cuda
    UnitTestCuda.cu
    )
  LIST( APPEND SOURCES
    UnitTestCuda.cpp
  )
  LIST( APPEND LIBRARIES
    kokkosdevicecuda cuda cudart
  )
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  UnitTest
  SOURCES ${SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS ${LIBRARIES}
  )
  
