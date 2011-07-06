INCLUDE(PackageAddExecutableAndTest)   
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINRARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
IF(Kokkos_ENABLE_CUDA)
  CUDA_ADD_LIBRARY(
    unit_cuda
    SOURCES UnitTestCuda.cu
    DEPLIBS kokkosdevicecuda
    )
ELSE()
  PACKAGE_ADD_LIBRARY(
    unit_cuda
    SOURCES UnitTestCudaStub.cpp
    )
ENDIF()
IF(Kokkos_ENABLE_TPI)
  PACKAGE_ADD_LIBRARY(
    unit_tpi
    SOURCES UnitTestTPI.cpp
    DEPLIBS kokkosdevicetpi
  )
ELSE()
    PACKAGE_ADD_LIBRARY(
    unit_tpi
    SOURCES UnitTestTPIStub.cpp
  )
ENDIF()
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  UnitTest
  SOURCES UnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS unit_cuda unit_tpi kokkos kokkosdevicehost kokkosimpl 
  )
  
