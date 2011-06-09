INCLUDE(PackageAddExecutableAndTest)   


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINRARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(Kokkos_ENABLE_CUDA)
  PACKAGE_ADD_CUDA_LIBRARY(
    perf_cudalib
    SOURCES PerfTestCuda.cu
    DEPLIBS kokkosdevicecuda
    )
ELSE()
  PACKAGE_ADD_LIBRARY(
    perf_cudalib
    SOURCES PerfTestCudaStub.cpp
    )
ENDIF()
IF(Kokkos_ENABLE_TPI)
  PACKAGE_ADD_LIBRARY(
    perf_tpilib
    SOURCES PerfTestTPI.cpp
    DEPLIBS kokkosdevicetpi
    )
ELSE()
  PACKAGE_ADD_LIBRARY(
    perf_tpilib
    SOURCES PerfTestTPIStub.cpp
    )
ENDIF()
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  PerfTest
  SOURCES PerfTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS OUTPUT
  DEPLIBS perf_cudalib perf_tpilib kokkosdevicehost kokkosimpl 
  )
 
