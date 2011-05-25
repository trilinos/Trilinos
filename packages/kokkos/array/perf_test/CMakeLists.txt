INCLUDE(PackageAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINRARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CUDA_SOURCES "")


APPEND_SET(CUDA_SOURCES
  PerfTestMain.cu
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/impl/Kokkos_ArrayBounds.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/impl/Kokkos_ViewTracker.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/impl/Kokkos_MemoryInfo.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/DeviceCuda/Kokkos_DeviceCuda.cu
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/DeviceHost/Kokkos_DeviceHost.cpp
  )
IF(Kokkos_ENABLE_CUDA)
  CUDA_COMPILE(PerfTestMain.cu)
  CUDA_ADD_EXECUTABLE(
    perf_cuda
    ${CUDA_SOURCES}
    )
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  PerfTest
  SOURCES PerfTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS OUTPUT
  DEPLIBS kokkos kokkosdevicetpi
  )


