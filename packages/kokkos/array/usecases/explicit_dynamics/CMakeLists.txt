INCLUDE(TribitsAddExecutable)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../TPL)

SET(SOURCES explicit_main.cpp explicit_test_host.cpp)
SET(LIBRARIES kokkos kokkosdevicehost kokkosimpl)

IF(Kokkos_ENABLE_DeviceCuda)

  add_definitions(-DTEST_KOKKOS_CUDA )

  PACKAGE_ADD_CUDA_LIBRARY(
    lib_cuda_explicit
    SOURCES explicit_test_cuda.cu
    )

  LIST( APPEND LIBRARIES
    kokkosdevicecuda lib_cuda_explicit
  )

ENDIF()

TRIBITS_ADD_EXECUTABLE(
  explicit_dynamics
  SOURCES ${SOURCES}
  COMM serial mpi
  DEPLIBS ${LIBRARIES}
  )

