INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES kokkosexample_fixture)

IF( Kokkos_ENABLE_Cuda )

  # add MPI include dir to nvcc
  IF( TPL_ENABLE_MPI )
    CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
  ENDIF()

  TRIBITS_ADD_LIBRARY(
    feint_cuda
    SOURCES feint_cuda.cu
    CUDALIBRARY
    )

  LIST( APPEND LIBRARIES kokkosexample_fixture_cuda )

ENDIF()

TRIBITS_ADD_EXECUTABLE(
  feint
  SOURCES ${SOURCES}
  COMM serial mpi
  DEPLIBS ${LIBRARIES}
  )

