INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES kokkoscore_devicehost kokkoscore_devicethreads kokkoscore_impl)

IF( TPL_ENABLE_CUDA AND TPL_ENABLE_CUSPARSE )

  SET(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)

  # add MPI include dir to nvcc
  IF( TPL_ENABLE_MPI )
    CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
  ENDIF()

  TRIBITS_ADD_LIBRARY(
    test_fixture_cuda
    SOURCES TestFixture.cu
    CUDALIBRARY
    )

  LIST( APPEND LIBRARIES
    kokkoscore_devicecuda test_fixture_cuda
  )

ENDIF()

TRIBITS_ADD_EXECUTABLE(
  test_fixture
  SOURCES ${SOURCES}
  COMM serial mpi
  DEPLIBS ${LIBRARIES}
  )

