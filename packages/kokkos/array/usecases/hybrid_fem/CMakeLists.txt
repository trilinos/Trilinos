INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../common)

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp ../common/*.cpp)

SET(LIBRARIES kokkosarray_devicehost kokkosarray_impl)
  

IF( KokkosArray_ENABLE_CUDA AND KokkosArray_ENABLE_CUSPARSE )

  # add MPI include dir to nvcc
  IF( TPL_ENABLE_MPI )
    LIST( APPEND TPL_CUDA_INCLUDE_DIRS ${MPI_BASE_DIR}/include )
  ENDIF()

  PACKAGE_ADD_CUDA_LIBRARY(
    hybrid_fem_cuda
    SOURCES TestCuda.cu
    )

  LIST( APPEND LIBRARIES
    kokkosarray_devicecuda hybrid_fem_cuda ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY}
  )

ENDIF()

TRIBITS_ADD_EXECUTABLE(
  proxyapp
  SOURCES ${SOURCES}
  COMM serial mpi
  DEPLIBS ${LIBRARIES}
  )

