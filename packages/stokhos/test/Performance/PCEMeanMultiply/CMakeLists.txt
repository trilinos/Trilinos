
ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(CUDA_SOURCES "")
IF (Stokhos_ENABLE_CUDA AND HAVE_TPETRA_CUDA)
  SET(CUDA_SOURCES TestMeanMultiply_Cuda.cpp)
ENDIF()
TRIBITS_ADD_EXECUTABLE(
  KokkosMeanMultiplyTest
  SOURCES TestMeanMultiply.hpp
          TestMeanMultiply.cpp
          TestMeanMultiply_Host.cpp
          ${CUDA_SOURCES}
  COMM serial mpi
  )
