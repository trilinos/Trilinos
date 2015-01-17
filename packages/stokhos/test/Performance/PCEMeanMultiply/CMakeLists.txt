
ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_ADD_EXECUTABLE(
  KokkosMeanMultiplyTest
  SOURCES TestMeanMultiply.hpp
          TestMeanMultiply.cpp
          TestMeanMultiply_Host.cpp
          TestMeanMultiply_Cuda.cpp
  COMM serial mpi
  )
