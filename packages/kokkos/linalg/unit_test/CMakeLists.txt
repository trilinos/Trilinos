INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../TPL)

IF(Kokkos_ENABLE_Cuda)
  # mfh 01 Jan 2014: This is magic boilerplate I found in
  # kokkos/core/unit_test/CMakeLists.txt.  I have no idea
  # if it works.
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_COMPILE(OBJECTS_CUDA TestArithTraitsCuda.cu)
ENDIF()


SET(SOURCES
  UnitTestMain.cpp 
  TestCrsMatrix.cpp
  ../../TPL/gtest/gtest-all.cc
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrixUnitTest
  SOURCES ${SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  DEPLIBS kokkoscore kokkoscorelinalg
  )

IF(Kokkos_ENABLE_Cuda)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ArithTraits
    SOURCES test-arith-traits.cpp ${OBJECTS_CUDA}
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS kokkoscore kokkoscorelinalg
    )
ELSE()
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ArithTraits
    SOURCES test-arith-traits.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS kokkoscore kokkoscorelinalg
    )
ENDIF() # Kokkos_ENABLE_Cuda

  
