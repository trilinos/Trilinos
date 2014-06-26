INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(Kokkos_ENABLE_Cuda)
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( ${KokkosLinAlg_INCLUDE_DIRS} )
  CUDA_INCLUDE_DIRECTORIES( ${KokkosCore_INCLUDE_DIRS} )
  CUDA_INCLUDE_DIRECTORIES( ${Gtest_INCLUDE_DIRS} )
  CUDA_COMPILE(OBJECTS_CUDA TestArithTraitsCuda.cu)
ENDIF()


SET(SOURCES
  UnitTestMain.cpp 
  TestCrsMatrix.cpp
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

  
