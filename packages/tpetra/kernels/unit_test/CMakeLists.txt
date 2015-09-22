
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrixUnitTest
  SOURCES
    TestCrsMatrix.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ArithTraits
  SOURCES test-arith-traits.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  blas1_MV_unit
  SOURCES
    blas1_MV
    TpetraKernels_Blas1_MV_UnitTest_Cuda.cpp
    TpetraKernels_Blas1_MV_UnitTest_OpenMP.cpp
    TpetraKernels_Blas1_MV_UnitTest_Serial.cpp
    TpetraKernels_Blas1_MV_UnitTests.hpp
    TpetraKernels_Blas1_MV_UnitTest_Threads.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_quadmath)
IF (${PACKAGE_NAME}_ENABLE_quadmath)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    float128
    SOURCES float128.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF ()
  
