INCLUDE(TribitsAddExecutableAndTest)   
TRIBITS_ADD_EXECUTABLE_AND_TEST(  CrsMatrix_DefaultMultiply  
  SOURCES 
    CrsMatrix_DefaultMultiply.cpp  
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 

TRIBITS_ADD_EXECUTABLE_AND_TEST(  CrsMatrix_DefaultSolve  
  SOURCES 
    CrsMatrix_DefaultSolve.cpp  
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 

TRIBITS_ADD_EXECUTABLE_AND_TEST( CrsMatrix_UnitTests 
  SOURCES   
    CrsMatrix_UnitTests.cpp 
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 

IF (Kokkos_ENABLE_CUSPARSE) 
  TRIBITS_ADD_EXECUTABLE_AND_TEST( CrsMatrix_CUSPARSE
    SOURCES   
      CrsMatrix_CUSPARSE.cpp 
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  ) 
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST( DefaultSparseOps_TestSparseOps 
  SOURCES   
    DefaultSparseOps_TestSparseOps.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 
