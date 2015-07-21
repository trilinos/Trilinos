
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${PACKAGE_SOURCE_DIR}/LinAlg/test/CrsMatrix)

TRIBITS_ADD_EXECUTABLE_AND_TEST(  CrsMatrix_DefaultMultiply  
  SOURCES 
    CrsMatrix_DefaultMultiply.cpp  
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
) 

TRIBITS_ADD_EXECUTABLE_AND_TEST(  CrsMatrix_DefaultSolve  
  SOURCES 
    CrsMatrix_DefaultSolve.cpp  
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
) 

TRIBITS_ADD_EXECUTABLE_AND_TEST( CrsMatrix_UnitTests 
  SOURCES   
    CrsMatrix_UnitTests.cpp 
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
) 

TRIBITS_ADD_EXECUTABLE_AND_TEST( Raw_AddSparseMatrices
  SOURCES   
    Raw_AddSparseMatrices.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
) 


#disabling for PGI due to this being a known issue with the PGI compiler.
#see bug 4989 for a similar issue.
IF(NOT __LINUX_COMPILER_PGI)
  TRIBITS_ADD_EXECUTABLE_AND_TEST( DefaultSparseOps_TestSparseOps 
    SOURCES   
      DefaultSparseOps_TestSparseOps.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  ) 
ENDIF()

