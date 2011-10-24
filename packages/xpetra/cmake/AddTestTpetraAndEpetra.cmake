FUNCTION(PACKAGE_ADD_TEST_TPETRA_AND_EPETRA TEST_NAME)

  IF (${PACKAGE_NAME}_ENABLE_Tpetra)
    
    PACKAGE_ADD_TEST(
      ${TEST_NAME}
      NAME ${TEST_NAME}-Tpetra
      ARGS "--linAlgebra=1"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
    
  ENDIF()
  
  IF (${PACKAGE_NAME}_ENABLE_Epetra)
    
    PACKAGE_ADD_TEST(
      ${TEST_NAME}
      NAME ${TEST_NAME}-Epetra
      ARGS "--linAlgebra=0"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
    
  ENDIF()
  
ENDFUNCTION()
