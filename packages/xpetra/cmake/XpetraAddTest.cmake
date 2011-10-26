# Helper functions for tests using Xpetra::Parameters (--linAlgebra=0/1)

#
# Run test only if Tpetra available
#

FUNCTION(XPETRA_ADD_TEST_TPETRA TEST_NAME COMM NUM_MPI_PROCS)
  IF (${XPETRA_NAME}_ENABLE_Tpetra)
    
    PACKAGE_ADD_TEST(
      ${TEST_NAME}
      NAME ${TEST_NAME}-Tpetra
      ARGS "--linAlgebra=1"
      NUM_MPI_PROCS ${NUM_MPI_PROCS}
      COMM ${COMM}
      )
    
  ENDIF()

ENDFUNCTION()

#
# Run test only if Epetra available
#

FUNCTION(XPETRA_ADD_TEST_EPETRA TEST_NAME COMM NUM_MPI_PROCS)
  IF (${XPETRA_NAME}_ENABLE_Epetra)
    
    PACKAGE_ADD_TEST(
      ${TEST_NAME}
      NAME ${TEST_NAME}-Epetra
      ARGS "--linAlgebra=0"
      NUM_MPI_PROCS ${NUM_MPI_PROCS}
      COMM ${COMM}
      )
    
  ENDIF()

ENDFUNCTION()

#
# Run test for both Epetra and Tpetra (if available)
#

FUNCTION(XPETRA_ADD_TEST_TPETRA_AND_EPETRA TEST_NAME COMM NUM_MPI_PROCS)

  XPETRA_ADD_TEST_TPETRA(${TEST_NAME} ${COMM} ${NUM_MPI_PROCS})
  XPETRA_ADD_TEST_EPETRA(${TEST_NAME} ${COMM} ${NUM_MPI_PROCS})
 
ENDFUNCTION()

#
# Run parallel test for [Epetra,Tpetra] and [1,4] processors
#

FUNCTION(XPETRA_ADD_TEST_PARA_TPETRA_AND_EPETRA TEST_NAME)

  XPETRA_ADD_TEST_TPETRA_AND_EPETRA(${TEST_NAME} "mpi"        4)
  XPETRA_ADD_TEST_TPETRA_AND_EPETRA(${TEST_NAME} "serial mpi" 1)

ENDFUNCTION()

#
# Run sequential test for [Epetra,Tpetra] and 1 processors
#

FUNCTION(XPETRA_ADD_TEST_SEQ_TPETRA_AND_EPETRA TEST_NAME)

  # Also run seq tests with mpi builds, but with only one proc.
  XPETRA_ADD_TEST_TPETRA_AND_EPETRA(${TEST_NAME} "serial mpi" 1)

ENDFUNCTION()
