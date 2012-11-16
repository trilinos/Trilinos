# Helper functions for tests using Xpetra::Parameters (--linAlgebra=Epetra/Tpetra)

#
# Run test only if Tpetra available
#

FUNCTION(XPETRA_ADD_TEST_TPETRA TEST_NAME NUM_MPI_PROCS)
  IF (${PACKAGE_NAME}_ENABLE_Tpetra)

    TRIBITS_ADD_TEST(
      ${TEST_NAME}
      NAME ${TEST_NAME}-Tpetra
      ARGS "--linAlgebra=Tpetra"
      NUM_MPI_PROCS ${NUM_MPI_PROCS}
      COMM mpi serial
      )

  ENDIF()

ENDFUNCTION()

#
# Run test only if Epetra available
#

FUNCTION(XPETRA_ADD_TEST_EPETRA TEST_NAME NUM_MPI_PROCS)
  IF (${PACKAGE_NAME}_ENABLE_Epetra)

    TRIBITS_ADD_TEST(
      ${TEST_NAME}
      NAME ${TEST_NAME}-Epetra
      ARGS "--linAlgebra=Epetra"
      NUM_MPI_PROCS ${NUM_MPI_PROCS}
      COMM mpi serial
      )

  ENDIF()

ENDFUNCTION()

#
# Run test for both Epetra and Tpetra (if available)
#

FUNCTION(XPETRA_ADD_TEST NUM_MPI_PROCS)

  XPETRA_ADD_TEST_TPETRA(${TEST_NAME} ${NUM_MPI_PROCS})
  XPETRA_ADD_TEST_EPETRA(${TEST_NAME} ${NUM_MPI_PROCS})

ENDFUNCTION()
