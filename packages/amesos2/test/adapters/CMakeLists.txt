INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MultiVector_Adapter_UnitTests
  SOURCES
    Tpetra_MultiVector_Adapter_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_CrsMatrix_Adapter_UnitTests
  SOURCES
    Tpetra_CrsMatrix_Adapter_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Epetra_MultiVector_Adapter_UnitTests
    SOURCES
      Epetra_MultiVector_Adapter_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Epetra_RowMatrix_Adapter_UnitTests
    SOURCES
      Epetra_RowMatrix_Adapter_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  # This test requires EpetraExt for reading Matrix-Market files
  IF (${PACKAGE_NAME}_ENABLE_EpetraExt)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      CrsMatrix_Adapter_Consistency_Tests
      SOURCES
        CrsMatrix_Adapter_Consistency_Tests
        ${TEUCHOS_STD_UNIT_TEST_MAIN}
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      )
  ENDIF()
ENDIF()
