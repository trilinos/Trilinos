INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MultiVector_Adapter_UnitTests
  SOURCES
    Tpetra_MultiVector_Adapter_UnitTests
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Tpetra_CrsMatrix_Adapter_UnitTests
  SOURCES
    Tpetra_CrsMatrix_Adapter_UnitTests
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Epetra_MultiVector_Adapter_UnitTests
    SOURCES
      Epetra_MultiVector_Adapter_UnitTests
      ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Epetra_RowMatrix_Adapter_UnitTests
    SOURCES
      Epetra_RowMatrix_Adapter_UnitTests
      ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  # This test requires EpetraExt for reading Matrix-Market files
  IF (${PACKAGE_NAME}_ENABLE_EpetraExt)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      CrsMatrix_Adapter_Consistency_Tests
      SOURCES
        CrsMatrix_Adapter_Consistency_Tests
        ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      )
  ENDIF()
ENDIF()
