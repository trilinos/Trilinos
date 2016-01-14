

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTests
  SOURCES
    MultiVector_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTestsXpetraSpecific
  SOURCES
    MultiVector_UnitTestsXpetraSpecific
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


IF (HAVE_XPETRA_THYRA)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultiVector_UnitTestsThyraSpecific
    SOURCES
      MultiVector_UnitTestsThyraSpecific
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()

