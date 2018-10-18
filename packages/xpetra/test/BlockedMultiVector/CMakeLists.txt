
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockedMultiVector_UnitTests
  SOURCES
    BlockedMultiVector_UnitTests
    ../Xpetra_UnitTests
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_Thyra)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockedMultiVector_UnitTestsThyraSpecific
  SOURCES
    ThyraBlockedMultiVector_UnitTests
    ../Xpetra_UnitTests
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
ENDIF()
  