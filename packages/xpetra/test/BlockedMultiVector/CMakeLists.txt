
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockedMultiVector_UnitTests
  SOURCES
    BlockedMultiVector_UnitTests.cpp
    ../Xpetra_UnitTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_Thyra)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockedMultiVector_UnitTestsThyraSpecific
  SOURCES
    ThyraBlockedMultiVector_UnitTests.cpp
    ../Xpetra_UnitTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
ENDIF()
