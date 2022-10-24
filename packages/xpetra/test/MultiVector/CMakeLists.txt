

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTests
  SOURCES
    MultiVector_UnitTests.cpp
    ../Xpetra_UnitTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTestsXpetraSpecific
  SOURCES
    MultiVector_UnitTestsXpetraSpecific.cpp
    ../Xpetra_UnitTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


IF (HAVE_XPETRA_THYRA)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultiVector_UnitTestsThyraSpecific
    SOURCES
      MultiVector_UnitTestsThyraSpecific.cpp
      ../Xpetra_UnitTests.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()
