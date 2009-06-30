INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MultiVector_ComputeTests
  SOURCES MultiVector_ComputeTests.cpp
          ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkoslinalg
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTests
  SOURCES MultiVector_UnitTests.cpp
          ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkoslinalg
  )
