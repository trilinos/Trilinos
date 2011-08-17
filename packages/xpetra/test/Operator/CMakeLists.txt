
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Operator_UnitTests
  SOURCES
    Operator_UnitTests
    ${PACKAGE_SOURCE_DIR}/../../packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
