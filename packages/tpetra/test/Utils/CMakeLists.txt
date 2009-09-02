
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_UnitTests
  SOURCES
    TpetraUtils_UnitTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS ${ARGS}
  COMM serial
  STANDARD_PASS_OUTPUT
  )
