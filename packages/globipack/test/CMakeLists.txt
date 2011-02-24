INCLUDE(PackageAddExecutableAndTest)


INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_BINARY_DIR}/meritFuncs)
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/meritFuncs)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  VersionUnitTests
  SOURCES
    VersionUnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


ADD_SUBDIRECTORY(meritfuncs)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ArmijoPolyInterpLineSearch_UnitTests
  SOURCES
    ArmijoPolyInterpLineSearch_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS globipack_test_meritfuncs
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  GoldenQuadInterpBracket_UnitTests
  SOURCES
    GoldenQuadInterpBracket_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS globipack_test_meritfuncs
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Brents1DMinimization_UnitTests
  SOURCES
    Brents1DMinimization_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS globipack_test_meritfuncs
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  BrentsLineSearch_UnitTests
  SOURCES
    BrentsLineSearch_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS globipack_test_meritfuncs
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
