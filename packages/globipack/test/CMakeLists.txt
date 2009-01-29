INCLUDE(PackageAddExecutableAndTest)


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

