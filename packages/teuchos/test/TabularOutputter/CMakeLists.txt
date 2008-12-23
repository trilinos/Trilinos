
# Required modules
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TabularOutputterUnitTests
  COMM serial mpi  
  SOURCES
    TabularOutputter_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  )
