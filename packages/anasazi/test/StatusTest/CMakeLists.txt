
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  StatusTest_test
  SOURCES cxx_main.cpp
  COMM serial mpi
  ARGS
    "--verbose"
#    "--debug"
  STANDARD_PASS_OUTPUT
  )
