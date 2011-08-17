INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../utils)
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Xpetra_MatrixDriver
  SOURCES MueLu_MatrixDriver.cpp
  DEPLIBS muelu-xpetra muelu-utils
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Xpetra_MatrixDriverTests
  SOURCES MueLu_MatrixDriverTests.cpp
  DEPLIBS muelu-xpetra muelu-utils
  )
