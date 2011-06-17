INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../utils)
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Cthulhu_MatrixDriver
  SOURCES MueLu_MatrixDriver.cpp
  DEPLIBS muelu-cthulhu muelu-utils
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Cthulhu_MatrixDriverTests
  SOURCES MueLu_MatrixDriverTests.cpp
  DEPLIBS muelu-cthulhu muelu-utils
  )
