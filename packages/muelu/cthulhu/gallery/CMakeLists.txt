# MatrixDriver
ADD_SUBDIRECTORY(CMake/tpetra)
ADD_SUBDIRECTORY(CMake/epetra)

# MatrixDriverTests
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(MatrixDriverTests_SOURCES
  MueLu_MatrixDriverTests.cpp
  MueLu_MatrixFactory.hpp
  MueLu_MatrixTypes.hpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixDriverTests
  SOURCES ${MatrixDriverTests_SOURCES}
  )
