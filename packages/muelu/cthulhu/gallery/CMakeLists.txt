
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
SET(matrix_driver_SOURCES
  MueLu_MatrixDriver.cpp
  MueLu_MatrixFactory.hpp
  MueLu_MatrixTypes.hpp
  )
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  cthulhu_simple_example
  SOURCES ${matrix_driver_SOURCES}
  )
