
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
SET(matrix_driver_SOURCES
  matrix_driver.cpp
  matrix_factory.hpp
  matrix_types.hpp
  )
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  simple_example
  SOURCES ${matrix_driver_SOURCES}
  COMM serial
  )
