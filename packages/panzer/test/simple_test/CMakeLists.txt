
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
SET(simple_test_SOURCES
  simple_test.cpp
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple_test
  SOURCES ${simple_test_SOURCES}
  COMM serial mpi
  )
