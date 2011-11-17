
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
SET(simple_example_SOURCES
  simple_example.cpp
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  stk_simple_example
  SOURCES ${simple_example_SOURCES}
  COMM serial mpi
  )
