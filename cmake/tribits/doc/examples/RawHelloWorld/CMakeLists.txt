# Example of a simple project that uses raw CMake

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4 FATAL_ERROR)

PROJECT(RawHelloWorld)

# Build and install library

SET(HEADERS hello_world_lib.hpp)
SET(SOURCES hello_world_lib.cpp)

ADD_LIBRARY(hello_world_lib ${SOURCES})

INSTALL(TARGETS hello_world_lib DESTINATION lib)
INSTALL(FILES ${HEADERS} DESTINATION include)

# Build and install user executable

ADD_EXECUTABLE(hello_world hello_world_main.cpp)
TARGET_LINK_LIBRARIES(hello_world hello_world_lib)

INSTALL(TARGETS hello_world DESTINATION bin)

# Add some unit tests

ENABLE_TESTING()

# Test the executable
ADD_TEST(test ${CMAKE_CURRENT_BINARY_DIR}/hello_world.exe)
SET_TESTS_PROPERTIES(test PROPERTIES PASS_REGULAR_EXPRESSION "Hello World")

# Build and run some unit tests  
ADD_EXECUTABLE(unit_tests hello_world_unit_tests.cpp)
TARGET_LINK_LIBRARIES(unit_tests hello_world_lib)
ADD_TEST(unit_test ${CMAKE_CURRENT_BINARY_DIR}/unit_tests.exe)
SET_TESTS_PROPERTIES(unit_test
  PROPERTIES PASS_REGULAR_EXPRESSION "All unit tests passed")
