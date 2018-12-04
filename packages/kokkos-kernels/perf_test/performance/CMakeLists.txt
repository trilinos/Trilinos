INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Add the two performance tests for the Performance class
# located in Kokkos_Performance.hpp.
# performance_validate validates the behavior of the Performance class.
# performance_example is a simple example of using it.
ASSERT_DEFINED(TPL_ENABLE_yaml-cpp)
IF(TPL_ENABLE_yaml-cpp)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    performance_validate
    SOURCES performance_validate.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    performance_example
    SOURCES performance_example.cpp
    COMM serial mpi
    )

ENDIF()
