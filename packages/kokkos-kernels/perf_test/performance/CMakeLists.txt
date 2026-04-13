kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add the two performance tests for the Performance class
# located in Kokkos_Performance.hpp.
# performance_validate validates the behavior of the Performance class.
# performance_example is a simple example of using it.

#don't assert that this is defined anymore
#ASSERT_DEFINED(TPL_ENABLE_yamlcpp)

if(${PACKAGE_NAME}_ENABLE_yamlcpp)
  kokkoskernels_add_unit_test(performance_validate SOURCES performance_validate.cpp)
  kokkoskernels_add_unit_test(performance_example  SOURCES performance_example.cpp)
endif()
