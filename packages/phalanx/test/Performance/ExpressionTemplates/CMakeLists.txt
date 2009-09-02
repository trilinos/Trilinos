
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  performance_expression_templates
  SOURCES Performance_ExpressionTemplates.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )
