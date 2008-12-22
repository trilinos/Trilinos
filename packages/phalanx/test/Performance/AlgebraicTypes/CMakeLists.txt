
INCLUDE(Trilinos_Add_Executable_And_Test)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  performance_algebraic_types
  SOURCES Performance_AlgebraicTypes.cpp Special_AlgebraicTypes.hpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  performance_algebraic_types_multi
  SOURCES Performance_AlgebraicTypes_Multi.cpp Special_AlgebraicTypes.hpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )