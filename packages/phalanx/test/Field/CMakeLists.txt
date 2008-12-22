
INCLUDE(Trilinos_Add_Executable_And_Test)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../Utilities)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  field
  SOURCES Field.cpp
  DEPLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  mdfield
  SOURCES MDField.cpp
  DEPLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )