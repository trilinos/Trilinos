
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson02_init_map_vec
  SOURCES lesson02_init_map_vec.cpp
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson02_read_modify_vec
  SOURCES lesson02_read_modify_vec.cpp
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

