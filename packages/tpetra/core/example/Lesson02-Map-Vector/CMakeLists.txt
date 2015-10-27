
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson02_init_map_vec
  SOURCES lesson02_init_map_vec.cpp
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson02_read_modify_vec_no_kokkos
  SOURCES lesson02_read_modify_vec.cpp
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

ASSERT_DEFINED (Tpetra_INST_DOUBLE)
IF (Tpetra_INST_DOUBLE)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    lesson02_read_modify_vec_kokkos
    SOURCES lesson02_read_modify_vec_kokkos.cpp
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )
ENDIF ()
