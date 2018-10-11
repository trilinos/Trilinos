TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_initializing_tpetra_with_standalone_mpi
  SOURCES initializing_tpetra_with_standalone_mpi.cpp
  ARGS
  COMM mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_initializing_tpetra_with_teuchos_mpi
  SOURCES initializing_tpetra_with_teuchos_mpi.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_initializing_tpetra_with_teuchos_serial
  SOURCES initializing_tpetra_with_teuchos_serial.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_map_contiguous
  SOURCES map_contiguous.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_map_contiguous_and_uniform
  SOURCES map_contiguous_and_uniform.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_map_contiguous_no_global_num
  SOURCES map_contiguous_no_global_num.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_map_cyclic
  SOURCES map_cyclic.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_vector
  SOURCES vector.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_power_method_1
  SOURCES power_method_1.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_matrix_fill_1
  SOURCES matrix_fill_1.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_matrix_construct_heat2d_1
  SOURCES matrix_construct_heat2d_1.cpp
  ARGS
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_matrix_construct_heat2d_2
  SOURCES matrix_construct_heat2d_2.cpp
  ARGS
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  guide_data_redist_1
  SOURCES data_redist_1.cpp
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
)
