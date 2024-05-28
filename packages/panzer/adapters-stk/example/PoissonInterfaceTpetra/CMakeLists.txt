TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp)

TRIBITS_ADD_EXECUTABLE(
  PoissonInterfaceTpetra
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES})

TRIBITS_ADD_ADVANCED_TEST(
  PoissonInterfaceTpetra_2d_MPI_4
  EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
  TEST_0 EXEC PoissonInterfaceTpetra
    ARGS --nx=5 --nonlinear --test-jacobian
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
  TEST_1 EXEC PoissonInterfaceTpetra
    ARGS --nx=5 --linear --test-jacobian
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
  COMM mpi)

TRIBITS_ADD_ADVANCED_TEST(
  PoissonInterfaceTpetra_2d_diffsideids_MPI_1
  EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
  TEST_0 EXEC PoissonInterfaceTpetra
    ARGS --test-jacobian --mesh-filename=diffsideids.exo
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 1
  COMM serial mpi)

TRIBITS_ADD_ADVANCED_TEST(
  PoissonInterfaceTpetra_3d_MPI_4
  EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
  TEST_0 EXEC PoissonInterfaceTpetra
    ARGS --3d --nx=3 --nonlinear --test-jacobian
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
  TEST_1 EXEC PoissonInterfaceTpetra
    ARGS --3d --nx=6 --nonlinear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
  TEST_2 EXEC PoissonInterfaceTpetra
    ARGS --3d --nx=7 --linear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
  COMM mpi)

# Could add nx=40 to get to >1.99 instead of just >1.98, but no reason to add
# extra seconds to testing.
TRIBITS_ADD_ADVANCED_TEST(
  PoissonInterfaceTpetra_2d_MPI_4-ConvTest
  EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
  TEST_0 EXEC PoissonInterfaceTpetra
    ARGS --nx=10 --linear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
    OUTPUT_FILE PIE-ConvTest-10
  TEST_1 EXEC PoissonInterfaceTpetra
    ARGS --nx=20 --linear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
    OUTPUT_FILE PIE-ConvTest-20
  TEST_2 CMND python3
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
         PIE-ConvTest-
         10
         20
    PASS_REGULAR_EXPRESSION "Test Passed"
  COMM mpi)

TRIBITS_ADD_ADVANCED_TEST(
  PoissonInterfaceTpetra_3d_MPI_4-ConvTest
  EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
  TEST_0 EXEC PoissonInterfaceTpetra
    ARGS --3d --nx=4 --linear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
    OUTPUT_FILE PIE-ConvTest-4
  TEST_1 EXEC PoissonInterfaceTpetra
    ARGS --3d --nx=8 --linear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
    OUTPUT_FILE PIE-ConvTest-8
  TEST_2 CMND python3
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
         PIE-ConvTest-
         4
         8
    PASS_REGULAR_EXPRESSION "Test Passed"
  COMM mpi)

TRIBITS_COPY_FILES_TO_BINARY_DIR(diffsideids_files_tpetra
  SOURCE_FILES diffsideids.exo
  EXEDEPS PoissonInterfaceTpetra)
