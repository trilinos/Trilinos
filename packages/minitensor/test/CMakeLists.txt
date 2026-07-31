TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_Containers
  SOURCES test_containers.cc
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_LinearAlgebra
  SOURCES test_linear_algebra.cc
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_MatrixFunctions
  SOURCES test_matrix_functions.cc
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_Rotations
  SOURCES test_rotations.cc
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_GeometryMechanics
  SOURCES test_geometry_mechanics.cc
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_Solvers
  SOURCES test_solvers.cc
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE(
  Perf_Test_01
  SOURCES perf_test_01.cc
  COMM serial mpi
  )

# Each tu/tu_*.cc includes exactly one public header and nothing else,
# so building this executable proves every public header is
# self-contained (independently includable).
FILE(GLOB SELF_CONTAINED_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tu/*.cc)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Self_Contained
  SOURCES ${SELF_CONTAINED_SOURCES}
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )
