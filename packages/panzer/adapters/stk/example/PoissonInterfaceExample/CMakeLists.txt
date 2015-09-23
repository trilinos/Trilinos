INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp)
  
TRIBITS_ADD_EXECUTABLE(
  PoissonInterfaceExample 
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES})

TRIBITS_ADD_ADVANCED_TEST(
  PoissonInterfaceExample_MPI_4
  TEST_0 EXEC PoissonInterfaceExample
    ARGS --nx=5 --nonlinear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
  TEST_1 EXEC PoissonInterfaceExample
    ARGS --nx=5 --linear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
  COMM mpi)

# Could add nx=40 to get to >1.99 instead of just >1.98, but no reason to add
# extra seconds to testing.
TRIBITS_ADD_ADVANCED_TEST(
  PoissonInterfaceExample_MPI_4-ConvTest
  TEST_0 EXEC PoissonInterfaceExample
    ARGS --nx=10 --linear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
    OUTPUT_FILE PIE-ConvTest-10
  TEST_1 EXEC PoissonInterfaceExample
    ARGS --nx=20 --linear
    PASS_REGULAR_EXPRESSION "PASS BASICS"
    NUM_MPI_PROCS 4
    OUTPUT_FILE PIE-ConvTest-20
  TEST_2 CMND python
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../MixedPoissonExample/convergence_rate.py
         PIE-ConvTest-
         10
         20
    PASS_REGULAR_EXPRESSION "Test Passed"
  COMM mpi)
