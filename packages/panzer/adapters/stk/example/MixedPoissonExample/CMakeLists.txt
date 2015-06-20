
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp
  )
  
TRIBITS_ADD_EXECUTABLE(
  MixedPoissonExample 
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
  )

TRIBITS_ADD_ADVANCED_TEST(
  MixedPoissonExample
  TEST_0 EXEC MixedPoissonExample
    ARGS --use-epetra --x-elements=5 --y-elements=5 --z-elements=5
    PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
    NUM_MPI_PROCS 4
  TEST_1 EXEC MixedPoissonExample
    ARGS --use-tpetra --x-elements=5 --y-elements=5 --z-elements=5
    PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
    NUM_MPI_PROCS 4
  COMM mpi
  )

TRIBITS_ADD_ADVANCED_TEST(
  MixedPoissonExample-ConvTest
  TEST_0 EXEC MixedPoissonExample 
    ARGS --use-epetra --x-elements=5 --y-elements=5 --z-elements=5
    PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
    NUM_MPI_PROCS 4
    OUTPUT_FILE MPE-ConvTest-05
  TEST_1 EXEC MixedPoissonExample 
    ARGS --use-epetra --x-elements=10 --y-elements=10 --z-elements=10
    PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
    NUM_MPI_PROCS 4
    OUTPUT_FILE MPE-ConvTest-10
  TEST_2 EXEC MixedPoissonExample 
    ARGS --use-epetra --x-elements=20 --y-elements=20 --z-elements=20
    PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
    NUM_MPI_PROCS 4
    OUTPUT_FILE MPE-ConvTest-20
  TEST_3 CMND python
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
         MPE-ConvTest-
         5
         10
         20
    PASS_REGULAR_EXPRESSION "Test Passed"

  COMM mpi
  )

