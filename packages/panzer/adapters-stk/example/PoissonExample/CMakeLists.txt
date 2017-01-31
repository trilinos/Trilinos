
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  PoissonExample
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
  )
  
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  PoissonExample 
#  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
#  COMM serial mpi
#  NUM_MPI_PROCS 4
#  )

TRIBITS_ADD_ADVANCED_TEST(
  PoissonExample-ConvTest
  TEST_0 EXEC PoissonExample 
    ARGS --x-elements=5 --y-elements=5 
    PASS_REGULAR_EXPRESSION "ALL PASSED"
    NUM_MPI_PROCS 4
    OUTPUT_FILE MPE-ConvTest-05
  TEST_1 EXEC PoissonExample 
    ARGS --x-elements=10 --y-elements=10
    PASS_REGULAR_EXPRESSION "ALL PASSED"
    NUM_MPI_PROCS 4
    OUTPUT_FILE MPE-ConvTest-10
  TEST_2 EXEC PoissonExample 
    ARGS --x-elements=20 --y-elements=20 
    PASS_REGULAR_EXPRESSION "ALL PASSED"
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
