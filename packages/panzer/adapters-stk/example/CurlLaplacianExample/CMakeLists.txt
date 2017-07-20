
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp
  )
  
TRIBITS_ADD_EXECUTABLE(
  CurlLaplacianExample 
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
  )

TRIBITS_ADD_ADVANCED_TEST(
  CurlLaplacianExample
  TEST_0 EXEC CurlLaplacianExample
    ARGS --use-epetra
    PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
    NUM_MPI_PROCS 4
  TEST_1 EXEC CurlLaplacianExample
    ARGS --use-tpetra
    PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
    NUM_MPI_PROCS 4
  COMM mpi
  )

FOREACH( ORDER 1 2 3 4 )
  TRIBITS_ADD_ADVANCED_TEST(
    CurlLaplacianExample-ConvTest-Order-${ORDER}
    TEST_0 EXEC CurlLaplacianExample 
      ARGS --use-epetra --x-elements=4 --y-elements=4 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-${ORDER}-04
    TEST_1 EXEC CurlLaplacianExample 
      ARGS --use-epetra --x-elements=8 --y-elements=8 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-${ORDER}-08
    TEST_2 EXEC CurlLaplacianExample 
      ARGS --use-epetra --x-elements=16 --y-elements=16 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-${ORDER}-16
    TEST_3 EXEC CurlLaplacianExample 
      ARGS --use-epetra --x-elements=32 --y-elements=32 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-${ORDER}-32
    TEST_4 CMND python
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
           ${ORDER}
           MPE-ConvTest-${ORDER}-
           4
           8
           16
           32
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
    )
ENDFOREACH() 
