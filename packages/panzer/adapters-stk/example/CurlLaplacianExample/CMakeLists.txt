
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

#FOREACH( ORDER 1 2 3 4 )
FOREACH( ORDER 1 )
  TRIBITS_ADD_ADVANCED_TEST(
    CurlLaplacianExample-ConvTest-Order-${ORDER}
    TEST_0 EXEC CurlLaplacianExample 
      ARGS --use-epetra --x-elements=5 --y-elements=5 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-${ORDER}-05
    TEST_1 EXEC CurlLaplacianExample 
      ARGS --use-epetra --x-elements=10 --y-elements=10 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-${ORDER}-10
    TEST_2 EXEC CurlLaplacianExample 
      ARGS --use-epetra --x-elements=20 --y-elements=20  --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-${ORDER}-20
    TEST_3 CMND python
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
           ${ORDER}
           MPE-ConvTest-${ORDER}-
           5
           10
           20
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
    )
ENDFOREACH() 
