
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

## basis order 1-4
FOREACH( ORDER 1 )
  TRIBITS_ADD_ADVANCED_TEST(
    MixedPoissonExample-ConvTest-Order-${ORDER}
    TEST_0 EXEC MixedPoissonExample
      ARGS --x-elements=4 --y-elements=4 --z-elements=4 --hgrad-basis-order=${ORDER} --hdiv-basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-p${ORDER}-04
    TEST_1 EXEC MixedPoissonExample
      ARGS --x-elements=8 --y-elements=8 --z-elements=8 --hgrad-basis-order=${ORDER} --hdiv-basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 8
      OUTPUT_FILE MPE-ConvTest-p${ORDER}-08
    TEST_2 EXEC MixedPoissonExample
      ARGS --x-elements=16 --y-elements=16 --z-elements=16 --hgrad-basis-order=${ORDER} --hdiv-basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 16
      OUTPUT_FILE MPE-ConvTest-p${ORDER}-16
    TEST_3 EXEC MixedPoissonExample
      ARGS --x-elements=32 --y-elements=32 --z-elements=32 --hgrad-basis-order=${ORDER} --hdiv-basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 16
      OUTPUT_FILE MPE-ConvTest-p${ORDER}-32
    TEST_4 CMND python
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
         ${ORDER}
         MPE-ConvTest-p${ORDER}-
         4
         8
         16
         32
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
  )

ENDFOREACH()

