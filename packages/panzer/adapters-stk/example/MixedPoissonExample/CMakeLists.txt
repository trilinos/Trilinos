
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

## basis order 1-3
FOREACH( ORDER 1 2 3)
  SET(HDIV_ORDER ${ORDER})
  MATH(EXPR HGRAD_ORDER "${ORDER}")

  TRIBITS_ADD_ADVANCED_TEST(
    MixedPoissonExample-ConvTest-Hex-Order-${ORDER}
    # TEST_0 EXEC MixedPoissonExample
    #   ARGS --x-elements=2 --y-elements=2 --z-elements=2 --hgrad-basis-order=${HGRAD_ORDER} --hdiv-basis-order=${HDIV_ORDER}
    #   PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
    #   NUM_MPI_PROCS 2
    #   OUTPUT_FILE MPE-ConvTest-Hex-p${ORDER}-02
    TEST_0 EXEC MixedPoissonExample
      ARGS --x-elements=3 --y-elements=3 --z-elements=3 --hgrad-basis-order=${HGRAD_ORDER} --hdiv-basis-order=${HDIV_ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 3
      OUTPUT_FILE MPE-ConvTest-Hex-p${ORDER}-03
    TEST_1 EXEC MixedPoissonExample
      ARGS --x-elements=6 --y-elements=6 --z-elements=6 --hgrad-basis-order=${HGRAD_ORDER} --hdiv-basis-order=${HDIV_ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Hex-p${ORDER}-06
    TEST_2 EXEC MixedPoissonExample
      ARGS --x-elements=12 --y-elements=12 --z-elements=12 --hgrad-basis-order=${HGRAD_ORDER} --hdiv-basis-order=${HDIV_ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Hex-p${ORDER}-12
    TEST_3 CMND python
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
         ${ORDER}
         MPE-ConvTest-Hex-p${ORDER}-
         3
         6
         12
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
  )

ENDFOREACH()

