
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

## quadrilateral basis order 1-4
FOREACH( ORDER 1 2 3 4 )
  TRIBITS_ADD_ADVANCED_TEST(
    PoissonExample-ConvTest-Quad-Order-${ORDER}
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC PoissonExample
      ARGS --cell="Quad"  --x-elements=5 --y-elements=5 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-p${ORDER}-05
    TEST_1 EXEC PoissonExample
      ARGS --cell="Quad"  --x-elements=10 --y-elements=10 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-p${ORDER}-10
    TEST_2 EXEC PoissonExample
      ARGS --cell="Quad"  --x-elements=20 --y-elements=20 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-p${ORDER}-20
    TEST_3 EXEC PoissonExample
      ARGS --cell="Quad"  --x-elements=40 --y-elements=40 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-p${ORDER}-40
    TEST_4 EXEC PoissonExample
      ARGS --cell="Quad"  --x-elements=80 --y-elements=80 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-p${ORDER}-80
    TEST_5 CMND python
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
         ${ORDER}
         MPE-ConvTest-Quad-p${ORDER}-
         5
         10
         20
         40
         80
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
  )
ENDFOREACH()

## triangle basis order 1-4
FOREACH( ORDER 1 2 3 4 )
  TRIBITS_ADD_ADVANCED_TEST(
    PoissonExample-ConvTest-Tri-Order-${ORDER}
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC PoissonExample
      ARGS --cell="Tri"  --x-elements=5 --y-elements=5 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-p${ORDER}-05
    TEST_1 EXEC PoissonExample
      ARGS --cell="Tri"  --x-elements=10 --y-elements=10 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-p${ORDER}-10
    TEST_2 EXEC PoissonExample
      ARGS --cell="Tri"  --x-elements=20 --y-elements=20 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-p${ORDER}-20
    TEST_3 EXEC PoissonExample
      ARGS --cell="Tri"  --x-elements=40 --y-elements=40 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-p${ORDER}-40
    TEST_4 EXEC PoissonExample
      ARGS --cell="Tri"  --x-elements=80 --y-elements=80 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-p${ORDER}-80
    TEST_5 CMND python
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
         ${ORDER}
         MPE-ConvTest-Tri-p${ORDER}-
         5
         10
         20
         40
         80
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
  )
ENDFOREACH()
