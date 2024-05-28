
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  MixedCurlLaplacianExample
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
  )

IF(PANZER_HAVE_EPETRA_STACK)
  TRIBITS_ADD_ADVANCED_TEST(
    MixedCurlLaplacianExample
    EXCLUDE_IF_NOT_TRUE ${${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS}
    TEST_0 EXEC MixedCurlLaplacianExample
      ARGS --use-epetra --output-filename="base-curl-test-"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
      NUM_MPI_PROCS 4
    TEST_1 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --output-filename="base-curl-test-"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
    COMM mpi
    )
ELSE(PANZER_HAVE_EPETRA_STACK)
  TRIBITS_ADD_ADVANCED_TEST(
    MixedCurlLaplacianExample
    EXCLUDE_IF_NOT_TRUE ${${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS}
    TEST_0 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --output-filename="base-curl-test-"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
    COMM mpi
    )
ENDIF(PANZER_HAVE_EPETRA_STACK)

FOREACH( ORDER 1 2 )
  TRIBITS_ADD_ADVANCED_TEST(
    MixedCurlLaplacianExample-ConvTest-Quad-Order-${ORDER}
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Quad" --x-elements=4 --y-elements=4 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-${ORDER}-04
    TEST_1 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Quad" --x-elements=8 --y-elements=8 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-${ORDER}-08
    TEST_2 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Quad" --x-elements=16 --y-elements=16 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-${ORDER}-16
    TEST_3 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Quad" --x-elements=32 --y-elements=32 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-${ORDER}-32
    TEST_4 CMND python3
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
           ${ORDER}
           MPE-ConvTest-Quad-${ORDER}-
           4
           8
           16
           32
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
    )
ENDFOREACH()

FOREACH( ORDER 3 )
  TRIBITS_ADD_ADVANCED_TEST(
    MixedCurlLaplacianExample-ConvTest-Quad-Order-${ORDER}
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Quad" --x-elements=4 --y-elements=4 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-${ORDER}-04
    TEST_1 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Quad" --x-elements=8 --y-elements=8 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-${ORDER}-08
    TEST_2 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Quad" --x-elements=16 --y-elements=16 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Quad-${ORDER}-16
    TEST_4 CMND python3
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
           ${ORDER}
           MPE-ConvTest-Quad-${ORDER}-
           4
           8
           16
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
    )
ENDFOREACH()

FOREACH( ORDER 1 )
  TRIBITS_ADD_ADVANCED_TEST(
    MixedCurlLaplacianMultiblockExample-ConvTest-Quad-Order-${ORDER}
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --x-blocks=2 --cell="Quad" --x-elements=8 --y-elements=8 --z-elements=4 --basis-order=${ORDER} --output-filename="multiblock-"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-Multiblock-ConvTest-Quad-${ORDER}-08
    TEST_1 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --x-blocks=2 --cell="Quad" --x-elements=16 --y-elements=16 --z-elements=4 --basis-order=${ORDER} --output-filename="multiblock-"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-Multiblock-ConvTest-Quad-${ORDER}-16
    TEST_2 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --x-blocks=2 --cell="Quad" --x-elements=32 --y-elements=32 --z-elements=4 --basis-order=${ORDER} --output-filename="multiblock-"
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-Multiblock-ConvTest-Quad-${ORDER}-32
    TEST_3 CMND python3
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
           ${ORDER}
           MPE-Multiblock-ConvTest-Quad-${ORDER}-
           8
           16
           32
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
    )
ENDFOREACH()

FOREACH( ORDER 1 )
  TRIBITS_ADD_ADVANCED_TEST(
    MixedCurlLaplacianExample-ConvTest-Tri-Order-${ORDER}
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Tri" --x-elements=8 --y-elements=8 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-${ORDER}-08
    TEST_1 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Tri" --x-elements=16 --y-elements=16 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-${ORDER}-16
    TEST_2 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Tri" --x-elements=32 --y-elements=32 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-${ORDER}-32
    TEST_3 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Tri" --x-elements=64 --y-elements=64 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-${ORDER}-64
    TEST_4 CMND python3
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
           ${ORDER}
           MPE-ConvTest-Tri-${ORDER}-
           8
           16
           32
           64
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
    )
ENDFOREACH()

FOREACH( ORDER 2 )
  TRIBITS_ADD_ADVANCED_TEST(
    MixedCurlLaplacianExample-ConvTest-Tri-Order-${ORDER}
    EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
    TEST_0 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Tri" --x-elements=4 --y-elements=4 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-${ORDER}-04
    TEST_1 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Tri" --x-elements=8 --y-elements=8 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-${ORDER}-08
    TEST_2 EXEC MixedCurlLaplacianExample
      ARGS --use-tpetra --use-twod --cell="Tri" --x-elements=16 --y-elements=16 --z-elements=4 --basis-order=${ORDER}
      PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
      NUM_MPI_PROCS 4
      OUTPUT_FILE MPE-ConvTest-Tri-${ORDER}-16
    TEST_3 CMND python3
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
           ${ORDER}
           MPE-ConvTest-Tri-${ORDER}-
           4
           8
           16
      PASS_REGULAR_EXPRESSION "Test Passed"
    COMM mpi
    )
ENDFOREACH()

# FOREACH( ORDER 1 2 3 4 )
#   TRIBITS_ADD_ADVANCED_TEST(
#     CurlLaplacianExample-ConvTest-Hex-Order-${ORDER}
#     EXCLUDE_IF_NOT_TRUE ${PARENT_PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
#     TEST_0 EXEC CurlLaplacianExample
#       ARGS --use-epetra --use-threed --x-elements=4 --y-elements=4 --z-elements=4 --basis-order=${ORDER} --
#       PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
#       NUM_MPI_PROCS 4
#       OUTPUT_FILE MPE-ConvTest-Hex-${ORDER}-04
#     TEST_1 EXEC CurlLaplacianExample
#       ARGS --use-epetra --use-threed --x-elements=8 --y-elements=8 --z-elements=4 --basis-order=${ORDER}
#       PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
#       NUM_MPI_PROCS 4
#       OUTPUT_FILE MPE-ConvTest-Hex-${ORDER}-08
#     TEST_2 EXEC CurlLaplacianExample
#       ARGS --use-epetra --use-threed --x-elements=16 --y-elements=16 --z-elements=4 --basis-order=${ORDER}
#       PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
#       NUM_MPI_PROCS 4
#       OUTPUT_FILE MPE-ConvTest-Hex-${ORDER}-16
#     TEST_3 EXEC CurlLaplacianExample
#       ARGS --use-epetra --use-threed --x-elements=32 --y-elements=32 --z-elements=4 --basis-order=${ORDER}
#       PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
#       NUM_MPI_PROCS 4
#       OUTPUT_FILE MPE-ConvTest-Hex-${ORDER}-32
#     TEST_4 CMND python3
#       ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence_rate.py
#            ${ORDER}
#            MPE-ConvTest-Hex-${ORDER}-
#            4
#            8
#            16
#            32
#       PASS_REGULAR_EXPRESSION "Test Passed"
#     COMM mpi
#     )
# ENDFOREACH()
