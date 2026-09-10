
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  PoissonExample
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
  CATEGORIES BASIC PERFORMANCE
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
    TEST_5 CMND python3
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
    TEST_5 CMND python3
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

TRIBITS_ADD_ADVANCED_TEST(
  PoissonExampleCurvilinear-Q2MeshTest
  TEST_0 EXEC PoissonExample
    ARGS --mesh-filename=annulus_rad_intervals_3_Q2.exo --basis-order=2 --problem=annulus
    PASS_REGULAR_EXPRESSION "ALL PASSED"
    NUM_MPI_PROCS 2
    OUTPUT_FILE PoissonExampleCurvilinear_3_Q2.out
  COMM mpi
  TEST_1 EXEC PoissonExample
    ARGS --mesh-filename=annulus_rad_intervals_3_Q1.exo --basis-order=2 --problem=annulus
    PASS_REGULAR_EXPRESSION "ALL PASSED"
    NUM_MPI_PROCS 2
    OUTPUT_FILE PoissonExampleCurvilinear_3_Q1.out
  COMM mpi
  TEST_2 CMND python3
  ARGS   ${CMAKE_CURRENT_SOURCE_DIR}/compareWithGold_Curvilinear.py
         ${CMAKE_CURRENT_SOURCE_DIR}/PoissonExampleCurvilinear_3_Q2.gold
         ${CMAKE_CURRENT_SOURCE_DIR}/PoissonExampleCurvilinear_3_Q1.gold
         PoissonExampleCurvilinear_3_Q2.out
         PoissonExampleCurvilinear_3_Q1.out
    PASS_REGULAR_EXPRESSION "Test Passed"
  COMM mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Curvilinear_Poisson_mesh_copy
SOURCE_FILES annulus_rad_intervals_3_Q1.exo.2.0 annulus_rad_intervals_3_Q1.exo.2.1 annulus_rad_intervals_3_Q2.exo.2.0 annulus_rad_intervals_3_Q2.exo.2.1
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
EXEDEPS PoissonExample
)

# FE assembly (Tpetra::FECrsMatrix) must produce the same answer as the classic path.
# Both runs solve the identical problem; the only difference is whether the Jacobian is
# assembled into a shared FECrsMatrix migrated by endAssembly(), or into separate
# owned/ghosted matrices joined by an explicit export.
TRIBITS_ADD_ADVANCED_TEST(
  PoissonExample-FEAssembly
  TEST_0 EXEC PoissonExample
    ARGS --cell="Quad" --x-elements=10 --y-elements=10 --basis-order=2 --use-fe-assembly
    PASS_REGULAR_EXPRESSION "ALL PASSED"
    NUM_MPI_PROCS 4
    OUTPUT_FILE PoissonExample-FEAssembly-fe.out
  TEST_1 EXEC PoissonExample
    ARGS --cell="Quad" --x-elements=10 --y-elements=10 --basis-order=2 --use-classic-assembly
    PASS_REGULAR_EXPRESSION "ALL PASSED"
    NUM_MPI_PROCS 4
    OUTPUT_FILE PoissonExample-FEAssembly-classic.out
  TEST_2 CMND python3
    ARGS ${PACKAGE_SOURCE_DIR}/example/utilities/compareAssemblyModes.py
         PoissonExample-FEAssembly-fe.out
         PoissonExample-FEAssembly-classic.out
         "L2 Error"
         "H1 Error"
    PASS_REGULAR_EXPRESSION "Test Passed"
  COMM mpi
  )

# Performance tests: a mesh large enough to be worth timing, run through both
# Jacobian assembly paths so the FECrsMatrix path can be compared against the
# classic owned/ghosted path, and so assembly cost can be compared across
# platforms. The assembly is repeated to keep the timing out of the noise; the
# exodus write is off so file I/O does not swamp the measurement. RUN_SERIAL
# keeps the two runs off each other's cores (and off each other's output files).
TRIBITS_ADD_TEST(
  PoissonExample
  NAME PoissonExample-Performance-FEAssembly
  ARGS "--cell=Quad --x-elements=400 --y-elements=400 --basis-order=2 --use-fe-assembly --assembly-loops=10 --no-exodus-output --stacked-timer --test-name=\"Panzer PoissonExample FE Assembly\""
  PASS_REGULAR_EXPRESSION "ALL PASSED"
  COMM mpi
  NUM_MPI_PROCS 4
  CATEGORIES PERFORMANCE
  RUN_SERIAL
  )

TRIBITS_ADD_TEST(
  PoissonExample
  NAME PoissonExample-Performance-ClassicAssembly
  ARGS "--cell=Quad --x-elements=400 --y-elements=400 --basis-order=2 --use-classic-assembly --assembly-loops=10 --no-exodus-output --stacked-timer --test-name=\"Panzer PoissonExample Classic Assembly\""
  PASS_REGULAR_EXPRESSION "ALL PASSED"
  COMM mpi
  NUM_MPI_PROCS 4
  CATEGORIES PERFORMANCE
  RUN_SERIAL
  )
