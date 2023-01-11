
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  PoissonExampleCurvilinear
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
  )

TRIBITS_ADD_ADVANCED_TEST(
  PoissonExampleCurvilinear-Q2MeshTest
  TEST_0 EXEC PoissonExampleCurvilinear
    ARGS --mesh-filename=annulus_rad_intervals_3_Q2.exo --basis-order=2
    PASS_REGULAR_EXPRESSION "ALL PASSED"
    NUM_MPI_PROCS 2
    OUTPUT_FILE PoissonExampleCurvilinear_3_Q2.out
  COMM mpi
  TEST_1 EXEC PoissonExampleCurvilinear
    ARGS --mesh-filename=annulus_rad_intervals_3_Q1.exo --basis-order=2
    PASS_REGULAR_EXPRESSION "ALL PASSED"
    NUM_MPI_PROCS 2
    OUTPUT_FILE PoissonExampleCurvilinear_3_Q1.out
  COMM mpi
  TEST_2 CMND python
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/compareWithGold.py
         ${CMAKE_CURRENT_SOURCE_DIR}/PoissonExampleCurvilinear_3_Q2.gold
         ${CMAKE_CURRENT_SOURCE_DIR}/PoissonExampleCurvilinear_3_Q1.gold
         PoissonExampleCurvilinear_3_Q2.out
         PoissonExampleCurvilinear_3_Q1.out
    PASS_REGULAR_EXPRESSION "Test Passed"
  COMM mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Curvilinear_mesh_copy
SOURCE_FILES annulus_rad_intervals_3_Q1.exo.2.0 annulus_rad_intervals_3_Q1.exo.2.1 annulus_rad_intervals_3_Q2.exo.2.0 annulus_rad_intervals_3_Q2.exo.2.1
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
EXEDEPS PoissonExampleCurvilinear
)
