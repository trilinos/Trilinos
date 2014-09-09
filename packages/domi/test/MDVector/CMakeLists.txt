# -*- cmake -*-

INCLUDE(TribitsAddExecutableAndTest)
#INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Create the MDVector unit test executable
TRIBITS_ADD_EXECUTABLE(
  MDVector_UnitTests
  SOURCES
    MDVector_UnitTests.cpp
    ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_1D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_1D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_1D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=4,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_1_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_1_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_2_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_1_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_1_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_2_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_4_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_1_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_2_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=4,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_2_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_4_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_1_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=4,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_2_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=4,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_4_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=4,4,1"
  STANDARD_PASS_OUTPUT
  )

# Create the MDVector comm test executable
TRIBITS_ADD_EXECUTABLE(
  MDVector_CommTests
  SOURCES
    MDVector_CommTests.cpp
    ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_1D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --dims=10"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_1D_1_per
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --dims=10 --periodic=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_1D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --dims=10"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_1D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=10"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_1D_4_per
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=10 --periodic=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_1_1_per
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --periodic=0,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --commDims=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --commDims=2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_2_1_per
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --commDims=2,1 --periodic=1,0"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --commDims=1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --commDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_2_2_per
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --commDims=2 --periodic=0,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --commDims=4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --commDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_2_4_per
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --commDims=2 --periodic=1,0"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --commDims=4,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --commDims=4,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_4_4_per
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --commDims=4,4 --periodic=0,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --dims=8,10,12"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --dims=8,10,12 --commDims=1,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_1_2_per
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --dims=8,10,12 --commDims=1,1,2 --periodic=1,0,0"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --dims=10,12,8 --commDims=1,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_1_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --dims=12,8,10 --commDims=2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_1_1_per
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --dims=12,8,10 --commDims=2,1 --periodic=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=5,7,9 --commDims=1,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=7,9,5 --commDims=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_2_2_per
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=7,9,5 --commDims=1,2 --periodic=0,1,0"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=9,5,7 --commDims=1,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_1_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=9,10,11 --commDims=2,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_1_2_per
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=9,10,11 --commDims=2,1,2 --periodic=0,1,0"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_2_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=10,11,9 --commDims=2,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_1_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=11,9,10 --commDims=4,1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_1_1_per
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=11,9,10 --commDims=4,1,1 --periodic=0,0,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=9,10,8 --commDims=1,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=10,8,9 --commDims=1,4,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_4_2_per
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=10,8,9 --commDims=1,4,2 --periodic=0,0,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_1_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=6,9,12 --commDims=2,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_2_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=9,12,6 --commDims=2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_2_2_per
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=9,12,6 --commDims=2,2 --periodic=1,0,0"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_4_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=12,6,9 --commDims=2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_1_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=9,12,15 --commDims=4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_1_2_per
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=9,12,15 --commDims=4,1 --periodic=1,0,0"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_2_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=12,15,9 --commDims=4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=16,8,12 --commDims=1,4,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_4_4_per
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=16,8,12 --commDims=1,4,4 --periodic=0,1,0"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_2_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=4,4,8 --commDims=2,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_4_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=8,16,8 --commDims=2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_4_2_per
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=8,16,8 --commDims=2,4 --periodic=0,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_1_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=4,4,8 --commDims=4,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_2_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=8,16,8 --commDims=4,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_2_2_per
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=8,16,8 --commDims=4,2,2 --periodic=1,0,0"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_4_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=16,16,4 --commDims=4,4"
  STANDARD_PASS_OUTPUT
  )

# Create the MDVector slice test executable
TRIBITS_ADD_EXECUTABLE(
  MDVector_SliceTests
  SOURCES
    MDVector_SliceTests.cpp
    ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_1D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=1 --commDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_1D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=1 --commDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_1D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=1 --commDims=4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_2D_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_2D_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_2D_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_2D_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_2D_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_2D_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --commDims=4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_1_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_1_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_1_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_2_1_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_1_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_1_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_1_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=1,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_2_1_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_2_2_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=2,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_SliceTests
  NAME MDVector_SliceTests_3D_4_1_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --commDims=4,1,1"
  STANDARD_PASS_OUTPUT
  )
