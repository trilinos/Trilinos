# -*- cmake -*-

INCLUDE(TribitsAddExecutableAndTest)
#INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Create the MDMap unit test executable
TRIBITS_ADD_EXECUTABLE(
  MDMap_UnitTests
  SOURCES
    MDMap_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_1D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=1 --axisCommSizes=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_1D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=1 --axisCommSizes=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_1D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=1 --axisCommSizes=4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_1_8
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,8"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=4,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_8_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=8,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_1_16
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,16"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_2_8
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2,8"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=4,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_8_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=8,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_2D_16_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=16,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_1_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_1_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_2_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_4_1_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_1_8
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,8"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,4,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_8_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,8,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_1_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_2_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_4_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_4_1_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_4_2_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_8_1_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=8,1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_1_16
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,16"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_2_8
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2,8"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,4,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_8_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,8,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_1_16_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,16,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_1_8
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,1,8"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_2_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_4_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,4,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_2_8_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,8,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_4_1_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_4_2_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_4_4_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_8_1_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=8,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_8_2_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=8,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDMap_UnitTests
  NAME MDMap_UnitTests_3D_16_1_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=16,1,1"
  STANDARD_PASS_OUTPUT
  )
