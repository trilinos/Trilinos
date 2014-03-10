# -*- cmake -*-

INCLUDE(TribitsAddExecutableAndTest)
#INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Create the MDVector unit test executable
TRIBITS_ADD_EXECUTABLE(
  MDVector_UnitTests
  SOURCES
    MDVector_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_1D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=1 --axisCommSizes=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_1D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=1 --axisCommSizes=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_1D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=1 --axisCommSizes=4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=4,1"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_2D_1_8
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=4,2"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_2D_8_1
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=8,1"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_2D_1_16
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,16"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_2D_2_8
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_2D_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=4,4"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_2D_8_2
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=8,2"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_2D_16_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=16,1"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_1_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,1,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_1_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_2_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_1_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,1,1"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_1_1_8
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,4,2"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_1_8_1
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,8,1"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_1_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_2_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_4_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_1_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_2_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,2,1"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_8_1_1
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=8,1,1"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_1_1_16
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,16"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_1_2_8
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_1_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,4,4"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_1_8_2
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,8,2"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_1_16_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,16,1"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_2_1_8
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,1,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_2_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_2_4_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,4,2"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_2_8_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,8,1"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_1_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_2_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_UnitTests
  NAME MDVector_UnitTests_3D_4_4_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,4,1"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_8_1_2
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=8,1,2"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_8_2_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=8,2,1"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_UnitTests
#   NAME MDVector_UnitTests_3D_16_1_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=16,1,1"
#   STANDARD_PASS_OUTPUT
#   )

# Create the MDVector comm test executable
TRIBITS_ADD_EXECUTABLE(
  MDVector_CommTests
  SOURCES
    MDVector_CommTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
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
  NAME MDVector_CommTests_2D_1_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_1_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --axisCommSizes=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --axisCommSizes=2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --axisCommSizes=1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --axisCommSizes=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --axisCommSizes=4,1"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_2D_1_8
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --axisCommSizes=1,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --axisCommSizes=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --axisCommSizes=4,2"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_2D_8_1
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --axisCommSizes=4"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_2D_1_16
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --axisCommSizes=1,16"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_2D_2_8
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --axisCommSizes=2"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_2D_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --axisCommSizes=4,4"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_2D_8_2
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --axisCommSizes=8"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_2D_16_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --axisCommSizes=16"
#   STANDARD_PASS_OUTPUT
#   )

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
  ARGS "--teuchos-suppress-startup-banner --dims=8,10,12 --axisCommSizes=1,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_2_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --dims=10,12,8 --axisCommSizes=1,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_1_1
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --dims=12,8,10 --axisCommSizes=2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_1_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=5,7,9 --axisCommSizes=1,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_2_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=7,9,5 --axisCommSizes=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_4_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=9,5,7 --axisCommSizes=1,4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_1_2
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=9,10,11 --axisCommSizes=2,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_2_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=10,11,9 --axisCommSizes=2,2,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_1_1
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --dims=11,9,10 --axisCommSizes=4,1,1"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_1_1_8
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --dims=8,9,10 --axisCommSizes=1,1,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_2_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=9,10,8 --axisCommSizes=1,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_4_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=10,8,9 --axisCommSizes=1,4,2"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_1_8_1
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --dims=4,16,4 --axisCommSizes=1,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_1_4
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=6,9,12 --axisCommSizes=2,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_2_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=9,12,6 --axisCommSizes=2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_4_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=12,6,9 --axisCommSizes=2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_1_2
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=9,12,15 --axisCommSizes=4,1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_2_1
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --dims=12,15,9 --axisCommSizes=4"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_8_1_1
#   COMM mpi
#   NUM_MPI_PROCS 8
#   ARGS "--teuchos-suppress-startup-banner --dims=15,9,12 --axisCommSizes=8"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_1_1_16
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --dims=8,12,16 --axisCommSizes=1,1,16"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_1_2_8
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --dims=12,16,8 --axisCommSizes=1,2"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_1_4_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=16,8,12 --axisCommSizes=1,4,4"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_1_8_2
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --dims=2,16,4 --axisCommSizes=1,8"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_1_16_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --dims=2,32,2 --axisCommSizes=1,16"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_2_1_8
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --dims=4,2,16 --axisCommSizes=2,1,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_2_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=4,4,8 --axisCommSizes=2,2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_2_4_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=8,16,8 --axisCommSizes=2,4"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_2_8_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --dims=16,16,4 --axisCommSizes=2,8"
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_1_4
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=4,4,8 --axisCommSizes=4,1,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_2_2
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=8,16,8 --axisCommSizes=4,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDVector_CommTests
  NAME MDVector_CommTests_3D_4_4_1
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --dims=16,16,4 --axisCommSizes=4,4"
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_8_1_2
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --dims=8,16,8 --axisCommSizes=8,1,2"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_8_2_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --dims=16,16,4 --axisCommSizes=8,2"
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_TEST(
#   MDVector_CommTests
#   NAME MDVector_CommTests_3D_16_1_1
#   COMM mpi
#   NUM_MPI_PROCS 16
#   ARGS "--teuchos-suppress-startup-banner --dims=16,16,4 --axisCommSizes=16"
#   STANDARD_PASS_OUTPUT
#   )
