# -*- cmake -*-

INCLUDE(TribitsAddExecutableAndTest)
#INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Create the MDComm unit test executable
TRIBITS_ADD_EXECUTABLE(
  MDComm_UnitTests
  SOURCES
    MDComm_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

# Create the MDComm sub-comm test executable
TRIBITS_ADD_EXECUTABLE(
  MDComm_SubcommTest
  SOURCES
    MDComm_SubcommTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_1D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_1D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_1D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_1D_6
  COMM mpi
  NUM_MPI_PROCS 6
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_2D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_2D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=-1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_2D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2,-1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_2D_8
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_2D_16
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_3D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=3"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_3D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_3D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_3D_8
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_UnitTests
  NAME MDComm_UnitTests_3D_16
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_1D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_1D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_1D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_1D_6
  COMM mpi
  NUM_MPI_PROCS 6
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_2D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_2D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_2D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_2D_8
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_2D_16
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisCommSizes=4"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_3D_1
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=3"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_3D_2
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_3D_4
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=1,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_3D_8
  COMM mpi
  NUM_MPI_PROCS 8
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=2,2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MDComm_SubcommTest
  NAME MDComm_SubcommTest_3D_16
  COMM mpi
  NUM_MPI_PROCS 16
  ARGS "--teuchos-suppress-startup-banner --numDims=3 --axisCommSizes=4,2"
  STANDARD_PASS_OUTPUT
  )

#TRIBITS_COPY_FILES_TO_BINARY_DIR(target_name
#  SOURCE_FILES source_files
#  EXEDEPS exedeps
#  )
