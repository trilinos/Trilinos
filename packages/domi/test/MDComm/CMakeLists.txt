# -*- cmake -*-

INCLUDE(TribitsAddExecutableAndTest)
#INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Unit test the MDComm class
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDComm_UnitTests_1D_2
  SOURCES
    MDComm_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDComm_UnitTests_1D_6
  SOURCES
    MDComm_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 6
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDComm_UnitTests_2D_4
  SOURCES
    MDComm_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 4
  ARGS "--teuchos-suppress-startup-banner --numDims=2 --axisSizes=2,-1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDComm_SubcommTest_1D_1
  SOURCES
    MDComm_SubcommTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDComm_SubcommTest_1D_2
  SOURCES
    MDComm_SubcommTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--teuchos-suppress-startup-banner --numDims=1"
  STANDARD_PASS_OUTPUT
  )

#TRIBITS_COPY_FILES_TO_BINARY_DIR(target_name
#  SOURCE_FILES source_files
#  EXEDEPS exedeps
#  )
