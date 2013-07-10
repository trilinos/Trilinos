# -*- cmake -*-

INCLUDE(TribitsAddExecutableAndTest)
#INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Unit test the Slice class
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDComm_UnitTests
  SOURCES
    MDComm_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )

#TRIBITS_COPY_FILES_TO_BINARY_DIR(target_name
#  SOURCE_FILES source_files
#  EXEDEPS exedeps
#  )
