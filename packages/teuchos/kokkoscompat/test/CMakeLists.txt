

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# mfh 25 Jun 2013: Some hacks to get the test to compile.  Not sure if
# these indicate something wrong wth my configuration.  It seems to
# get the Teuchos and Kokkos Array header files just fine, but not the
# Kokkos Compat header files.
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  linkTest
  SOURCES
    linkTest.cpp 
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
