
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_ProblematicLookup
  SOURCES
    Map_ProblematicLookup
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_Bug5378
  SOURCES
    Map_Bug5378
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_UnitTests
  SOURCES
    Map_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# Tpetra bug 5401:
# https://software.sandia.gov/bugzilla/show_bug.cgi?id=5401
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NegativeBaseIndexTest
  SOURCES
  NegativeBaseIndexTest
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_OneToOne
  SOURCES
    Map_OneToOne
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

# Tpetra bug 5822, first test.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_Bug5822
  SOURCES
    Map_Bug5822
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

# Tpetra bug 5822, second test.
#
# mfh 23 Jan 2013: This test is disabled for now, because it appears
# to take a lot of memory and may crash the build machine.  Thus, we
# build it, but don't run it by default.
TRIBITS_ADD_EXECUTABLE(
  Map_Bug5822_2
  SOURCES
    Map_Bug5822_2
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

# Tpetra bug 5882
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_ExportTest_Bug5882
  SOURCES
    Map_ExportTest_Bug5882
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "passed!"
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_removeEmptyProcesses
  SOURCES
    Map_removeEmptyProcesses
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
  )
