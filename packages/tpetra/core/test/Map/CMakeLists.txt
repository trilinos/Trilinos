
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_ProblematicLookup
  SOURCES
    Map_ProblematicLookup
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

# kdd Test for bug 2431
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_Bug2431
  SOURCES
    Map_Bug2431
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

# tjf 20 Oct 2016 modified to use default LO, GO types
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_Bug5378
  SOURCES
    Map_Bug5378
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

IF (Tpetra_INST_INT_LONG_LONG OR Tpetra_INST_INT_LONG)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_Bug5399
  SOURCES
    Map_Bug5399
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )
ENDIF()

# This one uses the default LO and GO types.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_Bug6051
  SOURCES
    Map_Bug6051
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_LocalMap
  SOURCES
    Map_LocalMap
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# mfh 26 Sep 2015: Not sure if this one needs GO = int.
# It might be OK just to use the default GO type.
# tjf 20 Oct 2016 modified to use default LO, GO, NT types
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_OneToOne
  SOURCES
    Map_OneToOne
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

IF (Tpetra_INST_INT_LONG_LONG OR Tpetra_INST_INT_LONG)
  # Tpetra bug 5822, first test.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Map_Bug5822
    SOURCES
    Map_Bug5822
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    )

  # Tpetra bug 5822, second test.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Map_Bug5822_2
    SOURCES
    Map_Bug5822_2
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    )
ENDIF()

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
  Map_replaceCommWithSubset
  SOURCES
    Map_replaceCommWithSubset
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
  )

# tjf 20 Oct 2016 modified to use default LO, GO, NT types
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_removeEmptyProcesses
  SOURCES
    Map_removeEmptyProcesses
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
)

# tjf 20 Oct 2016 modified to use default LO, GO types
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_isOneToOne
  SOURCES
    isOneToOne
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
  )

# FIXME (mfh 20 Oct 2015) My recent commit
# 0abee4e8fffdc0975f329bde3b3671db916388bb reverts the fix for Bug
# 6412.  Thus, I am disabling the test for now.  This is a temporary
# measure until we have a chance to extract a test case from the
# application in question and use that to debug the fix.

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Map_Bug6412
#   SOURCES
#     Map_Bug6412.cpp
#     ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM mpi
#   NUM_MPI_PROCS 2-10
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_getMyGlobalIndices
  SOURCES
    Map_getMyGlobalIndices
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
  )

# This test needs to run on at least 2 MPI processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_isLocallyFitted
  SOURCES
    Map_isLocallyFitted.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2-10
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_getMinGlobalIndex
  SOURCES
    Map_getMinGlobalIndex
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_warn_on_post_finalize_destruction_issue_2372
  SOURCES
    Map_warn_on_post_finalize_destruction
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_no_warn_on_pre_finalize_destruction_issue_2372
  SOURCES
    Map_no_warn_on_pre_finalize_destruction
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_initializes_Kokkos_inner_scope
  SOURCES
    Map_initializes_Kokkos_inner_scope
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "FAILED:"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_does_not_double_initialize_Kokkos
  SOURCES
    Map_does_not_double_initialize_Kokkos
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "FAILED:"
  )

# Tpetra bug 5401:
# https://software.sandia.gov/bugzilla/show_bug.cgi?id=5401
#
# This test requires signed GlobalOrdinal type.
# It used to specify GlobalOrdinal=int, but this was unnecessary.
IF (Tpetra_INST_INT_LONG_LONG OR Tpetra_INST_INT_LONG OR Tpetra_INST_INT_INT)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    NegativeBaseIndexTest
    SOURCES
      NegativeBaseIndexTest
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    )
ENDIF ()


