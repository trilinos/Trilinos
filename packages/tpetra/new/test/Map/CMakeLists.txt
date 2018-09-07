TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_ProblematicLookup
  SOURCES
    Map2_ProblematicLookup
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_Bug2431
  SOURCES
    Map2_Bug2431
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_Bug5378
  SOURCES
    Map2_Bug5378
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

IF (Tpetra_INST_INT_LONG_LONG OR Tpetra_INST_INT_LONG)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Map2_Bug5399
    SOURCES
      Map2_Bug5399
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    NUM_MPI_PROCS 4
    STANDARD_PASS_OUTPUT
    )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_Bug6051
  SOURCES
    Map2_Bug6051
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  ARGS "--globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_UnitTests
  SOURCES
    Map2_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_LocalMap
  SOURCES
    Map2_LocalMap
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_OneToOne
  SOURCES
    Map2_OneToOne
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

IF (Tpetra_INST_INT_LONG_LONG OR Tpetra_INST_INT_LONG)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Map2_Bug5822
    SOURCES
    Map2_Bug5822
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Map2_Bug5822_2
    SOURCES
    Map2_Bug5822_2
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_replaceCommWithSubset
  SOURCES
    Map2_replaceCommWithSubset
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_removeEmptyProcesses
  SOURCES
    Map2_removeEmptyProcesses
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_isOneToOne
  SOURCES
    Map2_isOneToOne
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
#   Map2_Bug6412
#   SOURCES
#     Map2_Bug6412.cpp
#     ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM mpi
#   NUM_MPI_PROCS 2-10
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_EXECUTABLE(
  Map2_Bug6412
  SOURCES
    Map2_Bug6412.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_getMyGlobalIndices
  SOURCES
    Map2_getMyGlobalIndices
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
  )

# This test needs to run on at least 2 MPI processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_isLocallyFitted
  SOURCES
    Map2_isLocallyFitted.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2-10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_getMinGlobalIndex
  SOURCES
    Map2_getMinGlobalIndex
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_warn_on_post_finalize_destruction_issue_2372
  SOURCES
    Map2_warn_on_post_finalize_destruction
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_no_warn_on_pre_finalize_destruction_issue_2372
  SOURCES
    Map2_no_warn_on_pre_finalize_destruction
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_initializes_Kokkos_inner_scope
  SOURCES
    Map2_initializes_Kokkos_inner_scope
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "FAILED:"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map2_does_not_double_initialize_Kokkos
  SOURCES
    Map2_does_not_double_initialize_Kokkos
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
    Map2_NegativeBaseIndexTest
    SOURCES
      Map2_NegativeBaseIndexTest
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    )
ENDIF ()


