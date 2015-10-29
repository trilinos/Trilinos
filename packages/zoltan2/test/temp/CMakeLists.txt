
IF (${PROJECT_NAME}_ENABLE_Zoltan)
  # mfh 22 Oct 2015: If Tpetra is OFF, it's possible for
  # Tpetra_INST_INT_INT not to be defined at all.  It doesn't hurt to
  # check.  Otherwise, if Tpetra_INST_INT_INT is not defined, asking
  # "IF (Tpetra_INST_INT_INT)" will raise a CMake syntax error.
  IF (DEFINED Tpetra_INST_INT_INT)
    IF (Tpetra_INST_INT_INT)
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        multivectorTest
        SOURCES multivectorTest.cpp
        ARGS
        "401"
        NUM_MPI_PROCS 4
        COMM mpi serial
        )
    ENDIF()
  ENDIF()
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mapRemotes
  SOURCES mapRemotes.cpp
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE(
  testfail
  SOURCES testfail.cpp
  COMM mpi serial
)

TRIBITS_ADD_TEST(
  testfail
  NAME testfail1
  ARGS "FAIL"
  WILL_FAIL 
  FAIL_REGULAR_EXPRESSION "FAIL;BUMMER"
  COMM mpi serial
)

TRIBITS_ADD_TEST(
  testfail
  NAME testfail2
  ARGS "BUMMER"
  WILL_FAIL 
  FAIL_REGULAR_EXPRESSION "FAIL;BUMMER"
  COMM mpi serial
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  teuchosCommTest
  SOURCES teuchosCommTest.cpp
  ARGS
  "4000"
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  teuchosSubcommTest
  SOURCES teuchosSubcommTest.cpp
  ARGS
  "4"
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  absDefineTest
  SOURCES absdefinitiontest.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tpetraOneToOneMapTest
  SOURCES mapOneToOne.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tpetraVecWithCopiesTest
  SOURCES vecWithCopies.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  tpetraMapTest
#  SOURCES mapMemoryLeak.cpp
#  ARGS
#  "4"
#  NUM_MPI_PROCS 4
#  COMM mpi
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#)
