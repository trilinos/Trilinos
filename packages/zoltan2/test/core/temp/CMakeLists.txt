
# This test compares Zoltan with Tpetra and Epetra.  Thus, it should run
# only Tpetra and Zoltan with LO=GO=int, as Epetra does.
IF (${PACKAGE_NAME_UC}_HAVE_EPETRA AND Xpetra_ENABLE_Epetra)
  IF (Tpetra_INST_INT_INT AND Tpetra_INST_DOUBLE)
    IF (NOT Zoltan_ENABLE_ULONG_IDS AND NOT Zoltan_ENABLE_ULLONG_IDS)
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

# KDD:  10/29/15
# Disabling this test for until bug 6412 is fixed.
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  mapRemotes
#  SOURCES mapRemotes.cpp
#  COMM mpi
#  FAIL_REGULAR_EXPRESSION "FAIL"
#)

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

#IF (${PACKAGE_NAME_UC}_HAVE_EPETRA)
#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    XpetraEpetraMap
#    SOURCES XpetraEpetraMap.cpp
#    COMM serial mpi
#  )
#ENDIF()

ASSERT_DEFINED(Xpetra_ENABLE_Epetra)
IF (Xpetra_ENABLE_Epetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    XpetraEpetraMatrix
    SOURCES XpetraEpetraMatrix.cpp
    COMM serial mpi
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ddirectoryTest
  SOURCES ddirectoryTest.cpp
  ARGS
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
