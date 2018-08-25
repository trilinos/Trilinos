
# This test compares Zoltan with Tpetra and Epetra.  Thus, it should run
# only Tpetra and Zoltan with LO=GO=int, as Epetra does.
IF (${PROJECT_NAME}_ENABLE_Epetra)
  IF (Tpetra_INST_INT_INT)
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

#IF (${PROJECT_NAME}_ENABLE_Epetra)
#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    XpetraEpetraMap
#    SOURCES XpetraEpetraMap.cpp
#    COMM serial mpi
#  )
#ENDIF()

IF (${PROJECT_NAME}_ENABLE_Epetra)
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

#KDD 8/23/18  Disabling this test.
#It is not really a Zoltan2 test; it exercises Teuchos::Comm only.
#It hangs on waterman for some reason.  See github #3331
#Teuchos::Comm::createSubcommunicator appears to work in other contexts on
#waterman, so there's no harm in disabling this test.
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  teuchosSubcommTest
#  SOURCES teuchosSubcommTest.cpp
#  NUM_MPI_PROCS 4
#  COMM mpi
#)

# Tests demonstrating Teuchos problems causing issue #3331
TRIBITS_ADD_EXECUTABLE(
  waterman_mpi
  SOURCES waterman_mpi.cpp
  COMM mpi
)

TRIBITS_ADD_TEST(
  waterman_mpi
  NAME waterman_mpi_without_comm_duplicate
  ARGS "X"
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
)

TRIBITS_ADD_TEST(
  waterman_mpi
  NAME waterman_mpi_with_comm_duplicate
  ARGS "Y"
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
)

TRIBITS_ADD_EXECUTABLE(
  waterman_teuchoscomm
  SOURCES waterman_teuchoscomm.cpp
  COMM mpi serial
)

TRIBITS_ADD_TEST(
  waterman_teuchoscomm
  NAME waterman_teuchoscomm_with_comm_duplicate
  ARGS "Y"
  COMM mpi serial
  PASS_REGULAR_EXPRESSION "PASS"
)

TRIBITS_ADD_TEST(
  waterman_teuchoscomm
  NAME waterman_teuchoscomm_without_comm_duplicate
  ARGS "X"
  COMM mpi serial
  PASS_REGULAR_EXPRESSION "PASS"
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
