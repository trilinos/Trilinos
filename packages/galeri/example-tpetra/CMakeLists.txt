IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Xpetra)
  IF (Tpetra_INST_INT_INT)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Map_Tpetra
      SOURCES Map_Tpetra.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      CrsMatrix_Tpetra
      SOURCES CrsMatrix_Tpetra.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LinearProblem_Tpetra
      SOURCES LinearProblem_Tpetra.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

  ENDIF()

  IF (Tpetra_INST_INT_LONG_LONG)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Map_Tpetra_LL
      SOURCES Map_Tpetra.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      TARGET_DEFINES -DGALERI_TEST_USE_LONGLONG_GO
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      CrsMatrix_Tpetra_LL
      SOURCES CrsMatrix_Tpetra.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      TARGET_DEFINES -DGALERI_TEST_USE_LONGLONG_GO
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LinearProblem_Tpetra_LL
      SOURCES LinearProblem_Tpetra.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      TARGET_DEFINES -DGALERI_TEST_USE_LONGLONG_GO
      )
  ENDIF()
ENDIF()