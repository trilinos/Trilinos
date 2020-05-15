
IF(${PACKAGE_NAME}_ENABLE_Epetra)
  IF (NOT Trilinos_NO_32BIT_GLOBAL_INDICES)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Map
      SOURCES Map.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      CrsMatrix
      SOURCES CrsMatrix.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

   TRIBITS_ADD_EXECUTABLE_AND_TEST(
      VbrMatrix
      SOURCES VbrMatrix.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LinearProblem
      SOURCES LinearProblem.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

  ENDIF()

  IF (NOT Trilinos_NO_64BIT_GLOBAL_INDICES)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Map_LL
      SOURCES Map.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      TARGET_DEFINES -DGALERI_TEST_USE_LONGLONG_GO
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      CrsMatrix_LL
      SOURCES CrsMatrix.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      TARGET_DEFINES -DGALERI_TEST_USE_LONGLONG_GO
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LinearProblem_LL
      SOURCES LinearProblem.cpp
      NUM_MPI_PROCS 1
      COMM serial mpi
      TARGET_DEFINES -DGALERI_TEST_USE_LONGLONG_GO
      )

  ENDIF()
ENDIF()
