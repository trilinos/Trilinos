
IF (${PACKAGE_NAME}_ENABLE_Epetra)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    VerySimple
    SOURCES cxx_main.cpp
    ARGS -v
    COMM serial
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Xpetra)
  IF (Tpetra_INST_INT_INT)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      VerySimple_Tpetra
      SOURCES cxx_main_tpetra.cpp
      ARGS -v
      COMM serial mpi
      )

  ENDIF()

  IF (Tpetra_INST_INT_LONG_LONG)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      VerySimple_Tpetra_LL
      SOURCES cxx_main_tpetra.cpp
      ARGS -v
      COMM serial mpi
      TARGET_DEFINES -DGALERI_TEST_USE_LONGLONG_GO
      )

  ENDIF()
ENDIF()
