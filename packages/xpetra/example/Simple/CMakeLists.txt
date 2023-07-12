
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  IF ((Tpetra_INST_INT_INT) OR (Tpetra_INST_INT_LONG_LONG))
    IF ((Tpetra_INST_SERIAL) OR (Tpetra_INST_OPENMP))
      IF (Tpetra_INST_DOUBLE)
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          Laplace1D-Xpetra
          SOURCES Laplace1D-Xpetra.cpp
        )
      ENDIF()
    ENDIF()
  ENDIF()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Laplace1D-Tpetra
    SOURCES Laplace1D-Tpetra.cpp
    )
ENDIF()
