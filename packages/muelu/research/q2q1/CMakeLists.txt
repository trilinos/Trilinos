INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters AND ${PACKAGE_NAME}_ENABLE_Teko)
  IF ((${PACKAGE_NAME}_ENABLE_Tpetra                                      AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR
      (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos))

    TRIBITS_ADD_EXECUTABLE(
      Q2Q1
      SOURCES Q2Q1.cpp MueLu_Q2Q1PFactory.hpp
      COMM serial mpi
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(Q2Q1_cp
      SOURCE_FILES driver.xml status.py
      Q2Q1_9x9_A.dat   Q2Q1_9x9_B.dat
      Q2Q1_17x17_A.dat Q2Q1_17x17_B.dat
      BigA9x9.mm  Q2Q1_9x9_A.mm  Q2Q1_9x9_Bt.mm
      p2vMap9x9.mm PresCoords9x9.mm  Q2Q1_9x9_B.mm  VelCoords9x9.mm
      )
  ENDIF()

  #
  # Tests
  #

  INCLUDE(XpetraAddTest)

  IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2 AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters AND ${PACKAGE_NAME}_ENABLE_Teko)
    XPETRA_ADD_TEST_TPETRA(Q2Q1 1)
  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos  AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters AND ${PACKAGE_NAME}_ENABLE_Teko)
    # Disabling Epetra as Vanka is unavailable
    # XPETRA_ADD_TEST_EPETRA(Q2Q1 1)
  ENDIF()

ENDIF()
