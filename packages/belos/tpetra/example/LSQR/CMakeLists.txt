

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)

IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    LSQR_Tpetra_File_Ex
    SOURCES LSQRTpetraExFile.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyExampleLSQRFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockGmres
    SOURCE_FILES orsirr1_scaled.hb
    EXEDEPS LSQR_Tpetra_File_Ex
    )

  # TODO: SOLVE CIRCULAR DEPENDENCY between Belos and Ifpack2
  # ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack2)
  # IF(${PACKAGE_NAME}_ENABLE_Ifpack2)

  # TRIBITS_ADD_EXECUTABLE(
  #     PrecLSQR_Tpetra_File_Ex
  #     SOURCES PrecLSQRTpetraExFile.cpp
  #     COMM serial mpi
  #     )

  # ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack2)

ENDIF(${PACKAGE_NAME}_ENABLE_Tpetra)
