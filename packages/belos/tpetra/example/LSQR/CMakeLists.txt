
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

IF (${PACKAGE_NAME}_ENABLE_Triutils)
  TRIBITS_ADD_EXECUTABLE(
    LSQR_Tpetra_File_Ex
    SOURCES LSQRTpetraExFile.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyExampleLSQRFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES orsirr1_scaled.hb
    EXEDEPS LSQR_Tpetra_File_Ex
    )

ENDIF()
