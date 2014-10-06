ASSERT_DEFINED(Belos_ENABLE_Triutils)
ASSERT_DEFINED(Belos_ENABLE_EpetraExt)
IF (Belos_ENABLE_Triutils AND Belos_ENABLE_EpetraExt)
  TRIBITS_ADD_EXECUTABLE(
    LSQR_Epetra_File_Ex
    SOURCES LSQREpetraExFile.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleLSQRFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/LSQR
    SOURCE_FILES orsirr1_scaled.hb
    EXEDEPS LSQR_Epetra_File_Ex
    )

  ASSERT_DEFINED(Belos_ENABLE_Ifpack)
  IF(Belos_ENABLE_Ifpack)
    TRIBITS_ADD_EXECUTABLE(
      PrecLSQR_Epetra_File_Ex
      SOURCES PrecLSQREpetraExFile.cpp
      COMM serial mpi
      )
  ENDIF(Belos_ENABLE_Ifpack)

ENDIF(Belos_ENABLE_Triutils AND Belos_ENABLE_EpetraExt)
