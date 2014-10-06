ASSERT_DEFINED(Belos_ENABLE_Triutils)
IF (Belos_ENABLE_Triutils)

  ASSERT_DEFINED(Belos_ENABLE_Ifpack)
  IF (Belos_ENABLE_Ifpack)
    TRIBITS_ADD_EXECUTABLE(
      TFQMR_Epetra_File_Ex
      SOURCES TFQMREpetraExFile.cpp createEpetraProblem.cpp
      COMM serial mpi
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleTFQMRFiles
      SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockGmres
      SOURCE_FILES orsirr1.hb
      EXEDEPS TFQMR_Epetra_File_Ex
      )
  ENDIF(Belos_ENABLE_Ifpack)
ENDIF(Belos_ENABLE_Triutils)
