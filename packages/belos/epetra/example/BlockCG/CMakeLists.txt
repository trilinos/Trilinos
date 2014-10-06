

ASSERT_DEFINED(Belos_ENABLE_Triutils)
ASSERT_DEFINED(Belos_ENABLE_EpetraExt)
IF (Belos_ENABLE_Triutils AND Belos_ENABLE_EpetraExt)

  TRIBITS_ADD_EXECUTABLE(
    BlockCG_Epetra_File_Ex
    SOURCES BlockCGEpetraExFile.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    PseudoBlockCG_Epetra_File_Ex
    SOURCES PseudoBlockCGEpetraExFile.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleBlockCGFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockCG
    SOURCE_FILES bcsstk14.hb
    EXEDEPS BlockCG_Epetra_File_Ex
    )

  ASSERT_DEFINED(Belos_ENABLE_Ifpack)
  IF(Belos_ENABLE_Ifpack)

    TRIBITS_ADD_EXECUTABLE(
      BlockPrecCG_Epetra_File_Ex
      SOURCES BlockPrecCGEpetraExFile.cpp
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE(
      PseudoBlockPrecCG_Epetra_File_Ex
      SOURCES PseudoBlockPrecCGEpetraExFile.cpp
      COMM serial mpi
      )

  ENDIF(Belos_ENABLE_Ifpack)

ENDIF(Belos_ENABLE_Triutils AND Belos_ENABLE_EpetraExt)
