

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
    
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
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb 
  EXEDEPS BlockCG_Epetra_File_Ex
  )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

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

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
