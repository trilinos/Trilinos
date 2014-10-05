

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
    
  TRIBITS_ADD_EXECUTABLE(
    GCRODR_Epetra_File_Ex
    SOURCES GCRODREpetraExFile.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleGCRODRFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/GCRODR
    SOURCE_FILES sherman5.hb
    EXEDEPS GCRODR_Epetra_File_Ex
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    TRIBITS_ADD_EXECUTABLE(
      PrecGCRODR_Epetra_File_Ex
      SOURCES PrecGCRODREpetraExFile.cpp createEpetraProblem.cpp
      COMM serial mpi
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
