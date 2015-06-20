

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
    
  TRIBITS_ADD_EXECUTABLE(
    RCG_Epetra_Ex_File
    SOURCES RCGEpetraExFile.cpp createEpetraProblem.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleRCGFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/RCG
    SOURCE_FILES bcsstk14.hb
    EXEDEPS RCG_Epetra_Ex_File
    )

ENDIF()
