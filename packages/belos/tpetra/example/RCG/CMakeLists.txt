ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)

IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    RCG_Tpetra_Ex_File
    SOURCES RCGTpetraExFile.cpp
    COMM serial mpi
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    Tpetra_CopyExampleRCGFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES bcsstk14.hb
    EXEDEPS RCG_Tpetra_Ex_File
  )
ENDIF()
