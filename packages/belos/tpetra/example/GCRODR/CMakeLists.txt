

TRIBITS_ADD_EXECUTABLE(
  GCRODR_Tpetra_File_Ex
  SOURCES GCRODRTpetraExFile.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyExampleGCRODRFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES sherman5.hb
  EXEDEPS GCRODR_Tpetra_File_Ex
  )
