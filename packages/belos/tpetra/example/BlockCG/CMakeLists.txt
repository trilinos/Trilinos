TRIBITS_ADD_EXECUTABLE(
  BlockCG_Tpetra_File_Ex
  SOURCES BlockCGTpetraExFile.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  PseudoBlockCG_Tpetra_File_Ex
  SOURCES PseudoBlockCGTpetraExFile.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyExampleBlockCGFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS BlockCG_Tpetra_File_Ex
  )
