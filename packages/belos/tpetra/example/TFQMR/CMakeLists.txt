
TRIBITS_ADD_EXECUTABLE(
  TFQMR_Tpetra_File_Ex
  SOURCES TFQMRTpetraExFile.cpp 
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  Pseudo_Block_TFQMR_Tpetra_File_Ex
  SOURCES PseudoBlockTFQMRTpetraExFile.cpp 
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyExampleTFQMRFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/example/TFQMR
  SOURCE_FILES osrirr1.hb
  EXEDEPS TFQMR_Tpetra_File_Ex Pseudo_Block_TFQMR_Tpetra_File_Ex
)
