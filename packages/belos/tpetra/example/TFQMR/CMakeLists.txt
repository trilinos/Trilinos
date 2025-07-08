
TRIBITS_ADD_EXECUTABLE(
  TFQMR_Tpetra_File_Ex
  SOURCES TFQMRTpetraExFile.cpp 
  COMM serial mpi
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

TRIBITS_ADD_EXECUTABLE(
  Pseudo_Block_TFQMR_Tpetra_File_Ex
  SOURCES PseudoBlockTFQMRTpetraExFile.cpp 
  COMM serial mpi
)

ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyExampleTFQMRFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES osrirr1.hb
  EXEDEPS TFQMR_Tpetra_File_Ex 
)
