TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PowerMethodFileExample
  SOURCES Tpetra_Power_Method_From_File
  ARGS
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_PowerMethod_CopyFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockCG
  SOURCE_FILES bcsstk14.hb
  EXEDEPS PowerMethodFileExample
  # EXEDEPS Tpetra_BlockCG_hb_test Tpetra_BlockCG_hb_multiprec_test
  )

