
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_tfqmr_diag
  SOURCES test_tfqmr_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_tfqmr_hb
  SOURCES test_tfqmr_hb.cpp
  COMM serial mpi
  ARGS
    "--verbose --not-pseudo --filename=orsirr1_scaled.hb"
    "--verbose --not-pseudo --explicit --filename=orsirr1_scaled.hb"
    "--verbose --not-pseudo --recursive --filename=orsirr1_scaled.hb"
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_pseudo_tfqmr_hb
  SOURCES test_tfqmr_hb.cpp
  COMM serial mpi
  ARGS
    "--verbose --pseudo --filename=orsirr1_scaled.hb"
  STANDARD_PASS_OUTPUT
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyTestTFQMRFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES orsirr1.hb orsirr1_scaled.hb
  EXEDEPS Tpetra_tfqmr_hb
  )
