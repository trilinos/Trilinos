

TRIBITS_ADD_EXECUTABLE(
  Tpetra_FixedPoint_hb_test
  SOURCES test_fp_hb.cpp
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_FixedPoint_hb_CopyFiles1
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_FixedPoint_hb_test
  )


TRIBITS_ADD_TEST(
  Tpetra_FixedPoint_hb_test
  NAME "FixedPoint_tpetra_hb_test_no_precond"
  ARGS "--verbose --no-precond"
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  Tpetra_FixedPoint_hb_test
  NAME "FixedPoint_tpetra_hb_test_left_precond"
  ARGS "--verbose --use-precond --left"
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  Tpetra_FixedPoint_hb_test
  NAME "FixedPoint_tpetra_hb_test_right_precond"
  ARGS "--verbose --use-precond --right"
  COMM serial mpi
  )
