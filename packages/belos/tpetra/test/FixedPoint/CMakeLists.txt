

TRIBITS_ADD_EXECUTABLE(
  Tpetra_FixedPoint_hb_test
  SOURCES test_fp_hb.cpp
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_KokkosDense_FixedPoint_hb_test
  SOURCES test_kokkos_fp_hb.cpp
  ARGS
     "--verbose --no-precond"
     "--verbose --use-precond --left"
     "--verbose --use-precond --right"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_FixedPoint_hb_CopyFiles1
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_FixedPoint_hb_test Tpetra_KokkosDense_FixedPoint_hb_test
  )


TRIBITS_ADD_TEST(
  Tpetra_FixedPoint_hb_test
  NAME "FixedPoint_tpetra_hb_test_no_precond"
  ARGS "--verbose --no-precond"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_FixedPoint_hb_test
  NAME "FixedPoint_tpetra_hb_test_left_precond"
  ARGS "--verbose --use-precond --left"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_FixedPoint_hb_test
  NAME "FixedPoint_tpetra_hb_test_right_precond"
  ARGS "--verbose --use-precond --right"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
