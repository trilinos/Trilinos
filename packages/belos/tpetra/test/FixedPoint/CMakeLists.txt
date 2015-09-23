

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_FixedPoint_hb_test
  SOURCES test_fp_hb.cpp
  ARGS
  COMM serial mpi
  )

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_FixedPoint_hb_CopyFiles1
  SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_FixedPoint_hb_test
  )

