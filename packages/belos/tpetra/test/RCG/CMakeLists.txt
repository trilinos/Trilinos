TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_rcg_hb_test
  SOURCES test_rcg_hb.cpp
  ARGS
      "--verbose --tol=1e-6 --filename=bcsstk14.hb --num-rhs=3 --max-subspace=100 --recycle=10 --max-iters=4000"
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyTestRCGFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_rcg_hb_test
  )
