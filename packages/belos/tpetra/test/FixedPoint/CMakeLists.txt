BELOS_ADD_EXECUTABLE_AND_TEUCHOS_AND_KOKKOS_TESTS(
  Tpetra_FixedPoint_hb_test
  SOURCES test_fp_hb.cpp
  ARGS
     "--verbose --no-precond"
     "--verbose --use-precond --left"
     "--verbose --use-precond --right"
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_FixedPoint_hb_CopyFiles1
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_FixedPoint_hb_test
  )
