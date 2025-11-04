TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BiCGStab_diag
  SOURCES test_bicgstab_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BiCGStab_hb
  SOURCES test_bicgstab_hb.cpp
  COMM serial mpi
  ARGS "--verbose --max-iters=10"
       "--verbose --max-iters=10 --tol=1e-10"
       "--verbose --max-iters=10 --num-rhs=2"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyTestBiCGSTABFiles
  SOURCE_DIR ${${PACKAGE_NAME}_SOURCE_DIR}/testmatrices
  SOURCE_FILES cage4.hb
  EXEDEPS Tpetra_BiCGStab_hb
  )
