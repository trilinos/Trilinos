TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_minres_diag
  SOURCES test_minres_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_minres_indefinite
  SOURCES test_minres_indefinite.cpp
  COMM serial mpi
  ARGS
    "--debug --verbose --tol=1e-5"
  STANDARD_PASS_OUTPUT 
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_minres_hb
  SOURCES test_minres_hb.cpp 
  COMM serial mpi
  ARGS
    "--verbose --filename=bcsstk14.hb --tol=1e-5"
    "--verbose --filename=bcsstk14.hb --num-rhs=2 --tol=1e-5"
  STANDARD_PASS_OUTPUT 
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyTestMinresFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_minres_hb
)
