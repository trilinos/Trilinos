
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_pseudo_gmres_multidiag
  SOURCES test_pseudo_gmres_multidiag.cpp BelosLinearMultiShiftProblem.hpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_pseudo_gmres_multi_hb
  SOURCES test_pseudo_gmres_multi_hb.cpp BelosLinearMultiShiftProblem.hpp
  COMM serial mpi
  ARGS "--verbose --filename=orsirr1.hb"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyTestGmresMultiFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES orsirr1.hb
  EXEDEPS Tpetra_pseudo_gmres_multi_hb
)
