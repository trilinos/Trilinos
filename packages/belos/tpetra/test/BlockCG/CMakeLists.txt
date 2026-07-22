
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_pseudo_cg_indefinite
  SOURCES test_pseudo_cg_indefinite.cpp   
  COMM serial mpi
  ARGS
    "--debug --verbose --tol=1e-5"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockCG_hb_test
  SOURCES test_bl_cg_hb.cpp
  COMM serial mpi
  ARGS
    "--verbose"
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_PseudoBlockCG_hb_test
  SOURCES test_pseudo_bl_cg_hb.cpp
  ARGS "--verbose"
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_PseudoBlockStochasticCG_hb_test
  SOURCES test_pseudo_stochastic_cg_hb.cpp
  ARGS
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_pseudo_cg_hb
  SOURCES test_pseudo_cg_hb.cpp
  COMM serial mpi
  ARGS "--verbose --filename=bcsstk14.hb"
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_resolve_cg_hb
  SOURCES test_resolve_cg_hb.cpp
  COMM serial mpi
  ARGS "--verbose --filename=bcsstk14.hb"
  STANDARD_PASS_OUTPUT
)

ASSERT_DEFINED(Tpetra_INST_COMPLEX_DOUBLE)
IF (Tpetra_INST_COMPLEX_DOUBLE)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_BlockCG_complex_hb_test
    SOURCES test_bl_cg_complex_hb.cpp
    ARGS "--verbose"
    COMM serial mpi
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    Tpetra_BlockCG_complex_hb_CopyFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS Tpetra_BlockCG_complex_hb_test
  )

ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockCG_hb_CopyFiles2
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_BlockCG_hb_test
)
