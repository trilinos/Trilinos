ASSERT_DEFINED(Anasazi_SOURCE_DIR)

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
  ARGS "--verbose"
  COMM serial mpi
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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(Tpetra_INST_COMPLEX_DOUBLE)
ASSERT_DEFINED(Tpetra_INST_FLOAT)

IF (${PACKAGE_NAME}_ENABLE_Triutils)

  IF (Tpetra_INST_FLOAT)
  
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tpetra_BlockCG_DoubleFloat_test
      SOURCES test_bl_cg_hb_df.cpp
      ARGS "--verbose"
      COMM serial mpi
      )

  ENDIF()

  IF (Tpetra_INST_COMPLEX_DOUBLE)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tpetra_BlockCG_complex_hb_test
      SOURCES test_bl_cg_complex_hb.cpp
      ARGS "--verbose"
      COMM serial mpi
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockCG_complex_hb_CopyFiles
      SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
      SOURCE_FILES mhd1280b.cua
      EXEDEPS Tpetra_BlockCG_complex_hb_test
      )

  ENDIF ()

ENDIF()

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockCG_hb_CopyFiles2
  SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_BlockCG_hb_test
  )

