

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockCG_hb_test
  SOURCES test_bl_cg_hb.cpp
  ARGS "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_PseudoBlockCG_hb_test
  SOURCES test_pseudo_bl_cg_hb.cpp
  ARGS
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_PseudoBlockStochasticCG_hb_test
  SOURCES test_pseudo_stochastic_cg_hb.cpp
  ARGS
  COMM serial mpi
  )

# Tpetra matrix file reader does not read complex-valued matrices 6/28/2018
#ASSERT_DEFINED(Tpetra_INST_COMPLEX_DOUBLE)
#IF(Tpetra_INST_COMPLEX_DOUBLE)

#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    Tpetra_BlockCG_complex_hb_test
#    SOURCES test_bl_cg_complex_hb.cpp
#    ARGS
#    COMM serial mpi
#    )

#  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockCG_complex_hb_CopyFiles
#    SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
#    SOURCE_FILES mhd1280b.cua
#    EXEDEPS Tpetra_BlockCG_complex_hb_test
#    )

#ENDIF()

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockCG_hb_CopyFiles2
  SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_BlockCG_hb_test
  )

