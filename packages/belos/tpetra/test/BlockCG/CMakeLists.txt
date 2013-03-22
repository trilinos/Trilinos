
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockCG_hb_test
  SOURCES test_bl_cg_hb.cpp 
  ARGS 
  COMM serial mpi
  DEPLIBS tpetra
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_PseudoBlockStochasticCG_hb_test
  SOURCES test_pseudo_stochastic_cg_hb.cpp 
  ARGS 
  COMM serial mpi
  DEPLIBS tpetra
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  Tpetra_BlockCG_complex_hb_test
#  SOURCES test_bl_cg_complex_hb.cpp 
#  ARGS 
#  COMM serial mpi
#  )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Tpetra_BlockCG_hb_multiprec_test
#   SOURCES test_bl_cg_hb_multiprec.cpp 
#   ARGS 
#   COMM serial mpi
#   )

#TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockCG_complex_hb_CopyFiles
#  SOURCE_DIR ${Belos_SOURCE_DIR}/test/MVOPTester
#  SOURCE_FILES mhd1280b.cua
#  EXEDEPS Tpetra_BlockCG_complex_hb_test
#  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockCG_hb_CopyFiles1
  SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/test/BlockCG
  SOURCE_FILES bcsstk17.rsa bcsstk18.rsa
  EXEDEPS Tpetra_BlockCG_hb_test 
  # EXEDEPS Tpetra_BlockCG_hb_test Tpetra_BlockCG_hb_multiprec_test
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockCG_hb_CopyFiles2
  SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockCG
  SOURCE_FILES bcsstk14.hb bcsstk15.hb
  EXEDEPS Tpetra_BlockCG_hb_test 
  # EXEDEPS Tpetra_BlockCG_hb_test Tpetra_BlockCG_hb_multiprec_test
  )
