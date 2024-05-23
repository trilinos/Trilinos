
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockGMRES_diag_test
  SOURCES test_bl_gmres_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(Tpetra_INST_FLOAT)

IF (${PACKAGE_NAME}_ENABLE_Triutils)
  IF (Tpetra_INST_FLOAT)

   TRIBITS_ADD_EXECUTABLE_AND_TEST(
     Tpetra_BlockGMRES_hb_DoubleFloat_test
     SOURCES test_bl_gmres_hb_df.cpp
     ARGS "--verbose"
     COMM serial mpi
     )

  ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_resolve_gmres_hb
  SOURCES test_resolve_gmres_hb.cpp
  COMM serial mpi
  ARGS
    "--verbose --filename=orsirr1.hb"
    "--verbose --filename=orsirr1.hb --pseudo"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_pseudo_gmres_hb
  SOURCES test_pseudo_gmres_hb.cpp
  COMM serial mpi
  ARGS "--verbose --filename=orsirr1.hb"
  STANDARD_PASS_OUTPUT
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockGMRES_hb_test
  SOURCES test_bl_gmres_hb.cpp
  ARGS 
    "--verbose --ortho-type=DGKS"
    "--verbose --ortho-type=ICGS"
    "--verbose --ortho-type=IMGS"
  COMM serial mpi
  )

ASSERT_DEFINED(Belos_ENABLE_TSQR)
IF (Belos_ENABLE_TSQR)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockGMRES_hb_TSQR_test
  SOURCES test_bl_gmres_hb.cpp
  ARGS "--verbose --num-rhs=3 --block-size=3 --ortho-type=TSQR"
  COMM serial mpi
  )
ENDIF()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockFGMRES_hb_test
  SOURCES test_bl_fgmres_hb.cpp
  ARGS "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_HybridGMRES_hb_test
  SOURCES test_hybrid_gmres_hb.cpp
  ARGS "--verbose --max-degree=10 --poly-type=Arnoldi"
       "--verbose --max-degree=10 --poly-type=Gmres"
       "--verbose --max-degree=10 --poly-type=Roots"
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockGMRES_hb_CopyFiles
  SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_BlockGMRES_hb_test
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockGMRES_pseudo_CopyFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/example/BlockGmres
  SOURCE_FILES orsirr1.hb
  EXEDEPS Tpetra_BlockGMRES_hb_test
)

