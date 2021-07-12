ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_tfqmr_hb
    SOURCES test_tfqmr_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1_scaled.hb"
      "--verbose --pseudo --filename=orsirr1_scaled.hb"
      "--verbose --explicit --filename=orsirr1_scaled.hb"
      "--verbose --recursive --filename=orsirr1_scaled.hb"
    STANDARD_PASS_OUTPUT 
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyTestTFQMRFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockGmres
    SOURCE_FILES orsirr1.hb orsirr1_scaled.hb
    EXEDEPS Tpetra_tfqmr_hb
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
