

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tfqmr_hb
    SOURCES test_tfqmr_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1_scaled.hb"
    STANDARD_PASS_OUTPUT 
  )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ptfqmr_hb
    SOURCES test_ptfqmr_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --left-prec"
      "--verbose --right-prec"
    STANDARD_PASS_OUTPUT 
  )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestTFQMRFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockGmres
    SOURCE_FILES orsirr1.hb orsirr1_scaled.hb
    EXEDEPS tfqmr_hb
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
