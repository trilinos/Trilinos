
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tfqmr_diag
  SOURCES test_tfqmr_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pseudo_tfqmr_hb
    SOURCES test_pseudo_tfqmr_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1_scaled.hb"
    STANDARD_PASS_OUTPUT 
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tfqmr_hb
    SOURCES test_tfqmr_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1_scaled.hb"
      "--verbose --explicit --filename=orsirr1_scaled.hb"
      "--verbose --recursive --filename=orsirr1_scaled.hb"
    STANDARD_PASS_OUTPUT 
  )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ptfqmr_hb
    SOURCES test_ptfqmr_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --left-prec"
      "--verbose --left-prec --num-rhs=10"
      "--verbose --right-prec"
      "--verbose --right-prec --num-rhs=10"
    STANDARD_PASS_OUTPUT 
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pseudo_ptfqmr_hb
    SOURCES test_pseudo_ptfqmr_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --left-prec"
      "--verbose --left-prec --num-rhs=10"
      "--verbose --right-prec"
      "--verbose --right-prec --num-rhs=10"
    STANDARD_PASS_OUTPUT 
  )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestTFQMRFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES orsirr1.hb orsirr1_scaled.hb
    EXEDEPS tfqmr_hb
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
