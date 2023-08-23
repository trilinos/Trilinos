
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_tfqmr_diag
  SOURCES test_tfqmr_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_tfqmr_hb
    SOURCES test_tfqmr_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --not-pseudo --filename=orsirr1_scaled.hb"
      "--verbose --not-pseudo --explicit --filename=orsirr1_scaled.hb"
      "--verbose --not-pseudo --recursive --filename=orsirr1_scaled.hb"
    STANDARD_PASS_OUTPUT 
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_pseudo_tfqmr_hb
    SOURCES test_tfqmr_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --pseudo --filename=orsirr1_scaled.hb"
    STANDARD_PASS_OUTPUT 
  )

  #ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack2)
  #IF (${PACKAGE_NAME}_ENABLE_Ifpack2)

  #TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #  Tpetra_ptfqmr_hb
  #  SOURCES test_ptfqmr_hb.cpp 
  #  COMM serial mpi
  #  ARGS
  #    "--verbose --not-pseudo --left-prec"
  #    "--verbose --not-pseudo --left-prec --num-rhs=10"
  #    "--verbose --not-pseudo --right-prec"
  #    "--verbose --not-pseudo --right-prec --num-rhs=10"
  #  STANDARD_PASS_OUTPUT 
  #)

  #TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #  Tpetra_pseudo_ptfqmr_hb
  #  SOURCES test_ptfqmr_hb.cpp 
  #  COMM serial mpi
  #  ARGS
  #    "--verbose --pseudo --left-prec"
  #    "--verbose --pseudo --left-prec --num-rhs=10"
  #    "--verbose --pseudo --right-prec"
  #    "--verbose --pseudo --right-prec --num-rhs=10"
  #  STANDARD_PASS_OUTPUT 
  #)

  #ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack2)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyTestTFQMRFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockGmres
    SOURCE_FILES orsirr1.hb orsirr1_scaled.hb
    EXEDEPS Tpetra_tfqmr_hb
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
