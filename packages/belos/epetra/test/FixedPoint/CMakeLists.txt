
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE(
    fp_hb
    SOURCES test_fp_hb.cpp
    )

  ASSERT_DEFINED(Anasazi_SOURCE_DIR)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestFPFiles
    SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
    SOURCE_FILES bcsstk14.hb
    EXEDEPS fp_hb
  )

  TRIBITS_ADD_TEST(
    fp_hb
    NAME "FixedPoint_epetra_hb_test_no_precond"
    ARGS "--verbose --filename=bcsstk14.hb --no-precond"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_TEST(
    fp_hb
    NAME "FixedPoint_epetra_hb_test_left_precond"
    ARGS "--verbose --filename=bcsstk14.hb --use-precond --left"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_TEST(
    fp_hb
    NAME "FixedPoint_epetra_hb_test_right_precond"
    ARGS "--verbose --filename=bcsstk14.hb --use-precond --right"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
