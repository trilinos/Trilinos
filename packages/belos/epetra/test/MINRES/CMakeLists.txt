
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  minres_diag
  SOURCES test_minres_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  minres_indefinite
  SOURCES test_minres_indefinite.cpp
  COMM serial mpi
  ARGS
    "--debug --verbose --tol=1e-5"
  STANDARD_PASS_OUTPUT 
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    minres_hb
    SOURCES test_minres_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --filename=bcsstk14.hb --tol=1e-5"
      "--verbose --filename=bcsstk14.hb --num-rhs=2 --tol=1e-5"
    STANDARD_PASS_OUTPUT 
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestMinresFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES bcsstk14.hb
    EXEDEPS minres_hb
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      pminres_hb
      SOURCES test_pminres_hb.cpp
      COMM serial mpi
      ARGS
        "--verbose --filename=bcsstk14.hb --left-prec --max-iters=100"
        "--verbose --filename=bcsstk14.hb --right-prec --max-iters=100"
      STANDARD_PASS_OUTPUT
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
