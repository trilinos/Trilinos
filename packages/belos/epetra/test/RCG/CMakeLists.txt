
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    rcg_hb
    SOURCES test_rcg_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --tol=1e-6 --filename=bcsstk14.hb --num-rhs=10 --max-subspace=100 --recycle=10 --max-iters=4000"
    STANDARD_PASS_OUTPUT 
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestRCGFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/RCG
    SOURCE_FILES bcsstk14.hb gr_30_30.hb
    EXEDEPS rcg_hb
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
    prec_rcg_hb
    SOURCES test_prec_rcg_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --tol=1e-6 --filename=bcsstk14.hb --num-rhs=10 --max-subspace=100 --recycle=5 --max-iters=500 --left-prec"
      "--verbose --tol=1e-6 --filename=bcsstk14.hb --num-rhs=10 --max-subspace=100 --recycle=5 --max-iters=500 --right-prec"
    STANDARD_PASS_OUTPUT 
    )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF()
