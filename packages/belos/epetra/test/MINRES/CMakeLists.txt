
INCLUDE(TribitsAddExecutableAndTest)

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
    SOURCES test_minres_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=bcsstk14.hb --tol=1e-5"
    STANDARD_PASS_OUTPUT 
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestMinresFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockCG
    SOURCE_FILES bcsstk14.hb
    EXEDEPS minres_hb
    )
 
ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
