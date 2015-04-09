
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    fp_hb
    SOURCES test_fp_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=bcsstk14.hb"
    STANDARD_PASS_OUTPUT 
    )


  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestFPFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockCG
    SOURCE_FILES bcsstk14.hb
    EXEDEPS fp_hb
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
