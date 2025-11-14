ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Experimental)
ASSERT_DEFINED(Tpetra_INST_COMPLEX_DOUBLE)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_gcrodr_hb_test
  SOURCES test_gcrodr_hb.cpp
  ARGS "--verbose --filename=sherman5.hb --tol=1e-4 --num-rhs=2 --max-subspace=61 --recycle=23 --max-cycles=75"
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyTestGCRODRFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES sherman5.hb
    EXEDEPS Tpetra_gcrodr_hb_test
    )

IF (Tpetra_INST_COMPLEX_DOUBLE)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tpetra_gcrodr_complex_hb
      SOURCES test_gcrodr_complex_hb.cpp
      ARGS
        "--verbose --filename=mhd1280b.cua"
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTpetraGCRODRComplexFiles
      SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
      SOURCE_FILES mhd1280b.cua 
      EXEDEPS Tpetra_gcrodr_complex_hb
    )

ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Experimental)
  # At the moment, all the tests in this subdirectory depend on
  # experimental code.

  # Didn't work / failed, because issue in the experimental code.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_BlockGCRODR
    SOURCES test_block_gcrodr.cpp 
    ARGS 
    COMM serial mpi
  )

ENDIF() # ${PACKAGE_NAME}_ENABLE_Experimental
