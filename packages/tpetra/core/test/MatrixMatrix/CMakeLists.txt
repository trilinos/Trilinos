TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_UnitTests
  SOURCES MatrixMatrix_UnitTests.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS "--matnames-file=\"matrixsystems.xml\" --v"
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(FECrs_MatrixMatrix_UnitTests SOURCES FECrs_MatrixMatrix_UnitTests.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN})

IF(KokkosKernels_ENABLE_Experimental)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MatrixMatrixKernels_UnitTests
    SOURCES MatrixMatrixKernels_UnitTests.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS "--matnames-file=\"matrixsystems.xml\" --v"
    NUM_MPI_PROCS 1
  )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  TpetraMatMatFileCopies
  DEST_FILES cdt_d.mtx
             cdt_m.mtx
             cdt_tce.mtx
             denserATa.mtx
             denserATc.mtx
             rapA.mtx
             rapA_colmap.mtx
             rapA_domainmap.mtx
             rapA_rangemap.mtx
             rapA_rowmap.mtx
             rapAc.mtx
             rapAc_colmap.mtx
             rapAc_domainmap.mtx
             rapAc_rangemap.mtx
             rapAc_rowmap.mtx
             rapP.mtx
             rapP_colmap.mtx
             rapP_domainmap.mtx
             rapP_rangemap.mtx
             rapP_rowmap.mtx
             rapR.mtx
             rapR_colmap.mtx
             rapR_domainmap.mtx
             rapR_rangemap.mtx
             rapR_rowmap.mtx
             rectAPa.mtx
             rectAPa_v2.mtx
             rectAPb.mtx
             rectAPc.mtx
             rectAPc_v2.mtx
             rectAtransIa.mtx
             rectAtransIb.mtx
             rectAtransIc.mtx
             rectQPa.mtx
             rectQPb.mtx
             rectQPc.mtx
             rectQRa.mtx
             rectQRb.mtx
             rectQRc.mtx
             rectRAa.mtx
             rectRAb.mtx
             rectRAc.mtx
             rectRPa.mtx
             rectRPb.mtx
             rectRPc.mtx
             rectTransAdda.mtx
             rectTransAddb.mtx
             rectTransAddc.mtx
             smallAdda.mtx
             smallAddb.mtx
             smallAddc.mtx
             smallATa.mtx
             smallATb.mtx
             smallATc.mtx
             smallBTa.mtx
             smallBTb.mtx
             smallBTc.mtx
             smallRectAIa.mtx
             smallRectAIb.mtx
             smallRectAIc.mtx
             syn1a.mtx
             syn1a_v2.mtx
             syn1b.mtx
             syn1c.mtx
             syn1c_v2.mtx
             synadd1a.mtx
             synadd1b.mtx
             synadd1c.mtx
             synATa.mtx
             synATb.mtx
             synATc.mtx
             synBTa.mtx
             synBTb.mtx
             synBTc.mtx
             synRecta.mtx
             synRectb.mtx
             synRectc.mtx
             synTa.mtx
             synTb.mtx
             synTc.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix/matrices
  DEST_DIR ${PACKAGE_BINARY_DIR}/test/MatrixMatrix/matrices
  EXEDEPS MatrixMatrix_UnitTests
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  TpetraMatMatSystemFileCopy
  DEST_FILES matrixsystems.xml
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix
  EXEDEPS MatrixMatrix_UnitTests
)
