TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_UnitTests
  SOURCES
  MatrixMatrix_UnitTests.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS "--matnames-file=\"matrixsystems.xml\" --v"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FECrs_MatrixMatrix_UnitTests
  SOURCES
  FECrs_MatrixMatrix_UnitTests.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )

IF (KokkosKernels_ENABLE_Experimental)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrixKernels_UnitTests
  SOURCES
  MatrixMatrixKernels_UnitTests.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS "--matnames-file=\"matrixsystems.xml\" --v"
  NUM_MPI_PROCS 1
  )
ENDIF()


TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraMatMatFileCopies
DEST_FILES
cdt_m.mtx cdt_tce.mtx cdt_d.mtx
syn1a.mtx syn1b.mtx syn1c.mtx
syn1a_v2.mtx syn1c_v2.mtx
synTa.mtx synTb.mtx synTc.mtx
synBTa.mtx synBTb.mtx synBTc.mtx
synATa.mtx synATb.mtx synATc.mtx
synRecta.mtx synRectb.mtx synRectc.mtx
synadd1a.mtx synadd1b.mtx synadd1c.mtx
smallBTa.mtx smallBTb.mtx smallBTc.mtx
smallAdda.mtx smallAddb.mtx smallAddc.mtx
smallATa.mtx smallATb.mtx smallATc.mtx
rectAPa.mtx rectAPb.mtx rectAPc.mtx
rectAPa_v2.mtx rectAPc_v2.mtx
rectRAa.mtx rectRAb.mtx rectRAc.mtx
rectRPa.mtx rectRPb.mtx rectRPc.mtx
rectQPa.mtx rectQPb.mtx rectQPc.mtx
rectQRa.mtx rectQRb.mtx rectQRc.mtx
rectTransAdda.mtx rectTransAddb.mtx rectTransAddc.mtx
rectAtransIa.mtx rectAtransIb.mtx rectAtransIc.mtx
smallRectAIa.mtx smallRectAIb.mtx smallRectAIc.mtx
denserATa.mtx denserATc.mtx
rapR.mtx rapA.mtx rapP.mtx rapAc.mtx
rapR_domainmap.mtx rapA_domainmap.mtx rapP_domainmap.mtx rapAc_domainmap.mtx
rapR_rangemap.mtx rapA_rangemap.mtx rapP_rangemap.mtx rapAc_rangemap.mtx
rapR_rowmap.mtx rapA_rowmap.mtx rapP_rowmap.mtx rapAc_rowmap.mtx
rapR_colmap.mtx rapA_colmap.mtx rapP_colmap.mtx rapAc_colmap.mtx
SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix/matrices
DEST_DIR ${PACKAGE_BINARY_DIR}/test/MatrixMatrix/matrices
EXEDEPS MatrixMatrix_UnitTests
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraMatMatSystemFileCopy
DEST_FILES
matrixsystems.xml
SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix
EXEDEPS MatrixMatrix_UnitTests
)
