INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_UnitTests
  SOURCES 
    MatrixMatrix_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS "--matnames-file=\"matrixsystems.xml\" --v"
  DEPLIBS tpetraext tpetrainout
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraMatMatFileCopies
DEST_FILES
cdt_m.mtx cdt_tce.mtx cdt_d.mtx
syn1a.mtx syn1b.mtx syn1c.mtx 
synTa.mtx synTb.mtx synTc.mtx
synBTa.mtx synBTb.mtx synBTc.mtx
synATa.mtx synATb.mtx synATc.mtx
synRecta.mtx synRectb.mtx synRectc.mtx
synadd1a.mtx synadd1b.mtx synadd1c.mtx
smallBTa.mtx smallBTb.mtx smallBTc.mtx
smallAdda.mtx smallAddb.mtx smallAddc.mtx
smallATa.mtx smallATb.mtx smallATc.mtx
rectAPa.mtx rectAPb.mtx rectAPc.mtx
rectRAa.mtx rectRAb.mtx rectRAc.mtx
rectRPa.mtx rectRPb.mtx rectRPc.mtx
rectQPa.mtx rectQPb.mtx rectQPc.mtx
rectQRa.mtx rectQRb.mtx rectQRc.mtx
rectTransAdda.mtx rectTransAddb.mtx rectTransAddc.mtx
rectAtransIa.mtx rectAtransIb.mtx rectAtransIc.mtx
smallRectAIa.mtx smallRectAIb.mtx smallRectAIc.mtx
SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix/matrices
DEST_DIR ${PACKAGE_BINARY_DIR}/test/MatrixMatrix/matrices
EXEDEPS MatrixMatrix_UnitTests
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraMatMatSystemFileCopy
DEST_FILES
matrixsystems.xml
SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix
EXEDEPS MatrixMatrix_UnitTests
)



