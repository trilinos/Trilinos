INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
	#SimpleMatrixMatrix_test
	#SOURCES simplemain.cpp
	#COMM serial mpi
#)

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
  #MatrixMatrix_test
  #SOURCES main.cpp
  #ARGS "--matnames-file=\"matrixsystems.xml\" --v"
  #COMM serial mpi
#)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_UnitTests
  SOURCES 
    MatrixMatrix_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS "--matnames-file=\"matrixsystems.xml\" --v"
  COMM serial
  DEPLIBS tpetraext tpetrainout
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_TwoProcUnitTest
  SOURCES 
    TwoProc_UnitTest.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM mpi
  NUM_MPI_PROCS 2
  DEPLIBS tpetraext
)

IF(NOT TPL_ENABLE_MPI )
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
ENDIF(NOT TPL_ENABLE_MPI )



