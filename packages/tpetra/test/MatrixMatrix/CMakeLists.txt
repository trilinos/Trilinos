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
  COMM serial mpi
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_TwoProcUnitTest
  SOURCES 
    TwoProc_UnitTest.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM mpi
  NUM_MPI_PROCS 2
)

#PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraMatMatFileCopies
#DEST_FILES
#A.hb B.hb C.hb
#SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix
#EXEDEPS SimpleMatrixMatrix_test
#)

PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraMatMatFileCopies
DEST_FILES
cdt_m.hb cdt_tce.hb
gen1a.hb gen1b.hb
gen2a.hb gen2b.hb
gen3a.hb gen3b.hb
add2x4A.hb add2x4B.hb
bigAddA.hb bigAddB.hb
medAddA.hb medAddB.hb
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



