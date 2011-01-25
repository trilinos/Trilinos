INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
	#SimpleMatrixMatrix_test
	#SOURCES simplemain.cpp
	#COMM serial mpi
#)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_test
  SOURCES main.cpp
  ARGS "--matnames-file=\"matrixsystems.xml\" --v"
  COMM serial mpi
)

#PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraMatMatFileCopies
#DEST_FILES
#A.hb B.hb C.hb
#SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix
#EXEDEPS SimpleMatrixMatrix_test
#)

PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraMatMatFileCopies
DEST_FILES
cdt_d.hb cdt_m.hb cdt_tce.hb
wrong_d.hb wrong_m.hb wrong_tce.hb
gen1a.hb gen1b.hb gen1c.hb
SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix/matrices
DEST_DIR ${PACKAGE_BINARY_DIR}/test/MatrixMatrix/matrices
EXEDEPS MatrixMatrix_test
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraMatMatSystemFileCopy
DEST_FILES
matrixsystems.xml
SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix
EXEDEPS MatrixMatrix_test
)

