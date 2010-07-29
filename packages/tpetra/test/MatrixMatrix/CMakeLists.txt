INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
	SimpleMatrixMatrix_test
	SOURCES simplemain.cpp
	COMM serial mpi
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_test
  SOURCES main.cpp
  ARGS --matnames-file="mats.txt"
  COMM serial mpi
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraMatMatFileCopies
DEST_FILES
A.hb B.hb C.hb
SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix
EXEDEPS SimpleMatrixMatrix_test
)

