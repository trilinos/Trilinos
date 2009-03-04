INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  EpetraExtDiagScaledMatProdTransformer_UnitTests
  SOURCES
    EpetraExtDiagScaledMatProdTransformer_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --matrix-file=FourByFour.mtx
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(ThyraEpetraExtTransformerCopyFiles
  DEST_FILES FourByFour.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic
  SOURCE_PREFIX "In_"
  EXEDEPS EpetraExtDiagScaledMatProdTransformer_UnitTests
  )
