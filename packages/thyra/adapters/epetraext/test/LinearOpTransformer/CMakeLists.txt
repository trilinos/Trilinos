INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EpetraExtDiagScaledMatProdTransformer_UnitTests
  SOURCES
    EpetraExtDiagScaledMatProdTransformer_UnitTests.cpp
    LocalUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS "--matrix-file=FourByFour.mtx --matrix-file-2=FourByEight.mtx"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EpetraExtAddTransformer_UnitTests
  SOURCES
    EpetraExtAddTransformer_UnitTests.cpp
    LocalUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS "--matrix-file=FourByFour.mtx --matrix-file-2=FourByFour_2.mtx"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EpetraExtDiagScalingTransformer_UnitTests
  SOURCES
    EpetraExtDiagScalingTransformer_UnitTests.cpp
    LocalUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS "--matrix-file=FourByFour.mtx"
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraEpetraExtTransformerCopyFiles_4x4
  DEST_FILES FourByFour.mtx
  SOURCE_DIR ${PARENT_PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic
  SOURCE_PREFIX "In_"
  EXEDEPS EpetraExtDiagScaledMatProdTransformer_UnitTests
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraEpetraExtTransformerCopyFiles_4x8
  DEST_FILES FourByEight.mtx FourByFour_2.mtx
  SOURCE_DIR ${PARENT_PACKAGE_SOURCE_DIR}/adapters/epetraext/test/LinearOpTransformer
  SOURCE_PREFIX "In_"
  EXEDEPS EpetraExtDiagScaledMatProdTransformer_UnitTests EpetraExtAddTransformer_UnitTests
  )
