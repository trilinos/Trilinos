INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SimpleSolveExample
  SOURCES SimpleSolve
  ARGS
  COMM serial
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SimpleSolveFileExample
  SOURCES SimpleSolve_File
  ARGS
  COMM serial
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MultipleSolvesFileExample
  SOURCES MultipleSolves_File
  ARGS
  COMM serial
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SimpleSolveWithParameters
  SOURCES SimpleSolve_WithParameters
  ARGS
  COMM serial
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_PowerMethod_CopyFiles
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/epetra/example/BlockCG
  SOURCE_FILES bcsstk14.hb
  EXEDEPS SimpleSolveFileExample
  # EXEDEPS Tpetra_BlockCG_hb_test Tpetra_BlockCG_hb_multiprec_test
  )

