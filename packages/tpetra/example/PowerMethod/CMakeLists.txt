INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  PowerMethodExample
  SOURCES Tpetra_Power_Method
  ARGS 1000
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  PowerMethodFileExample
  SOURCES Tpetra_Power_Method_From_File
  ARGS
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_PowerMethod_CopyFiles
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/epetra/example/BlockCG
  SOURCE_FILES bcsstk14.hb
  EXEDEPS PowerMethodFileExample
  # EXEDEPS Tpetra_BlockCG_hb_test Tpetra_BlockCG_hb_multiprec_test
  )

