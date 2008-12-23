
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  example_rcb
  SOURCES example_rcb.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(ExampleGeometrixCopyFiles
  SOURCE_FILES simple.coords
  EXEDEPS example_rcb
  )
