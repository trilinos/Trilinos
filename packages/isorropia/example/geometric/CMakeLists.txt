


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_rcb
  SOURCES example_rcb.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(ExampleGeometrixCopyFiles
  SOURCE_FILES simple.coords
  EXEDEPS example_rcb
  )
