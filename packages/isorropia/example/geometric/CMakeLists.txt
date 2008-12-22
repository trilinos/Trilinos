
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  example_rcb
  SOURCES example_rcb.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )


COPY_FILES_TO_BINARY_DIR(ExampleGeometrixCopyFiles
  SOURCE_FILES simple.coords
  EXEDEPS example_rcb
  )
