
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  test_geometric
  SOURCES test_geometric.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  ARGS
    "--f=brack_3D.coords"
    "--f=hammond_2D.coords"
    "--f=simple.coords"
  )


COPY_FILES_TO_BINARY_DIR(TestGeometrixCopyFiles
  SOURCE_FILES brack_3D.coords hammond_2D.coords simple.coords
  EXEDEPS test_geometric
  )
