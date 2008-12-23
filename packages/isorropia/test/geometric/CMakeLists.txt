
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_geometric
  SOURCES test_geometric.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  ARGS
    "--f=brack_3D.coords"
    "--f=hammond_2D.coords"
    "--f=simple.coords"
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(TestGeometrixCopyFiles
  SOURCE_FILES brack_3D.coords hammond_2D.coords simple.coords
  EXEDEPS test_geometric
  )
