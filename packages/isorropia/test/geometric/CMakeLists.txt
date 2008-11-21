
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  test_geometric
  SOURCES test_geometric.cpp
  COMM serial mpi
  ARGS
    "--f=brack_3D.coords"
    "--f=hammond_2D.coords"
    "--f=simple.coords"
  )


SET(INPUT_FILES brack_3D.coords  hammond_2D.coords simple.coords)

COPY_FILES_TO_BINARY_DIR(TestGeometrixCopyFiles
  SOURCE_FILES ${INPUT_FILES}
  DEST_FILES ${INPUT_FILES}
  EXEDEPS test_geometric
  )
