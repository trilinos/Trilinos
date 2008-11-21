
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  example_rcb
  SOURCES example_rcb.cpp
  COMM serial mpi
  )


SET(INPUT_FILES simple.coords)

COPY_FILES_TO_BINARY_DIR(ExampleGeometrixCopyFiles
  SOURCE_FILES ${INPUT_FILES}
  DEST_FILES ${INPUT_FILES}
  EXEDEPS example_rcb
  )
