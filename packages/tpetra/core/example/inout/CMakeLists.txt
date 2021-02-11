TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NewReaderExample
  SOURCES NewReaderExample.cpp
  ARGS "--file=simple.mtx --distribution=2D "
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  NewReaderExample
  NAME NewReaderExample_rmat
  ARGS "--file=rmat_11.mtx --distribution=2D "
  COMM serial mpi
  )

# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraCore_inout_example_copy_files
  SOURCE_FILES
    simple.mtx
    rmat_11.mtx
    rmat_11_lap.mtx

  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../test/inout
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}	
)


