TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NewReaderExample
  SOURCES NewReaderExample.cpp
  ARGS "--file=simple.mtx --distribution=2D "
  COMM serial mpi
  )

# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraCore_inout_example_copy_files
  SOURCE_FILES
    simple.mtx
)


