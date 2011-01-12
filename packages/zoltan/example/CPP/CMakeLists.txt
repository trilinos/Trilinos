INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  exampleBLOCK
  SOURCES exampleBLOCK.cpp
  ARGS 
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(CPPExamplesCopyFile
  SOURCE_FILES objects.txt
 )
