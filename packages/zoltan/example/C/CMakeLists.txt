INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  simpleBLOCK
  SOURCES simpleBLOCK.c
  ARGS 
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  simpleRCB
  SOURCES simpleRCB.c
  ARGS 
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  simpleGRAPH
  SOURCES simpleGRAPH.c
  ARGS 
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  simplePHG
  SOURCES simplePHG.c
  ARGS 
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(CExamplesCopyFiles
  SOURCE_FILES hypergraph.txt mesh.txt objects.txt graph.txt
 )
