
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ml_XML
  SOURCES ml_XML.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Example_ML_XML_CopyFiles
  DEST_FILES
    ml_ParameterList.xml
  )
