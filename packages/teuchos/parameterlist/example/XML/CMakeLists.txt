
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XmlToParameterList
  SOURCES XmlToParameterList.cpp
  ARGS
    "--xml-in-file=Example1.xml"
    "--xml-in-file=Example1.xml --extra-xml-file=ExtraXML1.txt"
    "--extra-xml-file=ExtraXML2.txt"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  XmlToParameterList
  NAME XmlToParameterList_empty
  ARGS
    "--xml-in-file=Example2.xml"
    "--extra-xml-file=EmptyXML.txt"
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "EmptyXMLError"
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleXmlFiles
  DEST_FILES Example1.xml Example2.xml ExtraXML1.txt ExtraXML2.txt EmptyXML.txt
  SOURCE_PREFIX "_"
  EXEDEPS XmlToParameterList
  )
