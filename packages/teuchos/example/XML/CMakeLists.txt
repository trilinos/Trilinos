
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  XmlToParameterList
  SOURCES XmlToParameterList.cpp
  ARGS
    "--xml-in-file=Example1.xml"
    "--xml-in-file=Example1.xml --extra-xml=\"<ParameterList><ParameterList name='Direction'><ParameterList name='Newton'><ParameterList name='Linear Solver'><Parameter name='Tol' type='double' value='1e-2'/></ParameterList></ParameterList></ParameterList></ParameterList>\""
    "--extra-xml=\"<ParameterList/>\""
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_TEST(
  XmlToParameterList
  NAME XmlToParameterList_empty
  ARGS
    "--xml-in-file=Example2.xml"
    "--extra-xml=\" \""
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "EmptyXMLError"
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyExampleXmlFiles
  DEST_FILES Example1.xml Example2.xml
  SOURCE_PREFIX "_"
  EXEDEPS XmlToParameterList
  )
