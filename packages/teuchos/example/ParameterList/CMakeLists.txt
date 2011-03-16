
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ParameterList_example
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ParameterList_XML_Data_Types_example
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES xml_data_types.cpp
  )
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ParameterList_Custom_Validator_XML
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES custom_validator_xml.cpp
  )
