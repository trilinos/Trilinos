INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_LIBRARY(
  teuchos_xml_pl_test_helpers
  HEADERS Teuchos_XMLParameterListTestHelpers.hpp
  SOURCES Teuchos_XMLParameterListTestHelpers.cpp
  TESTONLY
  DEPLIBS teuchos
  )

SET(SerializationFiles 
    BadRootElement.xml BadParameterListElement.xml
    NoNameAttribute.xml NoTypeAttribute.xml NoValueAttribute.xml 
	  CantFindParameterEntryConverter.xml 
    DuplicateParameterIDs.xml )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ParameterList_test 
  SOURCES cxx_main.cpp 
  ARGS --verbose
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ParameterList_UnitTest
  SOURCES ParameterList_UnitTest.cpp Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ParameterEntry_UnitTest
  SOURCES ParameterList_UnitTest.cpp Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ParameterList_SerializationTest
  SOURCES ParameterList_SerializationTest.cpp Teuchos_StandardUnitTestMain.cpp
  DEPLIBS teuchos_xml_pl_test_helpers
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(ParameterListFileCopies
  DEST_FILES ${SerializationFiles}
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/ParameterList
  EXEDEPS ParameterList_SerializationTest
)
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ParameterList_UnitTest_Parallel
  SOURCES ParameterList_UnitTest_Parallel.cpp Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
  ParameterList_UnitTest_Parallel
  NAME_POSTFIX one
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyPLXmlFiles
  DEST_FILES input.xml
  SOURCE_PREFIX "_"
  EXEDEPS ParameterList_UnitTest_Parallel
  )

