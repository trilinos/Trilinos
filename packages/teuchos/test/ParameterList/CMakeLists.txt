INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StringIndexedOrderedValueObjectContainer_UnitTests
  SOURCES
    StringIndexedOrderedValueObjectContainer_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_LIBRARY(
  teuchos_xml_pl_test_helpers
  HEADERS Teuchos_XMLParameterListTestHelpers.hpp
  SOURCES Teuchos_XMLParameterListTestHelpers.cpp
  TESTONLY
  DEPLIBS teuchos
  )

SET(PLXMLTestFiles 
    BadRootElement.xml BadParameterListElement.xml
    NoNameAttribute.xml NoTypeAttribute.xml NoValueAttribute.xml 
	  CantFindParameterEntryConverter.xml 
    DuplicateParameterIDs.xml DuplicateSublists.xml )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ParameterList_test 
  SOURCES cxx_main.cpp 
  ARGS --verbose
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ParameterList_UnitTests
  SOURCES ParameterList_UnitTests.cpp
   ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ParameterEntry_UnitTest
  SOURCES ParameterEntry_UnitTest.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ParameterList_SerializationTest
  SOURCES ParameterList_SerializationTest.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  DEPLIBS teuchos_xml_pl_test_helpers
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ParameterList_XMLDuplicateSublists
  SOURCES ParameterList_XMLDuplicateSublists.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  ARGS --filename=DuplicateSublists.xml
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ParameterListFileCopies
  DEST_FILES ${PLXMLTestFiles}
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/ParameterList
  EXEDEPS ParameterList_SerializationTest
)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ParameterList_UnitTest_Parallel
  SOURCES ParameterList_UnitTest_Parallel.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  ParameterList_UnitTest_Parallel
  NAME_POSTFIX one
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyPLXmlFiles
  DEST_FILES input.xml
  SOURCE_PREFIX "_"
  EXEDEPS ParameterList_UnitTest_Parallel
  )

