INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

SET(XMLFiles 
    blah.txt MissingValidator.xml BadRootElement.xml BadParameterListElement.xml
    NoNameAttribute.xml NoTypeAttribute.xml NoValueAttribute.xml 
    StringToIntegralValidatorBadTag.xml ConflicitingValidatorIDs.xml
	  StringValidatorBadTag.xml CantFindParameterEntryConverter.xml 
    MissingPrototypeValidator.xml)

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

PACKAGE_COPY_FILES_TO_BINARY_DIR(ParameterListFileCopies
  DEST_FILES ${XMLFiles}
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/ParameterList
  EXEDEPS ParameterList_UnitTest
)
