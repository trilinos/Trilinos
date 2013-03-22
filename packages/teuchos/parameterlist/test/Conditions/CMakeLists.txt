INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

SET(SerializationFiles
  MissingParameterEntryDefinition.xml
  MissingValuesTag.xml
  )
  

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Conditions_test
  SOURCES cxx_main.cpp Condition_UnitTests.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Conditions_Serialization_test
  SOURCES cxx_main.cpp Condition_Serialization_UnitTests.cpp
  DEPLIBS teuchos_xml_pl_test_helpers
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ConditionFileCopies
  DEST_FILES ${SerializationFiles}
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/Conditions
  EXEDEPS Conditions_Serialization_test
)

