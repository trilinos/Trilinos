INCLUDE(PackageAddExecutableAndTest)

SET(SerializationFiles
  MissingDependee.xml
  MissingDependent.xml
  MissingDependeeTag.xml
  MissingDependentTag.xml
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Dependencies_Serialization_test
  SOURCES cxx_main.cpp Dependencies_SerializationTests.cpp
  DEPLIBS teuchos_xml_pl_test_helpers
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Dependencies_test
  SOURCES cxx_main.cpp Dependencies_UnitTests.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(DependencyFileCopies
  DEST_FILES ${SerializationFiles}
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/Dependencies
  EXEDEPS Dependencies_Serialization_test
)

