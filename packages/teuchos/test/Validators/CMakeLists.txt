INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

SET(SerializationFiles
  blah.txt
  testFile.txt
  testFile2.txt
  testFile3.txt
  MissingValidator.xml
  StringToIntegralValidatorBadTag.xml ConflictingValidatorIDs.xml
  StringValidatorBadTag.xml MissingPrototypeValidator.xml)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Validator_UnitTest
  SOURCES cxx_main.cpp Validator_UnitTest.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Validator_SerializationTest
  SOURCES cxx_main.cpp Validator_SerializationTest.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(ValidatorFileCopies
  DEST_FILES ${SerializationFiles}
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/Validators
  EXEDEPS Validator_UnitTest Validator_SerializationTest
  )
