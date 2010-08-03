INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Validator_test
  SOURCES main.cpp
  STANDARD_PASS_OUTPUT
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(ValidatorFileCopies
  DEST_FILES
    testFile.txt testFile2.txt testFile3.txt
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/Validators
  EXEDEPS Validator_test
  )
