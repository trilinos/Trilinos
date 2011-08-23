INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  BlockedCrsOperator_UnitTests
  SOURCES
    BlockedCrsOperator_UnitTests
    ${PACKAGE_SOURCE_DIR}/../../packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
  
PACKAGE_COPY_FILES_TO_BINARY_DIR(BlockedCrsOperator_UnitTests_copyfiles
  SOURCE_FILES 
  nsjac_test.mm
  nslhs_test.mm
  nsrhs_test.mm
  )