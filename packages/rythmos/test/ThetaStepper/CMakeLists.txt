
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(../UnitTest)
ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

FILE(GLOB SOURCES *.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ThetaStepper_UnitTest
  SOURCES ${SOURCES}
  DEPLIBS rythmos_test_models
  COMM serial
  STANDARD_PASS_OUTPUT
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(MatrixMatrixCopyFiles
  DEST_FILES modelParams.xml implicitEulerParams.xml
  )
