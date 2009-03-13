
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  levelSchedule
  SOURCES levelSchedule.cpp
  DEPLIBS isorropia_test_utils
  COMM serial
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(
  EXEDEPS levelSchedule
  )
