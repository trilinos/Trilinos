INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  3dTest
  SOURCES main.cpp
#  COMM serial mpi
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(3dTest_files
  SOURCE_FILES checkTest Box.dat README
#  EXEDEPS 3dTest
)
