INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  2dTest
  SOURCES main.cpp
#  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(2dTest_files
  SOURCE_FILES checkTest input.dat README
#  EXEDEPS 2dTest
)
