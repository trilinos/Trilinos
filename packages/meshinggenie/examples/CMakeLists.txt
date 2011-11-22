INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

<<<<<<< HEAD
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  3dTest
=======
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  2dTest
>>>>>>> 0cf1ba321b730416a92fc0be7153bb97e116788e
  SOURCES main.cpp
#  COMM serial mpi
)

<<<<<<< HEAD
PACKAGE_COPY_FILES_TO_BINARY_DIR(3dTest_files
  SOURCE_FILES checkTest Box.dat README
#  EXEDEPS 3dTest
=======
TRIBITS_COPY_FILES_TO_BINARY_DIR(2dTest_files
  SOURCE_FILES checkTest input.dat README
#  EXEDEPS 2dTest
>>>>>>> 0cf1ba321b730416a92fc0be7153bb97e116788e
)
