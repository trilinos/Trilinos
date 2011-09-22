INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  transpose_test  
  SOURCES transpose_test.cpp
  COMM mpi serial
  )
  
PACKAGE_COPY_FILES_TO_BINARY_DIR(transpose_test_example_copyfiles
  SOURCE_FILES Test.mat
  )