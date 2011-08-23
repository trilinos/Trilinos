
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TwoSquares
  NAME TwoSquares
  SOURCES TwoSquares.cpp
  ARGS 
  COMM serial
  PASS_REGULAR_EXPRESSION "Test passed!"
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(
   TwoSquaresCopyFiles
   DEST_FILES TwoSquares.grid
   SOURCE_DIR ${PACKAGE_SOURCE_DIR}/example/TwoSquares
   EXEDEPS TwoSquares
   )
