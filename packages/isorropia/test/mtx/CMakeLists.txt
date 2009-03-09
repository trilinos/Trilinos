
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_balanced_copy
  SOURCES test_create_balanced_copy.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  ARGS
    "--f=simple.mtx"
    "--f=simple_25X50.mtx"
    "--f=simple_50X25.mtx"
    "--f=small.mtx"
    "--f=felix.mtx"
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_simple
  SOURCES test_simple.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  ARGS
    "--f=simple.mtx"
    "--f=simple_25X50.mtx"
    "--f=simple_50X25.mtx"
    "--f=small.mtx"
    "--f=felix.mtx"
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_colorer
  SOURCES test_colorer.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  ARGS
    "--f=simple.mtx"
    "--f=simple_25X50.mtx"
    "--f=simple_50X25.mtx"
    "--f=small.mtx"
    "--f=felix.mtx"
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(TestCreateBalancedCopyCopyFiles
  SOURCE_FILES simple.mtx simple_25X50.mtx simple_50X25.mtx small.mtx felix.mtx
  EXEDEPS test_balanced_copy
  )
