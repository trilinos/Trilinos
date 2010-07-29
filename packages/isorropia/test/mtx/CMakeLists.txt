
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
    "--f=simple_50X25.mtx"
    "--f=small.mtx"
    "--f=felix.mtx"
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_unsym
  SOURCES test_unsym.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  ARGS
    "--f=west0067.mtx"
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_simple_methods
  SOURCES test_simple_methods.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
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
  SOURCE_FILES simple.mtx simple_25X50.mtx simple_50X25.mtx small.mtx felix.mtx arrow100.mtx west0067.mtx
  EXEDEPS test_balanced_copy
  )
