
INCLUDE(PackageAddExecutableAndTest)


ADD_SUBDIRECTORY(geometric)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  matrix_1
  SOURCES matrix_1.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
PACKAGE_ADD_EXECUTABLE(
  graphedge_weights
  SOURCES graphedge_weights.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
PACKAGE_ADD_EXECUTABLE(
  hgedge_weights
  SOURCES hgedge_weights.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
PACKAGE_ADD_EXECUTABLE(
  part_redist
  SOURCES part_redist.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
PACKAGE_ADD_EXECUTABLE(
  vert_weights
  SOURCES vert_weights.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
PACKAGE_ADD_EXECUTABLE(
  redist_map
  SOURCES redist_map.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
PACKAGE_ADD_EXECUTABLE(
  parallel_matching
  SOURCES parallel_matching.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
