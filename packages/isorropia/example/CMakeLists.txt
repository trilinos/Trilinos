


ADD_SUBDIRECTORY(geometric)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  matrix_1
  SOURCES matrix_1.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
TRIBITS_ADD_EXECUTABLE(
  graphedge_weights
  SOURCES graphedge_weights.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
TRIBITS_ADD_EXECUTABLE(
  hgedge_weights
  SOURCES hgedge_weights.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
TRIBITS_ADD_EXECUTABLE(
  part_redist
  SOURCES part_redist.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
TRIBITS_ADD_EXECUTABLE(
  vert_weights
  SOURCES vert_weights.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
TRIBITS_ADD_EXECUTABLE(
  redist_map
  SOURCES redist_map.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
TRIBITS_ADD_EXECUTABLE(
  parallel_matching
  SOURCES parallel_matching.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
TRIBITS_ADD_EXECUTABLE(
  probe_tridiag
  SOURCES probe_tridiag.cpp
  DEPLIBS isorropia_test_utils
  COMM mpi
  )
TRIBITS_ADD_EXECUTABLE(
  order_example
  SOURCES order_example.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(IsorropiaExampleFiles
  SOURCE_FILES orderer.xml simple.mtx
  )
