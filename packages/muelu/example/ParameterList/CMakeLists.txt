INCLUDE(XpetraAddTest)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
#
#

#
# MLParameterList
#

TRIBITS_ADD_EXECUTABLE(
  MLParameterList
  SOURCES MLParameterList.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(MLParameterList_cp SOURCE_FILES ml_ParameterList.xml ml_RebalancedParameterList.xml)

#
# Tests of MLParameterList
#

#  XPETRA_ADD_TEST_TPETRA(MLParameterList 1-4)

# TODO: add a test with the --xml option
