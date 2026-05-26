#
# Executable
#

#TODO: how to remove that?
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

TRIBITS_ADD_EXECUTABLE(
  Redirection
  SOURCES Redirection.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  KokkosKernels
  SOURCES KokkosKernels.cpp
  COMM serial mpi
)

#
# Tests
#

TRIBITS_ADD_TEST(
  Redirection
  NAME "Redirection_Tpetra"
  ARGS "--linAlgebra=Tpetra"
  NUM_MPI_PROCS 4
  COMM mpi # HAVE_MPI required
)
