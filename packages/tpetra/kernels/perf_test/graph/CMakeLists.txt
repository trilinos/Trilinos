INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_ADD_EXECUTABLE(
  coloring_test
  SOURCES coloring_test.cpp experiment_space.hpp KokkosKernelsGraphHelpers.hpp   
  )

TRIBITS_ADD_EXECUTABLE(
  gs_test
  SOURCES gs_test.cpp experiment_space.hpp pcg.hpp KokkosKernelsGraphHelpers.hpp
  )

TRIBITS_ADD_EXECUTABLE(
  write_binary_graph
  SOURCES write_binary_graph.cpp experiment_space.hpp KokkosKernelsGraphHelpers.hpp
  )

TRIBITS_ADD_EXECUTABLE(
  spemm_cusparse
  SOURCES SPMM.cpp KokkosKernelsGraphHelpers.hpp 
  )
  
TRIBITS_ADD_EXECUTABLE(
  symmetrize
  SOURCES symmetrize_test.cpp KokkosKernelsGraphHelpers.hpp experiment_space.hpp
  )  
  
#
# Make a trivial change here if you add a file either to this
# directory, or to the impl/ subdirectory.  Keep making changes
# whenever you add or remove a .cpp or .hpp file, and everything will
# be OK.
#

