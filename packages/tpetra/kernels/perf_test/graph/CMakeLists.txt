INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_ADD_EXECUTABLE(
  coloring_test
  SOURCES KokkosKernels_GraphColoringTest.cpp KokkosKernels_GraphHelpers.hpp   
  )

TRIBITS_ADD_EXECUTABLE(
  PCG_TEST
  SOURCES KokkosKernels_PCG_GS_TEST.cpp KokkosKernels_PCG.hpp KokkosKernels_GraphHelpers.hpp
  )

TRIBITS_ADD_EXECUTABLE(
  conMTX2BIN
  SOURCES KokkosKernels_WriteBinaryFromMTX.cpp experiment_space.hpp KokkosKernels_GraphHelpers.hpp
  )
 
TRIBITS_ADD_EXECUTABLE(
  spemm
  SOURCES KokkosKernels_SPMM_KK.cpp  KokkosKernels_GraphHelpers.hpp
  )
  
  
#TRIBITS_ADD_EXECUTABLE(
#  spemm_cuda
#  SOURCES KokkosKernels_SPMM_cuda.cpp KokkosKernels_GraphHelpers.hpp 
#  )
  
#TRIBITS_ADD_EXECUTABLE(
#  symmetrize
#  SOURCES KokkosKernels_SymmetrizeTest.cpp KokkosKernels_GraphHelpers.hpp experiment_space.hpp
#  )  
  
#TRIBITS_ADD_EXECUTABLE(
#  view_test
#  SOURCES view_test.cpp
#)  

#TRIBITS_ADD_EXECUTABLE(
#  simple_test
#  SOURCES mytest.cpp
#)  
  
#TRIBITS_ADD_EXECUTABLE(
#  mp_test
#  SOURCES memory_pool_test.cpp
#)
 
#TRIBITS_ADD_EXECUTABLE(
#  reduce_test
#  SOURCES myReduceTest.cpp
#)
  
TRIBITS_ADD_EXECUTABLE(
  pool_test 
  SOURCES pooltest.cpp
)
 
#TRIBITS_ADD_EXECUTABLE(
#  myMemoryPoolTest
#  SOURCES myMemoryPoolTest.cpp
#)

#
# Make a trivial change here if you add a file either to this
# directory, or to the impl/ subdirectory.  Keep making changes
# whenever you add or remove a .cpp or .hpp file, and everything will
# be OK.
#

