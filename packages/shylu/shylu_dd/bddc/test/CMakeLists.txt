
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

FILE(GLOB SOURCES *.cpp)

#removing due to dependence on a header not in Trilinos
#LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestDistributedIndexWithBulkData.cpp")

TRIBITS_ADD_EXECUTABLE(
    shylubddc_unit_tests
    SOURCES ${SOURCES}
    DEPLIBS
    COMM mpi
)

TRIBITS_ADD_TEST(
  shylubddc_unit_tests
  NAME bddc_standard_interface_test
  ARGS "--gtest_filter=SolverBDDC.Test1"
  COMM mpi
  NUM_MPI_PROCS 2
  FAIL_REGULAR_EXPRESSION "FAILED"
)		  

TRIBITS_ADD_TEST(
  shylubddc_unit_tests
  NAME bddc_simple_interface_test
  ARGS "--gtest_filter=SolverBDDCSimple.Test1"
  COMM mpi
  NUM_MPI_PROCS 2
  FAIL_REGULAR_EXPRESSION "FAILED"
)		  

TRIBITS_ADD_TEST(
  shylubddc_unit_tests
  NAME bddc_sparse_solver_test
  ARGS "--gtest_filter=SparseSolverBDDC.Test1"
  COMM mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "FAILED"
)		  

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUBDDCDriverFiles
    SOURCE_FILES sparseSolverTest.xml bddc.xml mueLu_SGS.xml elasticity3D.xml problemMaker.xml problemMakerSimple.xml
)

# Edit this comment and make trivial change if you add a file to this directory
