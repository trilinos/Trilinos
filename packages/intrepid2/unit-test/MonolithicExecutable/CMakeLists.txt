#FILE(GLOB ETI_SOURCES "ETI/*cpp")


#FILE(GLOB TEST_SOURCES "*.cpp" ${ETI_SOURCES} )
FILE(GLOB TEST_SOURCES "*.cpp")

INCLUDE_DIRECTORIES("../../assembly-examples")

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Intrepid2_Tests
  SOURCES ${TEST_SOURCES}
  ARGS 
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

# add single-group tests; allows for easier targeted builds and debugging (especially useful under CUDA)
# since each of these is redundant with Intrepid2_Tests above, we don't build this when all targets are built

# CellGeometry
TRIBITS_ADD_EXECUTABLE( CellGeometryTests NOEXESUFFIX SOURCES CellGeometryTests.cpp UnitTestMain.cpp )
SET_TARGET_PROPERTIES( Intrepid2_CellGeometryTests PROPERTIES EXCLUDE_FROM_ALL TRUE)

# DataTests
TRIBITS_ADD_EXECUTABLE( DataTests NOEXESUFFIX SOURCES DataTests.cpp UnitTestMain.cpp )
SET_TARGET_PROPERTIES( Intrepid2_DataTests PROPERTIES EXCLUDE_FROM_ALL TRUE)