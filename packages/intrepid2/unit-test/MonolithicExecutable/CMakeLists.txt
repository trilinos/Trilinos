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

# BasisEquivalenceTests
TRIBITS_ADD_EXECUTABLE( BasisEquivalenceTests NOEXESUFFIX SOURCES BasisEquivalenceTests_Wedge.cpp BasisEquivalenceTests_Tri.cpp BasisEquivalenceTests_Tet.cpp BasisEquivalenceTests_Quad.cpp BasisEquivalenceTests_Line.cpp BasisEquivalenceTests_Hypercube.cpp BasisEquivalenceTests_Hex.cpp UnitTestMain.cpp )
SET_TARGET_PROPERTIES( Intrepid2_BasisEquivalenceTests PROPERTIES EXCLUDE_FROM_ALL TRUE)

# BasisEquivalenceTests_Wedge
TRIBITS_ADD_EXECUTABLE( BasisEquivalenceTests_Wedge NOEXESUFFIX SOURCES BasisEquivalenceTests_Wedge.cpp UnitTestMain.cpp )
SET_TARGET_PROPERTIES( Intrepid2_BasisEquivalenceTests_Wedge PROPERTIES EXCLUDE_FROM_ALL TRUE)

# CellGeometry
TRIBITS_ADD_EXECUTABLE( CellGeometryTests NOEXESUFFIX SOURCES CellGeometryTests.cpp UnitTestMain.cpp )
SET_TARGET_PROPERTIES( Intrepid2_CellGeometryTests PROPERTIES EXCLUDE_FROM_ALL TRUE)

# DataTests
TRIBITS_ADD_EXECUTABLE( DataTests NOEXESUFFIX SOURCES DataTests.cpp UnitTestMain.cpp )
SET_TARGET_PROPERTIES( Intrepid2_DataTests PROPERTIES EXCLUDE_FROM_ALL TRUE)



