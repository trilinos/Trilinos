#FILE(GLOB ETI_SOURCES "ETI/*cpp")


#FILE(GLOB TEST_SOURCES "*.cpp" ${ETI_SOURCES} )
FILE(GLOB TEST_SOURCES "*.cpp")

TRIBITS_INCLUDE_DIRECTORIES("../../assembly-examples")

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Intrepid2_Tests
  SOURCES ${TEST_SOURCES}
  ARGS 
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

# add single-group tests; allows for easier targeted builds and debugging (especially useful under CUDA)
# since each of these is redundant with Intrepid2_Tests above, we don't build them when all targets are built

# BasisEquivalenceTests
TRIBITS_ADD_EXECUTABLE( BasisEquivalenceTests NOEXESUFFIX SOURCES BasisEquivalenceTests_Wedge.cpp BasisEquivalenceTests_Tri.cpp BasisEquivalenceTests_Tet.cpp BasisEquivalenceTests_Quad.cpp BasisEquivalenceTests_Line.cpp BasisEquivalenceTests_Hypercube.cpp BasisEquivalenceTests_Hex.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_BasisEquivalenceTests)
  SET_TARGET_PROPERTIES( Intrepid2_BasisEquivalenceTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

# BasisEquivalenceTests_Pyramid
TRIBITS_ADD_EXECUTABLE( BasisEquivalenceTests_Pyramid NOEXESUFFIX SOURCES BasisEquivalenceTests_Pyramid.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_BasisEquivalenceTests_Pyramid)
  SET_TARGET_PROPERTIES( Intrepid2_BasisEquivalenceTests_Pyramid PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

# BasisValuesTests_Pyramid
TRIBITS_ADD_EXECUTABLE( BasisValuesTests_Pyramid NOEXESUFFIX SOURCES BasisValuesTests_Pyramid.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_BasisValuesTests_Pyramid)
  SET_TARGET_PROPERTIES( Intrepid2_BasisValuesTests_Pyramid PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

# BasisEquivalenceTests_Wedge
TRIBITS_ADD_EXECUTABLE( BasisEquivalenceTests_Wedge NOEXESUFFIX SOURCES BasisEquivalenceTests_Wedge.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_BasisEquivalenceTests_Wedge)
  SET_TARGET_PROPERTIES( Intrepid2_BasisEquivalenceTests_Wedge PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

# CellGeometry
TRIBITS_ADD_EXECUTABLE( CellGeometryTests NOEXESUFFIX SOURCES CellGeometryTests.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_CellGeometryTests)
  SET_TARGET_PROPERTIES( Intrepid2_CellGeometryTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

# DataTests
TRIBITS_ADD_EXECUTABLE( DataTests NOEXESUFFIX SOURCES DataTests.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_DataTests)
  SET_TARGET_PROPERTIES( Intrepid2_DataTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

# StructuredVersusStandard integration
TRIBITS_ADD_EXECUTABLE( StructuredVersusStandard NOEXESUFFIX SOURCES StructuredIntegrationTests_StructuredVersusStandard.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_StructuredVersusStandard)
  SET_TARGET_PROPERTIES( Intrepid2_StructuredVersusStandard PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

# SubBasisInclusionTests
TRIBITS_ADD_EXECUTABLE( SubBasisInclusionTests NOEXESUFFIX SOURCES SubBasisInclusionTests.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_SubBasisInclusionTests)
  SET_TARGET_PROPERTIES( Intrepid2_SubBasisInclusionTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()