
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB SOURCES *.cpp)

LIST(APPEND COMMON_SOURCES
  Akri_StkMeshBuilder.cpp
  Akri_StkMeshBuilder.hpp
  Akri_UnitMeshUtils.cpp
  Akri_UnitMeshUtils.hpp
  Akri_UnitTestUtils.cpp
  Akri_UnitTestUtils.hpp
  Akri_Unit_BoundingBoxMesh.hpp
  Akri_Unit_CreateFacetedSphere.cpp
  Akri_Unit_CreateFacetedSphere.hpp
  Akri_Unit_InterfaceGeometry.cpp
  Akri_Unit_InterfaceGeometry.hpp
  Akri_Unit_LogRedirecter.cpp
  Akri_Unit_LogRedirecter.hpp
  Akri_Unit_Part_Decomposition_Fixture.cpp
  Akri_Unit_Part_Decomposition_Fixture.hpp
  Akri_Unit_RefinementFixture.hpp
  Akri_Unit_RefinementFixture_Tet.hpp
  Akri_Unit_Single_Element_Fixtures.cpp
  Akri_Unit_Single_Element_Fixtures.hpp
  Akri_Unit_main.cpp
)

LIST(APPEND UNIT_TESTS
  Analytic_CDMesh
  AnalyticSurface
  BoundingBoxDistance
  CDFEM_Parent_Edge
  CDMesh
  Closest_Point_Intersection
  Constrained_Redistance
  ContourElement
  CurvatureLeastSquares
  DecomposeWithSensitivities
  Eikonal
  Element
  Element_Cutter
  Explicit_Hamilton_Jacobi
  FastMarching
  Geometry
  LowerEnvelope
  MortonIndex
  ParallelErrorMessage
  PatchInterpolator
  Refine_Beam
  Refine_CDMesh
  Refine_Edge_Tet
  Refine_General
  Refine_Hex
  RefineInterval
  Refine_Quad
  Refine_Tet
  Refine_Tri
  SearchTree
  SemiLagrangian
  SideAttachedElements
  SignedDistanceFromSTL
  SkinFacets
  Smoothing
  Snap
  TriangleCalcs
  WindingNumber
)

FOREACH(test IN LISTS UNIT_TESTS)

  LIST(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Akri_Unit_${test}.cpp)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    krino_unit_${test}
    SOURCES Akri_Unit_${test}.cpp ${COMMON_SOURCES}
    ARGS "--gtest_filter=-*fuzz*:OutputUtils*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX
    )

ENDFOREACH()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
    krino_unit_remainder
    SOURCES ${SOURCES}
    ARGS "--gtest_filter=-*fuzz*:OutputUtils*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX
    )
