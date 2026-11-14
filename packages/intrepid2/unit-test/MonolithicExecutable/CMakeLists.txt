#FILE(GLOB ETI_SOURCES "ETI/*cpp")


#FILE(GLOB TEST_SOURCES "*.cpp" ${ETI_SOURCES} )
FILE(GLOB TEST_SOURCES "*.cpp")

# PAMatrix requires KokkosKernels for GEMMs on non-host devices (PAMatrix on host devices uses Teuchos::BLAS)
SET(Intrepid2_PAMATRIX_TESTS_SUPPORTED ON)
IF (NOT ${PACKAGE_NAME}_ENABLE_KokkosKernels)
  IF (Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL OR Kokkos_ENABLE_OPENMPTARGET OR Kokkos_ENABLE_OPENACC)
    SET(Intrepid2_PAMATRIX_TESTS_SUPPORTED OFF)
  ENDIF()
ENDIF()
IF (NOT Intrepid2_PAMATRIX_TESTS_SUPPORTED)
  MESSAGE(STATUS "Intrepid2: PAMatrix tests disabled (a device execution space is enabled, and Intrepid2 is built without KokkosKernels)")
  LIST(FILTER TEST_SOURCES EXCLUDE REGEX "PAMatrix")
ENDIF()

TRIBITS_INCLUDE_DIRECTORIES("../../assembly-examples")

TRIBITS_ADD_EXECUTABLE (
  Monolithic
  NOEXESUFFIX
  SOURCES ${TEST_SOURCES}
  )

# Basis Equivalence Tests -- can be expensive; divided according to element types

SET( ElemTypes
  Line
  Quad
  Tri
  Hexahedron
  Tet
  Wedge
  Pyramid
  Hypercube
)

FOREACH( ElemType IN LISTS ElemTypes)
  TRIBITS_ADD_TEST(
    Monolithic
    NOEXESUFFIX
    NAME BasisEquivalence${ElemType}
    ARGS "--group=BasisEquivalence* --test=${ElemType}*"
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )
ENDFOREACH()

# Projected Geometry -- Linear tests
TRIBITS_ADD_TEST(
  Monolithic
  NOEXESUFFIX
  NAME ProjectedLinearGeometry
  ARGS "--group=ProjectedGeometry* --test=Linear*"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
)

# Projected Geometry -- Circular tests, single-cell
TRIBITS_ADD_TEST(
  Monolithic
  NOEXESUFFIX
  NAME ProjectedCircularGeometry_SingleCell
  ARGS "--group=ProjectedGeometry* --test=CircularGeometryConvergesInP_SingleCell"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
)

# Projected Geometry -- Circular tests, multi-cell
TRIBITS_ADD_TEST(
  Monolithic
  NOEXESUFFIX
  NAME ProjectedCircularGeometry_MultiCell
  ARGS "--group=ProjectedGeometry* --test=CircularGeometryConvergesInP_MultiCell"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  RUN_SERIAL
)

# other tests are divided such that they can be handled using --group argument alone
SET( TestGroups
  AnalyticPolynomialsMatch
  BasisCardinality
  BasisOpDecomposition
  BasisValues
  CellGeometry
  CubatureTensor
  Data
  DirectSumBasis
  HostCopyTests.cpp
  IntegratedJacobi
  IntegratedLegendre
  Jacobi
  Legendre
  RealSpaceTools
  StructuredIntegration
  SubBasisInclusion
  TensorBasis
  TensorTopologyMap
  TensorViewFunctor
  TransformedBasisValues
  VectorData
  ViewIterator
)

FOREACH( TestGroup IN LISTS TestGroups )
  TRIBITS_ADD_TEST(
    Monolithic
    NOEXESUFFIX
    NAME ${TestGroup}
    ARGS "--group=${TestGroup}*"
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )
ENDFOREACH()

IF (Intrepid2_PAMATRIX_TESTS_SUPPORTED)
  TRIBITS_ADD_TEST(
    Monolithic
    NOEXESUFFIX
    NAME PAMatrix
    ARGS "--group=PAMatrix*"
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )
  # assemble() algorithm selections, and a small workspace budget (forcing small worksets), via the environment
  TRIBITS_ADD_TEST(
    Monolithic
    NOEXESUFFIX
    NAME PAMatrix_AssembleUnitColumns
    ARGS "--group=PAMatrix"
    NUM_MPI_PROCS 1
    ENVIRONMENT "INTREPID2_PAMATRIX_ASSEMBLE_ALGORITHM=unit-columns"
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )
  TRIBITS_ADD_TEST(
    Monolithic
    NOEXESUFFIX
    NAME PAMatrix_AssembleSparseLastSmallWorkspace
    ARGS "--group=PAMatrix"
    NUM_MPI_PROCS 1
    ENVIRONMENT "INTREPID2_PAMATRIX_ASSEMBLE_SPARSE_LAST=1" "INTREPID2_PAMATRIX_ASSEMBLE_WORKSPACE_MB=0.01"
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )
ENDIF()

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

# StructuredIntegrationTests_QuadratureSynthetic.cpp
TRIBITS_ADD_EXECUTABLE( StructuredQuadratureSynthetic NOEXESUFFIX SOURCES StructuredIntegrationTests_QuadratureSynthetic.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_StructuredQuadratureSynthetic )
  SET_TARGET_PROPERTIES( Intrepid2_StructuredQuadratureSynthetic PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

IF (Intrepid2_PAMATRIX_TESTS_SUPPORTED)
# PAMatrixTests
TRIBITS_ADD_EXECUTABLE( PAMatrixTests NOEXESUFFIX SOURCES PAMatrixTests_Apply.cpp PAMatrixTests_ApplySynthetic.cpp PAMatrixTests_Features.cpp PAMatrix_Impl_Tests.cpp PAMatrixTests_Orientations.cpp PAMatrixTests_PAMatrixVersusStandardIntegration.cpp PAMatrixTests_Serendipity.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_PAMatrixTests)
  SET_TARGET_PROPERTIES( Intrepid2_PAMatrixTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

TRIBITS_ADD_EXECUTABLE( PAMatrixApplyTests NOEXESUFFIX SOURCES PAMatrixTests_Apply.cpp  UnitTestMain.cpp )
IF (TARGET Intrepid2_PAMatrixApplyTests)
  SET_TARGET_PROPERTIES( Intrepid2_PAMatrixApplyTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

TRIBITS_ADD_EXECUTABLE( PAMatrixApplySyntheticTests NOEXESUFFIX SOURCES PAMatrixTests_ApplySynthetic.cpp  UnitTestMain.cpp )
IF (TARGET Intrepid2_PAMatrixApplySyntheticTests)
  SET_TARGET_PROPERTIES( Intrepid2_PAMatrixApplySyntheticTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

TRIBITS_ADD_EXECUTABLE( PAMatrixOrientationTests NOEXESUFFIX SOURCES PAMatrixTests_Orientations.cpp  UnitTestMain.cpp )
IF (TARGET Intrepid2_PAMatrixOrientationTests)
  SET_TARGET_PROPERTIES( Intrepid2_PAMatrixOrientationTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()

TRIBITS_ADD_EXECUTABLE( PAMatrixImplTests NOEXESUFFIX SOURCES PAMatrix_Impl_Tests.cpp  UnitTestMain.cpp )
IF (TARGET Intrepid2_PAMatrixImplTests)
  SET_TARGET_PROPERTIES( Intrepid2_PAMatrixImplTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
ENDIF()
ENDIF()

# TensorBasisTagTests
TRIBITS_ADD_EXECUTABLE( TensorBasisTagTests NOEXESUFFIX SOURCES TensorBasisTagTests.cpp UnitTestMain.cpp )
IF (TARGET Intrepid2_TensorBasisTagTests)
  SET_TARGET_PROPERTIES( Intrepid2_TensorBasisTagTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
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
