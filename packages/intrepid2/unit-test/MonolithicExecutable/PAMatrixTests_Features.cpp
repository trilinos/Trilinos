// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Tests of PAMatrix features: argument checking, constant-cell storage, matrix-valued point weights, assembly
            options, and GEMM statistics.
    \author Nathan V. Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_CellGeometryTestUtils.hpp>
#include <Intrepid2_HierarchicalBasisFamily.hpp>
#include <Intrepid2_SerendipityBasisFamily.hpp>
#include <Intrepid2_TestUtils.hpp>

#include "PAMatrixAssembly.hpp"

#include "StructuredIntegrationTests_TagDefs.hpp"
#include "StructuredIntegrationTests_Utils.hpp"

namespace
{
  using namespace Intrepid2;
  using DeviceType = DefaultTestDeviceType;
  using Scalar = double;
  using ExecutionSpace = typename DeviceType::execution_space;
  using HierarchicalFamily = HierarchicalBasisFamily<DeviceType>;

  template<int spaceDim>
  CellGeometry<Scalar,spaceDim,DeviceType> meshOfWidth(const int meshWidth, const AlgorithmChoice choice = Standard)
  {
    Kokkos::Array<int,spaceDim> gridCellCounts;
    for (int d=0; d<spaceDim; d++) gridCellCounts[d] = meshWidth;
    return getMesh<Scalar, spaceDim, DeviceType>(choice, gridCellCounts);
  }

  //! Copies an assembled matrix to a (C,F1,F2) device view (expanding constant-cell storage).
  ScalarView<Scalar,DeviceType> expandMatrix(const Data<Scalar,DeviceType> &matrix)
  {
    const int numCells = matrix.extent_int(0), F1 = matrix.extent_int(1), F2 = matrix.extent_int(2);
    ScalarView<Scalar,DeviceType> expanded("expanded matrix", numCells, F1, F2);
    Kokkos::parallel_for("expand matrix", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,F1,F2}),
    KOKKOS_LAMBDA(const int &c, const int &i, const int &j) { expanded(c,i,j) = matrix(c,i,j); });
    ExecutionSpace().fence();
    return expanded;
  }

  TEUCHOS_UNIT_TEST(PAMatrix, ArgumentChecks)
  {
    auto geometry = meshOfWidth<2>(2);
    Teuchos::RCP< Kokkos::Array<Scalar,2> > noWeight;
    auto paMatrix = constructPAMatrix<Scalar,HierarchicalFamily>(geometry, 2, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, noWeight,
                                                                 2, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, noWeight, false);
    const int numCells = geometry.numCells();
    const int F = paMatrix.leftBasisValues().numFields();

    // apply(): mismatched extents and insufficient workspace
    {
      auto workspace = paMatrix.allocateWorkspace(numCells);
      ScalarView<Scalar,DeviceType> input ("input",  numCells, F);
      ScalarView<Scalar,DeviceType> output("output", numCells, F);
      ScalarView<Scalar,DeviceType> wrongFields("wrong fields", numCells, F+1);
      ScalarView<Scalar,DeviceType> wrongCells ("wrong cells",  numCells+1, F);
      TEST_THROW(paMatrix.apply(output, wrongFields, workspace), std::invalid_argument);
      TEST_THROW(paMatrix.apply(wrongFields, input,  workspace), std::invalid_argument);
      TEST_THROW(paMatrix.apply(wrongCells, input,   workspace), std::invalid_argument);
      auto smallWorkspace = paMatrix.allocateWorkspace(1);
      TEST_THROW(paMatrix.apply(output, input, smallWorkspace), std::invalid_argument);
      TEST_NOTHROW(paMatrix.apply(output, input, smallWorkspace, false, 1));
    }

    // assemble() and extractDiagonal(): mismatched extents
    {
      Kokkos::View<Scalar***,DeviceType> wrong("wrong", numCells, F, F+1);
      Kokkos::Array<int,3> extents {numCells, F, F+1};
      Kokkos::Array<DataVariationType,3> variationTypes {GENERAL,GENERAL,GENERAL};
      Data<Scalar,DeviceType> wrongData(wrong, extents, variationTypes);
      TEST_THROW(paMatrix.assemble(wrongData), std::invalid_argument);

      Kokkos::View<Scalar**,DeviceType> wrongDiagonal("wrong diagonal", numCells, F+1);
      Kokkos::Array<int,2> diagonalExtents {numCells, F+1};
      Kokkos::Array<DataVariationType,2> diagonalVariationTypes {GENERAL,GENERAL};
      Data<Scalar,DeviceType> wrongDiagonalData(wrongDiagonal, diagonalExtents, diagonalVariationTypes);
      TEST_THROW(paMatrix.extractDiagonal(wrongDiagonalData), std::invalid_argument);
    }

    // orientations that are not aligned to the reference require the bases
    {
      ScalarView<Orientation,DeviceType> orientations("orientations", numCells);
      auto orientationsHost = Kokkos::create_mirror_view(orientations);
      ordinal_type edgeOrts[4] = {1,0,0,0};
      orientationsHost(0).setEdgeOrientation(4, edgeOrts);
      Kokkos::deep_copy(orientations, orientationsHost);
      using PAMatrixType = PAMatrix<DeviceType,Scalar>;
      auto construct = [&paMatrix](const ScalarView<Orientation,DeviceType> &orts) {
        PAMatrixType matrix(paMatrix.leftBasisValues(), paMatrix.cellMeasures(), paMatrix.rightBasisValues(), orts);
        return matrix.hasNontrivialOrientations();
      };
      TEST_THROW(construct(orientations), std::invalid_argument);
      // aligned orientations do not
      ScalarView<Orientation,DeviceType> alignedOrientations("aligned orientations", numCells);
      TEST_NOTHROW(construct(alignedOrientations));
    }
  }

  TEUCHOS_UNIT_TEST(PAMatrix, ConstantCellStorage)
  {
    // On a uniform grid, the matrix is the same for every cell, and allocateMatrixStorage() returns rank-2 storage.
    auto geometry = meshOfWidth<2>(3, Uniform);
    Teuchos::RCP< Kokkos::Array<Scalar,2> > noWeight;
    auto paMatrix = constructPAMatrix<Scalar,HierarchicalFamily>(geometry, 3, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, noWeight,
                                                                 3, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, noWeight, false);
    auto reference = paMatrix.allocateMatrixStorage();
    out << "underlying rank of matrix storage: " << reference.getUnderlyingViewRank() << "\n";
    TEST_EQUALITY(reference.getUnderlyingViewRank(), 2);
    referenceAssembly(paMatrix, reference);
    auto expected = expandMatrix(reference);

    using Algorithm = PAMatrixAssemblyOptions::Algorithm;
    for (const auto algorithm : {Algorithm::Automatic, Algorithm::UnitColumns, Algorithm::ApplyToColumns})
    {
      PAMatrixAssemblyOptions options;
      options.algorithm = algorithm;
      auto assembled = paMatrix.allocateMatrixStorage();
      paMatrix.assemble(assembled, options);
      out << "Comparing rank-2 assemble() (algorithm " << int(algorithm) << ")…\n";
      testFloatingEquality3(expected, expandMatrix(assembled), 1e-12, 1e-12, out, success, "reference", "assemble");
    }

    auto diagonal = paMatrix.allocateDiagonalStorage();
    TEST_EQUALITY(diagonal.getUnderlyingViewRank(), 1);
    paMatrix.extractDiagonal(diagonal);
    const int numCells = expected.extent_int(0), F = expected.extent_int(1);
    ScalarView<Scalar,DeviceType> expectedDiagonal("expected diagonal", numCells, F);
    Kokkos::parallel_for("expected diagonal", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,F}),
    KOKKOS_LAMBDA(const int &c, const int &i) { expectedDiagonal(c,i) = expected(c,i,i); });
    ExecutionSpace().fence();
    testFloatingEquality2(expectedDiagonal, diagonal, 1e-12, 1e-12, out, success, "expected diagonal", "extractDiagonal");
  }

  TEUCHOS_UNIT_TEST(PAMatrix, MatrixValuedWeights)
  {
    // (A grad u, grad v), with a symmetric, cell- and point-varying matrix A(C,P,D,D) with nonzero off-diagonal entries.
    const int spaceDim = 3;
    auto geometry = meshOfWidth<spaceDim>(2);
    const int numCells = geometry.numCells();
    const int p = 2;
    auto cellTopo = geometry.cellTopology();
    auto basis = getBasis<HierarchicalFamily>(cellTopo, FUNCTION_SPACE_HGRAD, p);
    auto cubature = DefaultCubatureFactory::create<DeviceType>(cellTopo, 2*p);
    auto points  = cubature->allocateCubaturePoints();
    auto weights = cubature->allocateCubatureWeights();
    cubature->getCubature(points, weights);
    auto gradients = basis->allocateBasisValues(points, OPERATOR_GRAD);
    basis->getValues(gradients, points, OPERATOR_GRAD);
    const int numPoints = points.extent_int(0);

    Kokkos::View<Scalar****,DeviceType> A("A", numCells, numPoints, spaceDim, spaceDim);
    Kokkos::parallel_for("fill A", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<4>>({0,0,0,0},{numCells,numPoints,spaceDim,spaceDim}),
    KOKKOS_LAMBDA(const int &c, const int &q, const int &i, const int &j) {
      A(c,q,i,j) = (i == j) ? 2.0 + 0.1 * c + 0.01 * q : 0.3 / (1.0 + i + j + 0.1 * c);
    });
    ExecutionSpace().fence();
    Kokkos::Array<int,4> extents {numCells, numPoints, spaceDim, spaceDim};
    Kokkos::Array<DataVariationType,4> variationTypes {GENERAL, GENERAL, GENERAL, GENERAL};
    Data<Scalar,DeviceType> AData(A, extents, variationTypes);

    // cell measures: the reference cubature weights, constant in the cell, as a (C,P) container
    auto weightsView = Kokkos::View<Scalar*,DeviceType>("weights", numPoints);
    {
      auto weightsHost = Kokkos::create_mirror_view(weightsView);
      for (int q=0; q<numPoints; q++)
      {
        // tensor-product weights: product over components (component 0 fastest)
        Scalar w = 1.0;
        int remainder = q;
        for (int r=0; r<weights.numTensorComponents(); r++)
        {
          auto component = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), weights.getTensorComponent(r).getUnderlyingView1());
          const int n = component.extent_int(0);
          w *= component(remainder % n);
          remainder /= n;
        }
        weightsHost(q) = w;
      }
      Kokkos::deep_copy(weightsView, weightsHost);
    }
    Kokkos::Array<int,2> measureExtents {numCells, numPoints};
    Kokkos::Array<DataVariationType,2> measureVariationTypes {CONSTANT, GENERAL};
    TensorData<Scalar,DeviceType> cellMeasures(Data<Scalar,DeviceType>(weightsView, measureExtents, measureVariationTypes));
    // the right transform is the identity, as a (C,P,D,D) matrix (PAMatrix requires matching ranks for vector-valued transforms)
    Kokkos::View<Scalar*,DeviceType> identityDiagonal("identity", spaceDim);
    Kokkos::deep_copy(identityDiagonal, 1.0);
    Kokkos::Array<DataVariationType,4> identityVariationTypes {CONSTANT, CONSTANT, BLOCK_PLUS_DIAGONAL, BLOCK_PLUS_DIAGONAL};
    Data<Scalar,DeviceType> identity(identityDiagonal, extents, identityVariationTypes);
    TransformedBasisValues<Scalar,DeviceType> left(AData, gradients);
    TransformedBasisValues<Scalar,DeviceType> right(identity, gradients);
    PAMatrix<DeviceType,Scalar> paMatrix(left, cellMeasures, right);

    auto reference = paMatrix.allocateMatrixStorage();
    referenceAssembly(paMatrix, reference);
    auto expected = expandMatrix(reference);

    using Algorithm = PAMatrixAssemblyOptions::Algorithm;
    for (const auto algorithm : {Algorithm::Automatic, Algorithm::ApplyToColumns})
    {
      PAMatrixAssemblyOptions options;
      options.algorithm = algorithm;
      auto assembled = paMatrix.allocateMatrixStorage();
      paMatrix.assemble(assembled, options);
      out << "Comparing assemble() with a matrix-valued weight (algorithm " << int(algorithm) << ")…\n";
      testFloatingEquality3(expected, expandMatrix(assembled), 1e-12, 1e-12, out, success, "reference", "assemble");
    }

    // apply() with worksets
    const int F = basis->getCardinality();
    ScalarView<Scalar,DeviceType> input ("input",  numCells, F, F);
    ScalarView<Scalar,DeviceType> output("output", numCells, F, F);
    Kokkos::parallel_for("set identity input", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,F}),
    KOKKOS_LAMBDA(const int &c, const int &j) { input(c,j,j) = 1.0; });
    ExecutionSpace().fence();
    const int worksetSize = 3;
    auto workspace = paMatrix.allocateWorkspace(worksetSize, F);
    paMatrix.apply(output, input, workspace, false, worksetSize);
    testFloatingEquality3(expected, output, 1e-12, 1e-12, out, success, "reference", "apply");
  }

  TEUCHOS_UNIT_TEST(PAMatrix, SparseLastStageOptions)
  {
    // The direct path's sparse and dense last stages, for a Serendipity basis.
    const int spaceDim = 3;
    auto geometry = meshOfWidth<spaceDim>(2);
    using SerendipityFamily = SerendipityBasisFamily<DeviceType>;
    Teuchos::RCP< Kokkos::Array<Scalar,spaceDim> > noWeight;
    auto paMatrix = constructPAMatrix<Scalar,SerendipityFamily>(geometry, 3, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, noWeight,
                                                                3, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, noWeight, false);
    using SparseLast = PAMatrixAssemblyOptions::SparseLastStage;
    PAMatrixAssemblyOptions options;
    options.sparseLastStage = SparseLast::Never;
    auto dense = paMatrix.allocateMatrixStorage();
    paMatrix.assemble(dense, options);
    for (const auto setting : {SparseLast::Always, SparseLast::Automatic})
    {
      options.sparseLastStage = setting;
      options.algorithm = PAMatrixAssemblyOptions::Algorithm::Automatic;
      auto assembled = paMatrix.allocateMatrixStorage();
      paMatrix.assemble(assembled, options);
      testFloatingEquality3(expandMatrix(dense), expandMatrix(assembled), 1e-12, 1e-12, out, success, "dense last stage", "sparse last stage");
    }
    options.algorithm = PAMatrixAssemblyOptions::Algorithm::UnitColumns;
    auto unitColumns = paMatrix.allocateMatrixStorage();
    paMatrix.assemble(unitColumns, options);
    testFloatingEquality3(expandMatrix(dense), expandMatrix(unitColumns), 1e-12, 1e-12, out, success, "direct", "unit columns");

    // Serendipity bases are not supported by apply() to columns
    options.algorithm = PAMatrixAssemblyOptions::Algorithm::ApplyToColumns;
    auto byApply = paMatrix.allocateMatrixStorage();
    TEST_THROW(paMatrix.assemble(byApply, options), std::invalid_argument);
    // nor by apply() or extractDiagonal()
    auto workspace = paMatrix.allocateWorkspace(geometry.numCells());
    const int F = paMatrix.leftBasisValues().numFields();
    ScalarView<Scalar,DeviceType> vector("vector", geometry.numCells(), F);
    TEST_THROW(paMatrix.apply(vector, vector, workspace), std::invalid_argument);
    auto diagonal = paMatrix.allocateDiagonalStorage();
    TEST_THROW(paMatrix.extractDiagonal(diagonal), std::invalid_argument);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, GEMMStatistics)
  {
    auto geometry = meshOfWidth<2>(2);
    Teuchos::RCP< Kokkos::Array<Scalar,2> > noWeight;
    auto paMatrix = constructPAMatrix<Scalar,HierarchicalFamily>(geometry, 2, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, noWeight,
                                                                 2, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, noWeight, false);
    using PAMatrixType = PAMatrix<DeviceType,Scalar>;
    const int numCells = geometry.numCells();
    const int F = paMatrix.leftBasisValues().numFields();
    auto workspace = paMatrix.allocateWorkspace(numCells);
    ScalarView<Scalar,DeviceType> input ("input",  numCells, F);
    ScalarView<Scalar,DeviceType> output("output", numCells, F);

    TEST_ASSERT(!PAMatrixType::gemmStatisticsEnabled());
    const double flopsBefore = PAMatrixType::gemmFlopCount();
    paMatrix.apply(output, input, workspace);
    TEST_EQUALITY(PAMatrixType::gemmFlopCount(), flopsBefore); // nothing recorded when disabled

    PAMatrixType::setGEMMStatisticsEnabled(true);
    const double baseFlops = PAMatrixType::gemmFlopCount();
    const double baseTime  = PAMatrixType::gemmTimeSeconds();
    paMatrix.apply(output, input, workspace);
    PAMatrixType::setGEMMStatisticsEnabled(false);
    const double flops = PAMatrixType::gemmFlopCount() - baseFlops;
    out << "recorded GEMM flops: " << flops << "; time: " << PAMatrixType::gemmTimeSeconds() - baseTime << " s\n";
    TEST_ASSERT(flops > 0);
    TEST_ASSERT(PAMatrixType::gemmTimeSeconds() >= baseTime);
    TEST_ASSERT(PAMatrixType::gemmThroughputGFlops(baseFlops, baseTime) >= 0);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, AssemblyOptionsFromEnvironment)
  {
    // Without the environment variables, the defaults are used.  (The variables are read once per process, so the
    // overrides themselves are tested by the ctest entries that set them.)
    const auto options = PAMatrixAssemblyOptions::defaultsWithEnvironmentOverrides();
    if (std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_ALGORITHM") == nullptr && std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_UNIT_COLUMNS") == nullptr)
    {
      TEST_ASSERT(options.algorithm == PAMatrixAssemblyOptions::Algorithm::Automatic);
    }
    if (std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_WORKSPACE_MB") == nullptr)
    {
      TEST_EQUALITY(options.workspaceMB, 512.0);
    }
  }
} // anonymous namespace
