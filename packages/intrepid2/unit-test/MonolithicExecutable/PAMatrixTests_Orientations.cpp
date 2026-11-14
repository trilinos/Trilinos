// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Tests of PAMatrix on meshes with mixed orientations (including rotated faces, whose orientation operators are
            not symmetric) and non-affine, cell-dependent geometry.
    \author Nathan V. Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_CellGeometryTestUtils.hpp>
#include <Intrepid2_HierarchicalBasisFamily.hpp>
#include <Intrepid2_OrientationTools.hpp>
#include <Intrepid2_TestUtils.hpp>

#include "PAMatrixAssembly.hpp"
#include "StandardAssembly.hpp"

#include "PAMatrixTests_Utils.hpp"
#include "StructuredIntegrationTests_TagDefs.hpp"
#include "StructuredIntegrationTests_Utils.hpp"

namespace
{
  using namespace Intrepid2;
  using DeviceType = DefaultTestDeviceType;
  using Scalar = double;
  using PointScalar = double;

  //! Orientations for numCells cells, cycling through: aligned to the reference; flipped edges and faces rotated by
  //! (nonsymmetric) quarter turns; and faces reflected.  Face orientations of hexahedra take values 0..7.
  ScalarView<Orientation,DeviceType> mixedOrientations(const int numCells, const int numEdges, const int numFaces)
  {
    ScalarView<Orientation,DeviceType> orientations("mixed orientations", numCells);
    auto orientationsHost = Kokkos::create_mirror_view(orientations);
    for (int cell=0; cell<numCells; cell++)
    {
      ordinal_type edgeOrts[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
      ordinal_type faceOrts[6]  = {0,0,0,0,0,0};
      const int variant = cell % 3;
      if (variant != 0)
      {
        for (int e=0; e<numEdges; e++) edgeOrts[e] = (e + cell) % 2;
        for (int f=0; f<numFaces; f++) faceOrts[f] = (variant == 1) ? (1 + 2 * ((f + cell) % 2)) // 1 or 3: rotations
                                                                    : (4 + ((f + cell) % 4));     // 4..7: reflections
      }
      Orientation orientation;
      if (numEdges > 0) orientation.setEdgeOrientation(numEdges, edgeOrts);
      if (numFaces > 0) orientation.setFaceOrientation(numFaces, faceOrts);
      orientationsHost(cell) = orientation;
    }
    Kokkos::deep_copy(orientations, orientationsHost);
    return orientations;
  }

  //! The orientation operator O, with shape (C,F,F): OrientationTools::modifyBasisByOrientation() applied to the identity.
  //! The oriented basis is phi~_i = sum_j O(i,j) phi_j.
  template<class BasisType>
  ScalarView<Scalar,DeviceType> orientationOperator(const ScalarView<Orientation,DeviceType> &orientations, const BasisType *basis)
  {
    const int numCells  = orientations.extent_int(0);
    const int numFields = basis->getCardinality();
    ScalarView<Scalar,DeviceType> identity("identity", numFields, numFields);
    auto identityHost = Kokkos::create_mirror_view(identity);
    for (int i=0; i<numFields; i++) identityHost(i,i) = 1.0;
    Kokkos::deep_copy(identity, identityHost);
    ScalarView<Scalar,DeviceType> O("orientation operator", numCells, numFields, numFields);
    OrientationTools<DeviceType>::modifyBasisByOrientation(O, identity, orientations, basis);
    return O;
  }

  //! Copies an assembled matrix to a (C,F1,F2) device view (expanding constant-cell storage).
  ScalarView<Scalar,DeviceType> expandMatrix(const Data<Scalar,DeviceType> &matrix)
  {
    const int numCells = matrix.extent_int(0), F1 = matrix.extent_int(1), F2 = matrix.extent_int(2);
    ScalarView<Scalar,DeviceType> expanded("expanded matrix", numCells, F1, F2);
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for("expand matrix", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,F1,F2}),
    KOKKOS_LAMBDA(const int &c, const int &i, const int &j) { expanded(c,i,j) = matrix(c,i,j); });
    ExecutionSpace().fence();
    return expanded;
  }

  //! Checks assemble() (with each algorithm), apply() (with and without worksets), and extractDiagonal() against expected.
  void checkPAMatrixOperations(const PAMatrix<DeviceType,Scalar> &paMatrix, const ScalarView<Scalar,DeviceType> &expected,
                               const double relTol, const double absTol, Teuchos::FancyOStream &out, bool &success)
  {
    using ExecutionSpace = typename DeviceType::execution_space;
    const int numCells = expected.extent_int(0), F1 = expected.extent_int(1), F2 = expected.extent_int(2);

    using Algorithm = PAMatrixAssemblyOptions::Algorithm;
    for (const auto algorithm : {Algorithm::Automatic, Algorithm::UnitColumns, Algorithm::ApplyToColumns})
    {
      for (const double workspaceMB : {512.0, 1e-3}) // the small budget forces one cell per workset
      {
        PAMatrixAssemblyOptions options;
        options.algorithm   = algorithm;
        options.workspaceMB = workspaceMB;
        auto assembled = paMatrix.allocateMatrixStorage();
        paMatrix.assemble(assembled, options);
        out << "Comparing assemble() (algorithm " << int(algorithm) << ", workspace " << workspaceMB << " MB)…\n";
        testFloatingEquality3(expected, expandMatrix(assembled), relTol, absTol, out, success, "expected", "assemble");
      }
    }

    // apply(), to all unit columns at once, with and without a workset size that does not divide the cell count
    for (const int worksetSize : {0, std::max(1, (numCells + 1) / 2 - 1), 1})
    {
      ScalarView<Scalar,DeviceType> input ("input",  numCells, F2, F2);
      ScalarView<Scalar,DeviceType> output("output", numCells, F1, F2);
      Kokkos::parallel_for("set identity input", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,F2}),
      KOKKOS_LAMBDA(const int &c, const int &j) { input(c,j,j) = 1.0; });
      ExecutionSpace().fence();
      auto workspace = paMatrix.allocateWorkspace((worksetSize > 0) ? worksetSize : numCells, F2);
      paMatrix.apply(output, input, workspace, false, worksetSize);
      out << "Comparing apply() (workset size " << worksetSize << ")…\n";
      testFloatingEquality3(expected, output, relTol, absTol, out, success, "expected", "apply");
    }

    if (F1 == F2)
    {
      auto diagonal = paMatrix.allocateDiagonalStorage();
      paMatrix.extractDiagonal(diagonal);
      ScalarView<Scalar,DeviceType> expectedDiagonal("expected diagonal", numCells, F1);
      Kokkos::parallel_for("expected diagonal", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,F1}),
      KOKKOS_LAMBDA(const int &c, const int &i) { expectedDiagonal(c,i) = expected(c,i,i); });
      ExecutionSpace().fence();
      out << "Comparing extractDiagonal()…\n";
      testFloatingEquality2(expectedDiagonal, diagonal, relTol, absTol, out, success, "expected diagonal", "extractDiagonal");
    }
  }

  /** \brief Compares PAMatrix with orientations to standard assembly (which orients the basis values before integration) on a
             perturbed mesh whose orientations are determined by permuted global vertex numbers; and, with orientations set
             explicitly to cycle through rotations and reflections, to O_L M O_R^T, computed from explicit orientation
             operators and the unoriented matrix M.
  */
  template<class BasisFamily, int spaceDim>
  void testOrientedPAMatrix(const int meshWidth,
                            const EFunctionSpace &fs1, const EOperator &op1, const int &p1,
                            const EFunctionSpace &fs2, const EOperator &op2, const int &p2,
                            const double relTol, const double absTol,
                            Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::Array<int,spaceDim> gridCellCounts;
    for (int d=0; d<spaceDim; d++) gridCellCounts[d] = meshWidth;
    auto geometry = perturbedMixedOrientationMesh<PointScalar,spaceDim,DeviceType>(gridCellCounts);
    const int numCells = geometry.numCells();
    Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > noWeight;

    // (a) versus standard assembly, with the mesh's orientations
    {
      auto paMatrix = constructPAMatrix<Scalar,BasisFamily>(geometry, p1, fs1, op1, noWeight, p2, fs2, op2, noWeight);
      double flopsIntegration = 0, flopsJacobian = 0;
      const int worksetSize = numCells;
      auto standard = performStandardAssembly<Scalar,BasisFamily>(geometry, worksetSize, p1, fs1, op1, noWeight, p2, fs2, op2, noWeight,
                                                                  flopsIntegration, flopsJacobian);
      out << "Mesh orientations (" << (paMatrix.hasNontrivialOrientations() ? "some" : "none") << " not aligned to the reference); comparing to standard assembly.\n";
      TEST_ASSERT(paMatrix.hasNontrivialOrientations() || (spaceDim == 1));
      checkPAMatrixOperations(paMatrix, standard, relTol, absTol, out, success);
    }

    // (b) versus O_L M O_R^T, with explicitly set orientations
    {
      const bool applyOrientations = false;
      auto unorientedMatrix = constructPAMatrix<Scalar,BasisFamily>(geometry, p1, fs1, op1, noWeight, p2, fs2, op2, noWeight,
                                                                    applyOrientations);
      auto cellTopo = geometry.cellTopology();
      auto orientations = mixedOrientations(numCells, (spaceDim > 1) ? cellTopo.getEdgeCount() : 0, (spaceDim == 3) ? cellTopo.getFaceCount() : 0);
      auto leftBasis  = getBasis<BasisFamily>(cellTopo, fs1, p1);
      auto rightBasis = getBasis<BasisFamily>(cellTopo, fs2, p2);
      PAMatrix<DeviceType,Scalar> orientedMatrix(unorientedMatrix.leftBasisValues(), unorientedMatrix.cellMeasures(),
                                                 unorientedMatrix.rightBasisValues(), orientations, leftBasis, rightBasis);

      auto M = unorientedMatrix.allocateMatrixStorage();
      referenceAssembly(unorientedMatrix, M);
      auto OL = orientationOperator(orientations, leftBasis.get());
      auto OR = orientationOperator(orientations, rightBasis.get());
      auto MHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), expandMatrix(M));
      auto OLHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), OL);
      auto ORHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), OR);
      const int F1 = leftBasis->getCardinality();
      const int F2 = rightBasis->getCardinality();
      ScalarView<Scalar,DeviceType> expected("expected", numCells, F1, F2);
      auto expectedHost = Kokkos::create_mirror_view(expected);
      std::vector<Scalar> MOT(size_t(F1) * F2);
      for (int c=0; c<numCells; c++)
      {
        for (int i=0; i<F1; i++)
          for (int j=0; j<F2; j++)
          {
            Scalar value = 0;
            for (int k=0; k<F2; k++) value += MHost(c,i,k) * ORHost(c,j,k);
            MOT[i + j * F1] = value;
          }
        for (int i=0; i<F1; i++)
          for (int j=0; j<F2; j++)
          {
            Scalar value = 0;
            for (int k=0; k<F1; k++) value += OLHost(c,i,k) * MOT[k + j * F1];
            expectedHost(c,i,j) = value;
          }
      }
      Kokkos::deep_copy(expected, expectedHost);

      auto referenceOriented = orientedMatrix.allocateMatrixStorage();
      referenceAssembly(orientedMatrix, referenceOriented);
      out << "Comparing referenceAssembly() (IntegrationTools + OrientationTools::modifyMatrixByOrientation()) to O_L M O_R^T…\n";
      testFloatingEquality3(expected, expandMatrix(referenceOriented), relTol, absTol, out, success, "O_L M O_R^T", "referenceAssembly");

      out << "Comparing PAMatrix operations with explicit orientations to O_L M O_R^T…\n";
      checkPAMatrixOperations(orientedMatrix, expected, relTol, absTol, out, success);
    }
  }

  using HierarchicalFamily = HierarchicalBasisFamily<DeviceType>;
  using NodalFamily        = DerivedNodalBasisFamily<DeviceType>;

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Hex_HGRAD_p3)
  {
    testOrientedPAMatrix<HierarchicalFamily,3>(2, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, 3, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, 3,
                                               1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Hex_HGRAD_p3_p4_values)
  {
    testOrientedPAMatrix<HierarchicalFamily,3>(2, FUNCTION_SPACE_HGRAD, OPERATOR_VALUE, 3, FUNCTION_SPACE_HGRAD, OPERATOR_VALUE, 4,
                                               1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Hex_HGRAD_Nodal_p3)
  {
    testOrientedPAMatrix<NodalFamily,3>(2, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, 3, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, 3,
                                        1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Hex_HCURL_p3)
  {
    testOrientedPAMatrix<HierarchicalFamily,3>(2, FUNCTION_SPACE_HCURL, OPERATOR_CURL, 3, FUNCTION_SPACE_HCURL, OPERATOR_CURL, 3,
                                               1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Hex_HCURL_values_p3)
  {
    testOrientedPAMatrix<HierarchicalFamily,3>(2, FUNCTION_SPACE_HCURL, OPERATOR_VALUE, 3, FUNCTION_SPACE_HCURL, OPERATOR_VALUE, 3,
                                               1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Hex_HDIV_p3)
  {
    testOrientedPAMatrix<HierarchicalFamily,3>(2, FUNCTION_SPACE_HDIV, OPERATOR_VALUE, 3, FUNCTION_SPACE_HDIV, OPERATOR_VALUE, 3,
                                               1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Hex_HDIV_DIV_p3)
  {
    testOrientedPAMatrix<HierarchicalFamily,3>(2, FUNCTION_SPACE_HDIV, OPERATOR_DIV, 3, FUNCTION_SPACE_HDIV, OPERATOR_DIV, 3,
                                               1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Hex_HDIV_HGRAD_p3)
  {
    testOrientedPAMatrix<HierarchicalFamily,3>(2, FUNCTION_SPACE_HDIV, OPERATOR_DIV, 3, FUNCTION_SPACE_HGRAD, OPERATOR_VALUE, 2,
                                               1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Hex_HCURL_HDIV_p2)
  {
    testOrientedPAMatrix<HierarchicalFamily,3>(2, FUNCTION_SPACE_HCURL, OPERATOR_CURL, 2, FUNCTION_SPACE_HDIV, OPERATOR_VALUE, 2,
                                               1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Quad_HCURL_p3)
  {
    testOrientedPAMatrix<HierarchicalFamily,2>(3, FUNCTION_SPACE_HCURL, OPERATOR_VALUE, 3, FUNCTION_SPACE_HCURL, OPERATOR_VALUE, 3,
                                               1e-11, 1e-11, out, success);
  }

  TEUCHOS_UNIT_TEST(PAMatrix, Orientations_Quad_HDIV_HGRAD_p3)
  {
    testOrientedPAMatrix<HierarchicalFamily,2>(3, FUNCTION_SPACE_HDIV, OPERATOR_VALUE, 3, FUNCTION_SPACE_HGRAD, OPERATOR_GRAD, 3,
                                               1e-11, 1e-11, out, success);
  }
} // anonymous namespace
