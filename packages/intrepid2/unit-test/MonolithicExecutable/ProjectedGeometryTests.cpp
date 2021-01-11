// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   ProjectedGeometryTests.cpp
    \brief  Tests against the ProjectedGeometry and ProjectedGeometryExamples classes.
    \author Nathan V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_ProjectedGeometry.hpp"
#include "Intrepid2_ProjectedGeometryExamples.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <Kokkos_Core.hpp>

namespace
{
  using namespace Intrepid2;

  template<typename PointScalar, typename ExecutionSpace>
  void testCircularGeometryProjectionConvergesInP(int meshWidth, Teuchos::FancyOStream &out, bool &success)
  {
    constexpr int spaceDim = 2;
    using BasisFamily = DerivedNodalBasisFamily<ExecutionSpace,PointScalar,PointScalar>;
    using ProjectedGeometry = ProjectedGeometry<spaceDim,PointScalar,ExecutionSpace>;
    
    constexpr PointScalar R = 1.0; // unit circle centered at origin
    constexpr PointScalar exactArea = M_PI * R * R;
    
    UnitSquareToCircle<PointScalar> circularGeometry;
    
    Kokkos::Array<PointScalar,spaceDim> origin;
    Kokkos::Array<PointScalar,spaceDim> extent;
    Kokkos::Array<int,spaceDim> cellCount;
    
    int numCells = 1;
    for (int d=0; d<spaceDim; d++)
    {
      origin[d] = -1.0;
      extent[d] =  2.0;
      cellCount[d] = meshWidth;
      numCells *= meshWidth;
    }
    
    CellGeometry<PointScalar, spaceDim> flatCellGeometry(origin,extent,cellCount);
    shards::CellTopology cellTopo = flatCellGeometry.cellTopology();
    
    auto cubature = Intrepid2::DefaultCubatureFactory::create<ExecutionSpace>(cellTopo,Intrepid2::Parameters::MaxOrder);
    auto cubatureWeights = cubature->allocateCubatureWeights();
    TensorPoints<PointScalar,ExecutionSpace> cubaturePoints  = cubature->allocateCubaturePoints();
    cubature->getCubature(cubaturePoints, cubatureWeights);
    
    const int numPoints = cubaturePoints.extent_int(0);
    
    PointScalar previousError = exactArea; // we check that for each p, its error is less than that for (p-1); initialize as 100% error
    for (int p=1; p<=Intrepid2::Parameters::MaxOrder; p++)
    {
      auto basisForNodes = getBasis<BasisFamily>(cellTopo, FUNCTION_SPACE_HGRAD, p);
      
      const int nodesPerCell = basisForNodes->getCardinality();
      
      ViewType<PointScalar> projectedNodes("projected nodes", numCells, nodesPerCell, spaceDim);
      
      ProjectedGeometry::projectOntoHGRADBasis(projectedNodes, basisForNodes, flatCellGeometry, circularGeometry, circularGeometry);

      CellGeometry<PointScalar, spaceDim> projectedGeometry(basisForNodes,projectedNodes);

      auto jacobian = projectedGeometry.allocateJacobianData(cubaturePoints);
      auto jacobianDet = CellTools<ExecutionSpace>::allocateJacobianDet(jacobian);
      
      // sanity checks on the allocations:
      TEST_EQUALITY(numCells,  jacobian.extent_int(0));
      TEST_EQUALITY(numPoints, jacobian.extent_int(1));
      TEST_EQUALITY(spaceDim,  jacobian.extent_int(2));
      TEST_EQUALITY(spaceDim,  jacobian.extent_int(3));
      TEST_EQUALITY(numCells,  jacobianDet.extent_int(0));
      TEST_EQUALITY(numPoints, jacobianDet.extent_int(1));
      
      projectedGeometry.setJacobian(jacobian, cubaturePoints);
      CellTools<ExecutionSpace>::setJacobianDet(jacobianDet, jacobian);
      
      auto cellMeasure = projectedGeometry.allocateCellMeasure(jacobianDet, cubatureWeights);
      projectedGeometry.computeCellMeasure(cellMeasure, jacobianDet, cubatureWeights);
          
      PointScalar actualArea = 0;
      Kokkos::RangePolicy<ExecutionSpace > reducePolicy(0, numCells);
      Kokkos::parallel_reduce( reducePolicy,
      KOKKOS_LAMBDA( const int &cellOrdinal, PointScalar &reducedValue )
      {
        PointScalar localSum = 0;
        for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
        {
          localSum += cellMeasure(cellOrdinal,pointOrdinal);
        }
        reducedValue += localSum;
      }, actualArea);
      
      PointScalar error = std::abs(actualArea - exactArea);
      std::cout << "error for p = " << p << ": " << error << std::endl;
      
      TEST_COMPARE(error, <, previousError);
    }
  }

  template<typename PointScalar, int spaceDim, typename ExecutionSpace>
  void testLinearGeometryIsExact(int meshWidth, const double &relTol, const double &absTol,
                                 Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::Array<PointScalar,spaceDim> origin;
    Kokkos::Array<PointScalar,spaceDim> extent;
    Kokkos::Array<int,spaceDim> cellCount;
    
    for (int d=0; d<spaceDim; d++)
    {
      origin[d] = 0.0;
      extent[d] = 1.0;
      cellCount[d] = meshWidth;
    }
    
    CellGeometry<PointScalar, spaceDim> flatCellGeometry(origin,extent,cellCount);
    auto linearBasis = flatCellGeometry.basisForNodes();
    
    ProjectedGeometryIdentityMap<PointScalar, spaceDim> exactGeometry;
    
    int numCells = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numCells *= meshWidth;
    }
    const int nodesPerCell = linearBasis->getCardinality();
    
    ViewType<PointScalar> projectedNodes("projected nodes", numCells, nodesPerCell, spaceDim);
    
    using ProjectedGeometry = ProjectedGeometry<spaceDim,PointScalar,ExecutionSpace>;
    ProjectedGeometry::projectOntoHGRADBasis(projectedNodes, linearBasis, flatCellGeometry, exactGeometry, exactGeometry);
    
    printFunctor3(projectedNodes, out, "projectedNodes");
    testFloatingEquality3(projectedNodes, flatCellGeometry, relTol, absTol, out, success, "projected geometry", "original CellGeometry");
  }

  template<typename PointScalar, int spaceDim, typename ExecutionSpace>
  void testLinearGeometryIsExact(int meshWidth, int polyOrderForBasis, const double &relTol, const double &absTol,
                                 Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::Array<PointScalar,spaceDim> origin;
    Kokkos::Array<PointScalar,spaceDim> extent;
    Kokkos::Array<int,spaceDim> cellCount;
    
    for (int d=0; d<spaceDim; d++)
    {
      origin[d] = 0.0;
      extent[d] = 1.0;
      cellCount[d] = meshWidth;
    }
    
    CellGeometry<PointScalar, spaceDim> flatCellGeometry(origin,extent,cellCount);
    
    using BasisPtr = Teuchos::RCP< Basis<ExecutionSpace,PointScalar,PointScalar> >;
    
    using BasisFamily = DerivedNodalBasisFamily<ExecutionSpace,PointScalar,PointScalar>;
    BasisPtr hgradBasisForProjection = getBasis<BasisFamily>(flatCellGeometry.cellTopology(), FUNCTION_SPACE_HGRAD, polyOrderForBasis);
    
    ProjectedGeometryIdentityMap<PointScalar, spaceDim> exactGeometry;
    
    int numCells = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numCells *= meshWidth;
    }
    const int nodesPerCell = hgradBasisForProjection->getCardinality();
    ViewType<PointScalar> projectedNodes("projected nodes", numCells, nodesPerCell, spaceDim);
    
    using ProjectedGeometry = ProjectedGeometry<spaceDim,PointScalar,ExecutionSpace>;
    ProjectedGeometry::projectOntoHGRADBasis(projectedNodes, hgradBasisForProjection, flatCellGeometry, exactGeometry, exactGeometry);
    
    auto cellTopo = flatCellGeometry.cellTopology();
    
    // get reference-element points at which we can evaluate the projected basis
    const int numNodesPerFlatCell = flatCellGeometry.numNodesPerCell();
    ScalarView<PointScalar,ExecutionSpace> refCellNodes("ref cell nodes", numNodesPerFlatCell, spaceDim);
    
    for (int node=0; node<numNodesPerFlatCell; node++)
    {
      auto nodeSubview = Kokkos::subdynrankview(refCellNodes, node, Kokkos::ALL());
      CellTools<ExecutionSpace>::getReferenceNode(nodeSubview, cellTopo, node);
    }
    
    auto hgradValuesAtNodes = hgradBasisForProjection->allocateOutputView(numNodesPerFlatCell);
    hgradBasisForProjection->getValues(hgradValuesAtNodes, refCellNodes, OPERATOR_VALUE);

    ViewType<PointScalar> evaluatedNodes("projected nodes evaluated on flat cell", numCells, numNodesPerFlatCell, spaceDim);
    
    auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,numNodesPerFlatCell});
    Kokkos::parallel_for("evaluate projected nodes", policy,
    KOKKOS_LAMBDA (const int &cellOrdinal, const int &nodeOrdinalInFlatCell) {
      for (int d=0; d<spaceDim; d++)
      {
        PointScalar nodeCoord = 0.0;
        for (int highOrderNode=0; highOrderNode<nodesPerCell; highOrderNode++)
        {
          const auto & basisWeight = projectedNodes(cellOrdinal,highOrderNode,d);
          const auto & basisValue  = hgradValuesAtNodes(highOrderNode,nodeOrdinalInFlatCell);
          nodeCoord += basisWeight * basisValue;
        }
        evaluatedNodes(cellOrdinal,nodeOrdinalInFlatCell,d) = nodeCoord;
      }
    });
    ExecutionSpace().fence();
    
    printFunctor3(projectedNodes, out, "projectedNodes");
    printFunctor3(evaluatedNodes, out, "evaluatedNodes");
    
    // need to use "classic" node ordering for the comparison below because the refCellNodes above use the classic node ordering
    auto subdivisionStrategy = CellGeometry<PointScalar, spaceDim>::NO_SUBDIVISION;
    auto nodeOrdering = CellGeometry<PointScalar, spaceDim>::HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS;
    
    CellGeometry<PointScalar, spaceDim> flatCellGeometryClassic(origin,extent,cellCount, subdivisionStrategy, nodeOrdering);
    
    testFloatingEquality3(evaluatedNodes, flatCellGeometryClassic, relTol, absTol, out, success, "projected geometry", "original CellGeometry");
  }

  TEUCHOS_UNIT_TEST( ProjectedGeometry, CircularGeometryConvergesInP_SingleCell )
  {
    using Scalar = double;
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    const int meshWidth = 1;
    testCircularGeometryProjectionConvergesInP<Scalar, ExecSpace>(meshWidth, out, success);
  }

  TEUCHOS_UNIT_TEST( ProjectedGeometry, CircularGeometryConvergesInP_MultiCell )
  {
    using Scalar = double;
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    const int meshWidth = 4;
    testCircularGeometryProjectionConvergesInP<Scalar, ExecSpace>(meshWidth, out, success);
  }

  TEUCHOS_UNIT_TEST( ProjectedGeometry, LinearGeometryIsExact_Line )
  {
    // 1D test that linear geometry is exactly recovered by the projection
    using PointScalar = double;
    const int spaceDim = 1;
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    
    const int meshWidth = 4;
    
    // these tolerances are fairly tight; we may need to loosen them to pass on all platforms
    const double relTol = 1e-15;
    const double absTol = 1e-15;
    testLinearGeometryIsExact<PointScalar, spaceDim, ExecSpace>(meshWidth, relTol, absTol, out, success);
    
    // now, test with a higher-order basis (but still linear geometry)
    for (int polyOrder=2; polyOrder<8; polyOrder++)
    {
      testLinearGeometryIsExact<PointScalar, spaceDim, ExecSpace>(meshWidth, polyOrder, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( ProjectedGeometry, LinearGeometryIsExact_Quad )
  {
    // 2D test that linear geometry is exactly recovered by the projection
    using PointScalar = double;
    const int spaceDim = 2;
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    
    const int meshWidth = 4;
    
    // these tolerances are fairly tight; we may need to loosen them to pass on all platforms
    const double relTol = 1e-14;
    const double absTol = 1e-14;
    testLinearGeometryIsExact<PointScalar, spaceDim, ExecSpace>(meshWidth, relTol, absTol, out, success);
    
    // now, test with a higher-order basis (but still linear geometry)
    for (int polyOrder=2; polyOrder<5; polyOrder++)
    {
      testLinearGeometryIsExact<PointScalar, spaceDim, ExecSpace>(meshWidth, polyOrder, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( ProjectedGeometry, LinearGeometryIsExact_Hex )
  {
    // 3D test that linear geometry is exactly recovered by the projection
    using PointScalar = double;
    const int spaceDim = 3;
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    
    const int meshWidth = 2;
    
    // these tolerances are fairly tight; we may need to loosen them to pass on all platforms
    const double relTol = 1e-14;
    const double absTol = 1e-14;
    testLinearGeometryIsExact<PointScalar, spaceDim, ExecSpace>(meshWidth, relTol, absTol, out, success);
    
    // now, test with a higher-order basis (but still linear geometry)
    for (int polyOrder=2; polyOrder<4; polyOrder++)
    {
      testLinearGeometryIsExact<PointScalar, spaceDim, ExecSpace>(meshWidth, polyOrder, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( ProjectedGeometryExamples, CircularGeometry )
  {
    using Scalar = double;
    UnitSquareToCircle<Scalar> circularMap;
    const int spaceDim = 2;
    
    const double tol = 1e-14;
    
    // we expect (-1,0), (0,-1), (1,0), and (0,1) to map to themselves
    using Point = Kokkos::Array<Scalar,2>;
    
    Point left   { -1,  0};
    Point bottom {  0, -1};
    Point right  {  1,  0};
    Point top    {  0,  1};
    
    std::vector< Point > points {left, bottom, right, top};
    for (const auto &point : points)
    {
      for (int d=0; d<spaceDim; d++)
      {
        auto coord = circularMap(point, d);
        TEST_FLOATING_EQUALITY(point[d], coord, tol);
      }
    }
    
    // the corners of the "unit" square should each map to ± sqrt(2)/2
    Point bottom_left  { -1, -1};
    Point bottom_right {  1, -1};
    Point top_right    {  1,  1};
    Point top_left     { -1,  1};
    
    points = std::vector<Point> {bottom_left, bottom_right, top_right, top_left};
    
    for (const auto &point : points)
    {
      for (int d=0; d<spaceDim; d++)
      {
        auto coord = circularMap(point, d);
        auto expected_coord = sqrt(0.5) * point[d]; // point[d] is unit length; gets us the right sign
        TEST_FLOATING_EQUALITY(expected_coord, coord, tol);
      }
    }
  }
} // namespace
