// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   CellGeometryTests.cpp
    \brief  Tests against CellGeometry
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_CellGeometryTestUtils.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ProjectedGeometry.hpp"
#include "Intrepid2_ProjectedGeometryExamples.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <Kokkos_Core.hpp>

using namespace Intrepid2;

namespace
{
  template< typename PointScalar, int spaceDim, typename DeviceType >
  void testJacobiansAgree(Intrepid2::CellGeometry<PointScalar,spaceDim,DeviceType> cellGeometry,
                          const double &relTol, const double &absTol, Teuchos::FancyOStream &out, bool &success)
  {
    // copy the nodes from CellGeometry into a raw View
    const int numCells = cellGeometry.extent_int(0);
    const int numNodes = cellGeometry.extent_int(1);
    ScalarView<PointScalar,DeviceType> cellNodes("cell nodes", numCells, numNodes, spaceDim);

    using ExecutionSpace = typename DeviceType::execution_space;
    auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,numNodes,spaceDim});
    
    // "workset"
    Kokkos::parallel_for("copy cell nodes from CellGeometry", policy,
    KOKKOS_LAMBDA (const int &cellOrdinal, const int &nodeOrdinal, const int &d) {
      cellNodes(cellOrdinal,nodeOrdinal,d) = cellGeometry(cellOrdinal,nodeOrdinal,d);
    });
    
    ExecutionSpace().fence();
    
    printView(cellNodes, out, "cellNodes");
    
    const int quadratureDegree = 4;
    DefaultCubatureFactory cub_factory;
    auto cellTopoKey = cellGeometry.cellTopology().getKey();
    auto quadrature = cub_factory.create<DeviceType, PointScalar, PointScalar>(cellTopoKey, quadratureDegree);
    ordinal_type numRefPoints = quadrature->getNumPoints();
    using WeightScalar = PointScalar;
    auto points  = getView<PointScalar,DeviceType>("quadrature points ref cell", numRefPoints, spaceDim);
    auto weights = getView<WeightScalar,DeviceType>("quadrature weights ref cell", numRefPoints);
    quadrature->getCubature(points, weights);
    
    auto basisForNodes = cellGeometry.basisForNodes();
    out << std::endl << "basisForNodes: " << basisForNodes->getName() << std::endl;

    ScalarView<PointScalar,DeviceType> expectedJacobians("cell Jacobians", numCells, numRefPoints, spaceDim, spaceDim);
    CellTools<DeviceType>::setJacobian(expectedJacobians, points, cellNodes, basisForNodes);
    
    TensorPoints<PointScalar,DeviceType> tensorPoints;
    TensorData<WeightScalar,DeviceType>  tensorWeights;
    
    using CubatureTensorType = CubatureTensor<DeviceType,PointScalar,WeightScalar>;
    CubatureTensorType* tensorQuadrature = dynamic_cast<CubatureTensorType*>(quadrature.get());

    if (tensorQuadrature)
    {
      tensorPoints  = tensorQuadrature->allocateCubaturePoints();
      tensorWeights = tensorQuadrature->allocateCubatureWeights();
      tensorQuadrature->getCubature(tensorPoints, tensorWeights);
    }
    else
    {
      std::vector<ViewType<PointScalar,DeviceType>> pointComponents {points};
      tensorPoints = TensorPoints<PointScalar,DeviceType>(pointComponents);
      Data<WeightScalar,DeviceType> weightData(weights);
      std::vector<Data<WeightScalar,DeviceType>> weightComponents {weightData};
      tensorWeights = TensorData<WeightScalar,DeviceType>(weightComponents);
    }
    
    Data<PointScalar,DeviceType> jacobians = cellGeometry.allocateJacobianData(tensorPoints);
    auto refData = cellGeometry.getJacobianRefData(tensorPoints);
    cellGeometry.setJacobian(jacobians, tensorPoints, refData);
    
    printFunctor4(jacobians, out, "jacobians from CellGeometry");
    testFloatingEquality4(expectedJacobians, jacobians, relTol, absTol, out, success, "expected Jacobians", "actual Jacobians from CellGeometry");
  }

  template<int spaceDim, class PointScalar>
  void testJacobians(const int &polyOrder, const int &meshWidth, Teuchos::FancyOStream &out, bool &success)
  {
    // tests Jacobian, Jacobian inverse, and Jacobian determinant on a uniform mesh
    
    const double relTol = 1e-12;
    const double absTol = 1e-12;
      
    using DeviceType = DefaultTestDeviceType;
    CellGeometry<PointScalar,spaceDim,DeviceType> cellNodes = uniformCartesianMesh<PointScalar,spaceDim,DeviceType>(1.0, meshWidth);
    
    shards::CellTopology cellTopo;
    if      (spaceDim == 1) cellTopo = shards::getCellTopologyData< shards::Line<>          >();
    else if (spaceDim == 2) cellTopo = shards::getCellTopologyData< shards::Quadrilateral<> >();
    else if (spaceDim == 3) cellTopo = shards::getCellTopologyData< shards::Hexahedron<>    >();
    
    auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
    const int numPoints = cubature->getNumPoints();
    
    // now, compute Jacobian values in the classic, expanded way
    ScalarView<double,DeviceType> cubaturePoints("cubature points",numPoints,spaceDim);
    ScalarView<double,DeviceType> cubatureWeights("cubature weights", numPoints);
    
    cubature->getCubature(cubaturePoints, cubatureWeights);
    
    const int numCells        = cellNodes.numCells();
    const int numNodesPerCell = cellNodes.numNodesPerCell();
    auto expandedCellNodes    = getView<PointScalar,DeviceType>("expanded cell nodes",numCells,numNodesPerCell,spaceDim);
    using ExecutionSpace = typename DeviceType::execution_space;
    auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,numNodesPerCell});
    Kokkos::parallel_for("fill expanded cell nodes", policy,
    KOKKOS_LAMBDA (const int &cellOrdinal, const int &nodeOrdinal)
    {
      for (int d=0; d<spaceDim; d++)
      {
        expandedCellNodes(cellOrdinal,nodeOrdinal,d) = cellNodes(cellOrdinal,nodeOrdinal,d);
      }
    });
    
    auto cellMeasures                = getView<PointScalar,DeviceType>("cell measures",        numCells, numPoints);
    auto expandedJacobianDeterminant = getView<PointScalar,DeviceType>("jacobian determinant", numCells, numPoints);
    auto expandedJacobian            = getView<PointScalar,DeviceType>("jacobian",             numCells, numPoints, spaceDim, spaceDim);
    auto expandedJacobianInverse     = getView<PointScalar,DeviceType>("jacobian inverse",     numCells, numPoints, spaceDim, spaceDim);
    
    using CellTools = Intrepid2::CellTools<DeviceType>;
    
    CellTools::setJacobian(expandedJacobian, cubaturePoints, expandedCellNodes, cellTopo);
    CellTools::setJacobianInv(expandedJacobianInverse, expandedJacobian);
    CellTools::setJacobianDet(expandedJacobianDeterminant, expandedJacobian);
    
    auto jacobian = cellNodes.allocateJacobianData(numPoints);
    auto jacobianDet = CellTools::allocateJacobianDet(jacobian);
    auto jacobianInv = CellTools::allocateJacobianInv(jacobian);
    cellNodes.setJacobian    (jacobian,  numPoints);
    CellTools::setJacobianDet(jacobianDet, jacobian);
    CellTools::setJacobianInv(jacobianInv, jacobian);
      
    //  printView(jacobianDet.getUnderlyingView(), std::cout, "new Jacobian det. data");
    
    testFloatingEquality4(expandedJacobian,            jacobian,    relTol, absTol, out, success, "standard Jacobian", "new Jacobian");
    testFloatingEquality4(expandedJacobianInverse,     jacobianInv, relTol, absTol, out, success, "standard Jacobian inverse", "new Jacobian inverse");
    testFloatingEquality2(expandedJacobianDeterminant, jacobianDet, relTol, absTol, out, success, "standard Jacobian det.", "new Jacobian det.");
    
    // test again, but now force CellGeometry to expand things
    // (this exercises a different code path within CellGeometry)
    bool copyAffineness = false;
    CellGeometry<PointScalar,spaceDim,DeviceType> nonAffineCellGeometry = getNodalCellGeometry(cellNodes, copyAffineness);
    
    jacobian = nonAffineCellGeometry.allocateJacobianData(cubaturePoints);
    jacobianDet = CellTools::allocateJacobianDet(jacobian);
    jacobianInv = CellTools::allocateJacobianInv(jacobian);
    
    auto refData = nonAffineCellGeometry.getJacobianRefData(cubaturePoints);
    nonAffineCellGeometry.setJacobian(jacobian, cubaturePoints, refData);
    CellTools::setJacobianDet(jacobianDet, jacobian);
    CellTools::setJacobianInv(jacobianInv, jacobian);
    
    testFloatingEquality4(expandedJacobian,            jacobian,    relTol, absTol, out, success, "standard Jacobian", "new Jacobian (non-affine path)");
    testFloatingEquality4(expandedJacobianInverse,     jacobianInv, relTol, absTol, out, success, "standard Jacobian inverse", "new Jacobian inverse  (non-affine path)");
    testFloatingEquality2(expandedJacobianDeterminant, jacobianDet, relTol, absTol, out, success, "standard Jacobian det.", "new Jacobian det.  (non-affine path)");
  }

  TEUCHOS_UNIT_TEST( CellGeometry, AllocateCellMeasure_AffineQuad )
  {
    const int spaceDim = 2;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,1};
    Kokkos::Array<int,spaceDim> gridCellCounts{2,2};
    
    using Geometry = CellGeometry<PointScalar, spaceDim, DeviceType>;
    Geometry::SubdivisionStrategy   subdivisionStrategy = Geometry::NO_SUBDIVISION;
    Geometry::HypercubeNodeOrdering nodeOrdering        = Geometry::HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS;
    Geometry cellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy, nodeOrdering);
    
    const int polyOrder = 3;
    
    auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(cellGeometry.cellTopology(),polyOrder*2);
    auto tensorCubatureWeights = cubature->allocateCubatureWeights();
    auto tensorCubaturePoints  = cubature->allocateCubaturePoints();
    
    const int numCells  = cellGeometry.numCells();
    const int numPoints = tensorCubaturePoints.extent_int(0);
    
    Data<PointScalar,DeviceType> jacobian = cellGeometry.allocateJacobianData(tensorCubaturePoints);
    Data<PointScalar,DeviceType> jacobianDet = CellTools<DeviceType>::allocateJacobianDet(jacobian);
    
    // jacobian should have shape (C,P,D,D)
    TEST_EQUALITY(4,         jacobian.rank());
    TEST_EQUALITY(numCells,  jacobian.extent_int(0));
    TEST_EQUALITY(numPoints, jacobian.extent_int(1));
    TEST_EQUALITY(spaceDim,  jacobian.extent_int(2));
    TEST_EQUALITY(spaceDim,  jacobian.extent_int(3));
    
    // jacobianDet should have shape (C,P)
    TEST_EQUALITY(2,         jacobianDet.rank());
    TEST_EQUALITY(numCells,  jacobianDet.extent_int(0));
    TEST_EQUALITY(numPoints, jacobianDet.extent_int(1));
    
    cubature->getCubature(tensorCubaturePoints, tensorCubatureWeights);
    
    TensorData<PointScalar,DeviceType> cellMeasures = cellGeometry.allocateCellMeasure(jacobianDet, tensorCubatureWeights);
    
    // cellMeasures should have shape (C,P)
    TEST_EQUALITY(2,         cellMeasures.rank());
    TEST_EQUALITY(numCells,  cellMeasures.extent_int(0));
    TEST_EQUALITY(numPoints, cellMeasures.extent_int(1));
    
    // given that cellGeometry is affine, cellMeasures should have non-trivial tensor structure, with a separated cell dimension:
    const int expectedNumTensorComponents = spaceDim + 1;
    TEST_EQUALITY(expectedNumTensorComponents, cellMeasures.numTensorComponents());
    TEST_EQUALITY(true, cellMeasures.separateFirstComponent());
  }

  TEUCHOS_UNIT_TEST( CellGeometry, AllocateCellMeasure_GeneralQuad )
  {
    const int spaceDim = 2;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,1};
    Kokkos::Array<int,spaceDim> gridCellCounts{2,2};
    
    using Geometry = CellGeometry<PointScalar, spaceDim, DeviceType>;
    Geometry::SubdivisionStrategy   subdivisionStrategy = Geometry::NO_SUBDIVISION;
    Geometry::HypercubeNodeOrdering nodeOrdering        = Geometry::HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS;
    Geometry uniformCellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy, nodeOrdering);
    
    const bool copyAffineness = false;
    Geometry cellGeometry = getNodalCellGeometry(uniformCellGeometry, copyAffineness);
    
    const int polyOrder = 3;
    
    auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(cellGeometry.cellTopology(),polyOrder*2);
    auto tensorCubatureWeights = cubature->allocateCubatureWeights();
    auto tensorCubaturePoints  = cubature->allocateCubaturePoints();
    
    const int numCells  = cellGeometry.numCells();
    const int numPoints = tensorCubaturePoints.extent_int(0);
    
    Data<PointScalar,DeviceType> jacobian = cellGeometry.allocateJacobianData(tensorCubaturePoints);
    Data<PointScalar,DeviceType> jacobianDet = CellTools<DeviceType>::allocateJacobianDet(jacobian);
    
    // jacobian should have shape (C,P,D,D)
    TEST_EQUALITY(4,         jacobian.rank());
    TEST_EQUALITY(numCells,  jacobian.extent_int(0));
    TEST_EQUALITY(numPoints, jacobian.extent_int(1));
    TEST_EQUALITY(spaceDim,  jacobian.extent_int(2));
    TEST_EQUALITY(spaceDim,  jacobian.extent_int(3));
    
    // jacobianDet should have shape (C,P)
    TEST_EQUALITY(2,         jacobianDet.rank());
    TEST_EQUALITY(numCells,  jacobianDet.extent_int(0));
    TEST_EQUALITY(numPoints, jacobianDet.extent_int(1));
    
    cubature->getCubature(tensorCubaturePoints, tensorCubatureWeights);
    
    TensorData<PointScalar,DeviceType> cellMeasures = cellGeometry.allocateCellMeasure(jacobianDet, tensorCubatureWeights);
    
    // cellMeasures should have shape (C,P)
    TEST_EQUALITY(2,         cellMeasures.rank());
    TEST_EQUALITY(numCells,  cellMeasures.extent_int(0));
    TEST_EQUALITY(numPoints, cellMeasures.extent_int(1));
    
    // given that cellGeometry is not affine, cellMeasures should have trivial tensor structure, and should not have separated first (cell) component:
    const int expectedNumTensorComponents = 1;
    TEST_EQUALITY(expectedNumTensorComponents, cellMeasures.numTensorComponents());
    TEST_EQUALITY(false, cellMeasures.separateFirstComponent());
  }

  TEUCHOS_UNIT_TEST( CellGeometry, CellNodesOrder_Uniform_Quads )
  {
    const int spaceDim = 2;
    const double tol=1e-15;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    // unit quad should be defined as shards does: counter-clockwise, starting at the origin vertex
    Kokkos::Array<PointScalar,spaceDim> origin{0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,1};
    Kokkos::Array<int,spaceDim> gridCellCounts{1,1};
    
    using Geometry = CellGeometry<PointScalar, spaceDim, DeviceType>;
    Geometry::SubdivisionStrategy   subdivisionStrategy = Geometry::NO_SUBDIVISION;
    Geometry::HypercubeNodeOrdering nodeOrdering        = Geometry::HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS;
    Geometry cellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy, nodeOrdering);
    
    const int cellOrdinal = 0;
    const int numNodes = 4;
    
    using std::vector;
    using Point = vector<PointScalar>;
    vector< Point > expectedNodes {Point{0,0}, Point{1,0}, Point{1,1}, Point{0,1}};
    for (int node=0; node<numNodes; node++)
    {
      for (int d=0; d<spaceDim; d++)
      {
        TEST_FLOATING_EQUALITY(expectedNodes[node][d], cellGeometry(cellOrdinal,node,d), tol);
      }
    }
    
    // now, test the tensor node ordering
    expectedNodes = vector< Point >{Point{0,0}, Point{1,0}, Point{0,1}, Point{1,1}};
    nodeOrdering = Geometry::HYPERCUBE_NODE_ORDER_TENSOR;
    cellGeometry = Geometry(origin, domainExtents, gridCellCounts, subdivisionStrategy, nodeOrdering);
    
    for (int node=0; node<numNodes; node++)
    {
      for (int d=0; d<spaceDim; d++)
      {
        TEST_FLOATING_EQUALITY(expectedNodes[node][d], cellGeometry(cellOrdinal,node,d), tol);
      }
    }
  }

  TEUCHOS_UNIT_TEST( CellGeometry, NodeBased_Constructor )
  {
    const int spaceDim = 3;
  
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    using ExecutionSpace = DeviceType::execution_space;
    
    const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >());

    using Geometry = CellGeometry<PointScalar, spaceDim, DeviceType>;
    Geometry::HypercubeNodeOrdering nodeOrdering = Geometry::HYPERCUBE_NODE_ORDER_TENSOR;

    ScalarView<int,DeviceType> cellToNodes;
    ScalarView<PointScalar,DeviceType> nodes("nodes", 1, 8, spaceDim);
    Kokkos::parallel_for("fill nodes for CellGeometry", Kokkos::RangePolicy<ExecutionSpace>(0,nodes.extent(1)), 
    KOKKOS_LAMBDA (const int i) {
      nodes(0,i,0) = (i%2 == 1) ?     1. : 0.;
      nodes(0,i,1) = ((i/2)%2 == 1) ? 1. : 0.;
      nodes(0,i,2) = (i/4 == 1) ?     1. : 0.;
    });

    const bool claimAffine = false;
    CellGeometry<PointScalar, spaceDim, DeviceType> cellGeometry(cellTopo, cellToNodes, nodes, claimAffine, nodeOrdering);
    
    TEST_EQUALITY(1,         cellGeometry.numCells());
    TEST_EQUALITY(8,         cellGeometry.numNodesPerCell());
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Uniform_Pyramids )
  {
    const int spaceDim = 3;
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,2,4};
    Kokkos::Array<int,spaceDim> gridCellCounts{2,4,2};
    
    using CellGeometryType = CellGeometry<PointScalar, spaceDim, DeviceType>;
    
    std::vector<CellGeometryType::SubdivisionStrategy> subdivisionStrategies = {CellGeometryType::SIX_PYRAMIDS};
    for (const auto & subdivisionStrategy : subdivisionStrategies)
    {
      CellGeometry<PointScalar, spaceDim, DeviceType> cellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy);
      testJacobiansAgree(cellGeometry, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Uniform_Quads )
  {
    const int spaceDim = 2;
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,2};
    Kokkos::Array<int,spaceDim> gridCellCounts{4,8};
    CellGeometry<PointScalar, spaceDim, DeviceType>::SubdivisionStrategy subdivisionStrategy = CellGeometry<PointScalar, spaceDim, DeviceType>::NO_SUBDIVISION;
    CellGeometry<PointScalar, spaceDim, DeviceType> cellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy);
    
    testJacobiansAgree(cellGeometry, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Uniform_Hexahedra )
  {
    const int spaceDim = 3;
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,2,3};
    Kokkos::Array<int,spaceDim> gridCellCounts{2,4,8};
    CellGeometry<PointScalar, spaceDim, DeviceType>::SubdivisionStrategy subdivisionStrategy = CellGeometry<PointScalar, spaceDim, DeviceType>::NO_SUBDIVISION;
    CellGeometry<PointScalar, spaceDim, DeviceType> cellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy);
    
    testJacobiansAgree(cellGeometry, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Uniform_Triangles )
  {
    const int spaceDim = 2;
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,2};
    Kokkos::Array<int,spaceDim> gridCellCounts{2,4};
    
    using CellGeometryType = CellGeometry<PointScalar, spaceDim, DeviceType>;
    
    std::vector<CellGeometryType::SubdivisionStrategy> subdivisionStrategies = {CellGeometryType::TWO_TRIANGLES_RIGHT, CellGeometryType::TWO_TRIANGLES_LEFT, CellGeometryType::FOUR_TRIANGLES};
    for (const auto & subdivisionStrategy : subdivisionStrategies)
    {
      CellGeometry<PointScalar, spaceDim, DeviceType> cellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy);
      testJacobiansAgree(cellGeometry, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Uniform_Tetrahedra )
  {
    const int spaceDim = 3;
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,2,4};
    Kokkos::Array<int,spaceDim> gridCellCounts{2,4,2};
    
    using CellGeometryType = CellGeometry<PointScalar, spaceDim, DeviceType>;
    
    std::vector<CellGeometryType::SubdivisionStrategy> subdivisionStrategies = {CellGeometryType::FIVE_TETRAHEDRA}; // TODO: add CellGeometryType::SIX_TETRAHEDRA once CellGeometry supports that
    for (const auto & subdivisionStrategy : subdivisionStrategies)
    {
      CellGeometry<PointScalar, spaceDim, DeviceType> cellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy);
      testJacobiansAgree(cellGeometry, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CellGeometry, JacobiansUniform_1D )
  {
    using PointScalar = double; // it appears that CubatureFactory does not support DFad
    
    const int spaceDim = 1;
    const int polyOrder = 1; // this just affects the number of quadrature points
    const int meshWidth = 10;
    testJacobians<spaceDim, PointScalar>(polyOrder, meshWidth, out, success);
  }

  TEUCHOS_UNIT_TEST( CellGeometry, JacobiansUniform_2D )
  {
    using PointScalar = double; // it appears that CubatureFactory does not support DFad
    
    const int spaceDim = 2;
    const int polyOrder = 1; // this just affects the number of quadrature points
    const int meshWidth = 4;
    testJacobians<spaceDim,PointScalar>(polyOrder, meshWidth, out, success);
  }

  TEUCHOS_UNIT_TEST( CellGeometry, JacobiansUniform_3D )
  {
    using PointScalar = double; // it appears that CubatureFactory does not support DFad
    
    const int spaceDim = 3;
    const int polyOrder = 1; // this just affects the number of quadrature points
    const int meshWidth = 2;
    testJacobians<spaceDim,PointScalar>(polyOrder, meshWidth, out, success);
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Linear_Pyramids )
  {
    const int spaceDim = 3;
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,2,4};
    Kokkos::Array<int,spaceDim> gridCellCounts{2,4,2};
    CellGeometry<PointScalar, spaceDim, DeviceType>::SubdivisionStrategy subdivisionStrategy = CellGeometry<PointScalar, spaceDim, DeviceType>::SIX_PYRAMIDS;
    CellGeometry<PointScalar, spaceDim, DeviceType> uniformGridCellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy);
    
    const std::vector<bool> copyAffineValues {false, true};
    
    for (auto copyAffine : copyAffineValues)
    {
      if (copyAffine) out << "testing affine path.\n";
      else            out << "testing non-affine path (with affine geometry).\n";
      auto nodalCellGeometry = getNodalCellGeometry(uniformGridCellGeometry,copyAffine);
      testJacobiansAgree(nodalCellGeometry, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Linear_Quads )
  {
    const int spaceDim = 2;
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,2};
    Kokkos::Array<int,spaceDim> gridCellCounts{4,8};
    CellGeometry<PointScalar, spaceDim, DeviceType>::SubdivisionStrategy subdivisionStrategy = CellGeometry<PointScalar, spaceDim, DeviceType>::NO_SUBDIVISION;
    CellGeometry<PointScalar, spaceDim, DeviceType> uniformGridCellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy);
    
    const std::vector<bool> copyAffineValues {false, true};
    
    for (auto copyAffine : copyAffineValues)
    {
      if (copyAffine) out << "testing affine path.\n";
      else            out << "testing non-affine path (with affine geometry).\n";
      auto nodalCellGeometry = getNodalCellGeometry(uniformGridCellGeometry,copyAffine);
      testJacobiansAgree(nodalCellGeometry, relTol, absTol, out, success);
    }
    
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Linear_Hexahedra )
  {
    const int spaceDim = 3;
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,2,3};
    Kokkos::Array<int,spaceDim> gridCellCounts{2,4,8};
    CellGeometry<PointScalar, spaceDim, DeviceType>::SubdivisionStrategy subdivisionStrategy = CellGeometry<PointScalar, spaceDim, DeviceType>::NO_SUBDIVISION;
    CellGeometry<PointScalar, spaceDim, DeviceType> uniformGridCellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy);
    
    const std::vector<bool> copyAffineValues {false, true};
    
    for (auto copyAffine : copyAffineValues)
    {
      if (copyAffine) out << "testing affine path.\n";
      else            out << "testing non-affine path (with affine geometry).\n";
      auto nodalCellGeometry = getNodalCellGeometry(uniformGridCellGeometry,copyAffine);
      testJacobiansAgree(nodalCellGeometry, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Linear_Triangles )
  {
    const int spaceDim = 2;
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using PointScalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,2};
    Kokkos::Array<int,spaceDim> gridCellCounts{2,4};
    
    using CellGeometryType = CellGeometry<PointScalar, spaceDim, DeviceType>;
    
    std::vector<CellGeometryType::SubdivisionStrategy> subdivisionStrategies = {CellGeometryType::TWO_TRIANGLES_RIGHT, CellGeometryType::TWO_TRIANGLES_LEFT, CellGeometryType::FOUR_TRIANGLES};
    for (const auto & subdivisionStrategy : subdivisionStrategies)
    {
      CellGeometry<PointScalar, spaceDim, DeviceType> uniformGridCellGeometry(origin, domainExtents, gridCellCounts, subdivisionStrategy);
      
      // "affine" and "non-affine" paths for linear triangles should be identical in terms of e.g. the number of Jacobian evaluations (the affine structure is detected rather than enforced because of the claimAffine argument).  We test both possibilities anyway.
      const std::vector<bool> copyAffineValues {false, true};
      for (auto copyAffine : copyAffineValues)
      {
        if (copyAffine) out << "testing affine path.\n";
        else            out << "testing non-affine path (with affine geometry).\n";
        auto nodalCellGeometry = getNodalCellGeometry(uniformGridCellGeometry,copyAffine);
        testJacobiansAgree(nodalCellGeometry, relTol, absTol, out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( CellGeometry, Nodes_Uniform_Pyramids )
  {
    const int spaceDim = 3;
    using PointScalar = double;
    using HostDeviceType = Kokkos::HostSpace::device_type;
    
    Kokkos::Array<PointScalar,spaceDim> origin{0,0,0};
    Kokkos::Array<PointScalar,spaceDim> domainExtents{1,1,1};
    Kokkos::Array<int,spaceDim> gridCellCounts{1,1,1};
    
    using CellGeometryType = CellGeometry<PointScalar, spaceDim, HostDeviceType>;
    CellGeometryType cellGeometry(origin, domainExtents, gridCellCounts, CellGeometryType::SIX_PYRAMIDS);
    
    ordinal_type numGridCells = gridCellCounts[0] * gridCellCounts[1] * gridCellCounts[2];
    ordinal_type expectedNumCells = 6 * numGridCells;
    ordinal_type numCells = cellGeometry.numCells();
    TEST_EQUALITY(numCells, expectedNumCells);
    
//    ordinal_type expectedNumNodes = 9; // 1 interior vertex, 8 vertices for the hexahedron
    // TODO: check global node numbers.  We expect to have interior nodes listed first (there's just one per grid cell), followed by a Cartesian numbering of the vertices with x as the fastest-moving index.
    
  }

//TODO: write this test, using ProjectedGeometry.
//  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Quadratic_Quads )
//  {
//    success = false;
//    out << "Still need to write this test!\n";
//  }
//
//  TEUCHOS_UNIT_TEST( CellGeometry, Jacobian_Quadratic_Triangles )
//  {
//    // TODO: write this test.  Use ProjectedGeometry.  (Should also write a few tests directly against ProjectedGeometry, not to mention Data::storeMatVec().)
//    success = false;
//    out << "Still need to write this test!\n";
//  }

} // namespace
