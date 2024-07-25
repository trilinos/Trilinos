// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   BasisValuesTests.cpp
    \brief  Tests to verify that the version of Basis::getValues() that takes the new BasisValues and TensorPoints arguments produces the same results as the one that takes raw Kokkos DynRankViews as arguments.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <Kokkos_Core.hpp>

using namespace Intrepid2;

namespace
{
  using namespace Intrepid2;

  void testTensorPointCubature(const shards::CellTopology &cellTopo, const int &quadratureDegree,
                               const double &relTol, const double &absTol, Teuchos::FancyOStream &out, bool &success)
  {
    using DeviceType = DefaultTestDeviceType;
    using PointScalar = double;
    using WeightScalar = double;
    using CubatureType   = Cubature<DeviceType,PointScalar,WeightScalar>;
    using PointViewType  = typename CubatureType::PointViewTypeAllocatable;
    using WeightViewType = typename CubatureType::WeightViewTypeAllocatable;
    
    DefaultCubatureFactory cub_factory;
    auto cellTopoKey = cellTopo.getKey();
    auto quadrature = cub_factory.create<DeviceType, PointScalar, WeightScalar>(cellTopoKey, quadratureDegree);
    ordinal_type numRefPoints = quadrature->getNumPoints();
    const int spaceDim = cellTopo.getDimension();
    
    PointViewType points("quadrature points ref cell", numRefPoints, spaceDim);
    WeightViewType weights("quadrature weights ref cell", numRefPoints);
    
    quadrature->getCubature(points, weights);
    
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
    
    printView(points, out, "Points being tested");
    printFunctor2(tensorPoints, out, "tensorPoints");
    
    // points and tensorPoints should be identical, no roundoff: use 0.0 tolerances
    testFloatingEquality2(points,  tensorPoints,    0.0,    0.0, out, success, "points",  "tensorPoints");
    testFloatingEquality1(weights,tensorWeights, relTol, absTol, out, success, "weights", "tensorWeights");
  }

  TEUCHOS_UNIT_TEST( CubatureTensor, Line )
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >());
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      testTensorPointCubature(cellTopo, polyOrder, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CubatureTensor, Quadrilateral )
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >());
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      testTensorPointCubature(cellTopo, polyOrder, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CubatureTensor, Hexahedron )
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >());

    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      testTensorPointCubature(cellTopo, polyOrder, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CubatureTensor, Tetrahedron )
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<> >());
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      testTensorPointCubature(cellTopo, polyOrder, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CubatureTensor, Triangle )
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<> >());
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      testTensorPointCubature(cellTopo, polyOrder, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( CubatureTensor, Wedge )
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Wedge<> >());
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      testTensorPointCubature(cellTopo, polyOrder, relTol, absTol, out, success);
    }
  }
} // namespace
