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
