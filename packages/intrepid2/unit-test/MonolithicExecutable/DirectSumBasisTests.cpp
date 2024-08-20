// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   DirectSumBasisTests.cpp
    \brief  Tests against the DirectSumBasis class.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <Kokkos_Core.hpp>

namespace
{
  using namespace Intrepid2;
// #pragma mark DirectSumBasis group
  TEUCHOS_UNIT_TEST( DirectSumBasis, AllocateBasisValues )
  {
    using Basis = HierarchicalBasisFamily<DefaultTestDeviceType>::HDIV_HEX;
    
    Basis basis(1,1,1);
    
    int quadratureDegree = 1;
    using DeviceType = DefaultTestDeviceType;
    using PointScalar = double;
    using WeightScalar = double;
    DefaultCubatureFactory cub_factory;
    using CubatureType   = Cubature<DeviceType,PointScalar,WeightScalar>;
    using PointViewType  = typename CubatureType::PointViewTypeAllocatable;
    using WeightViewType = typename CubatureType::WeightViewTypeAllocatable;
    auto cellTopo = basis.getBaseCellTopology();
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
    
    // test OPERATOR_DIV
    auto basisValues = basis.allocateBasisValues(tensorPoints, OPERATOR_DIV);
    
    const int basisCardinality = basis.getCardinality();
    TEST_EQUALITY(basisCardinality, basisValues.extent_int(0));
    TEST_EQUALITY(numRefPoints, basisValues.extent_int(1));
    TEST_EQUALITY(2, basisValues.rank());
    
    int family1Fields = basisValues.numFieldsInFamily(0);
    int family2Fields = basisValues.numFieldsInFamily(1);
    int family3Fields = basisValues.numFieldsInFamily(2);
    
    TEST_EQUALITY(basisCardinality, family1Fields + family2Fields + family3Fields);
    
    // test OPERATOR_VALUE
    basisValues = basis.allocateBasisValues(tensorPoints, OPERATOR_VALUE);
    TEST_EQUALITY(basisCardinality, basisValues.extent_int(0));
    TEST_EQUALITY(numRefPoints, basisValues.extent_int(1));
    TEST_EQUALITY(spaceDim, basisValues.extent_int(2));
    TEST_EQUALITY(3, basisValues.rank());
  }
} // namespace
