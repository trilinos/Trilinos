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
