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

/** \file   HostCopyTests.cpp
    \brief  Tests against various copy-like constructors that create a host-accessible version of a device-accessible container.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_BasisValues.hpp"
#include "Intrepid2_Data.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedVectorData.hpp"
#include "Intrepid2_VectorData.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <Kokkos_Core.hpp>

using namespace Intrepid2;

namespace
{
  using namespace Intrepid2;

// #pragma mark HostCopy group
  TEUCHOS_UNIT_TEST(HostCopy, BasisValues)
  {
    // first, construct a BasisValues object with scalar tensor data
    
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using Scalar = double;
    const double value = 3.0;
    Data<Scalar,ExecutionSpace> data(value, Kokkos::Array<int,2>{1,1}); // (F,P)
    
    std::vector< Data<Scalar,ExecutionSpace> > tensorComponents {data, data};
    TensorData<Scalar,ExecutionSpace> tensorData(tensorComponents);
    
    BasisValues<Scalar,ExecutionSpace> basisValues(tensorData);
    
    using HostExecSpace = Kokkos::HostSpace::execution_space;
    BasisValues<Scalar,HostExecSpace> basisValuesHost(basisValues);
    
    TEST_EQUALITY(value*value, basisValuesHost(0,0));
    
    // now, construct BasisValues with vector data
    TensorData<Scalar,ExecutionSpace> emptyTensorData; // represents 0 in the vector
    std::vector< TensorData<Scalar,ExecutionSpace> > vectorComponents { emptyTensorData, tensorData };
    std::vector< std::vector< TensorData<Scalar,ExecutionSpace> > > families { vectorComponents };
    
    VectorData<Scalar,ExecutionSpace> vectorData(families);
    
    basisValues      = BasisValues<Scalar,ExecutionSpace>(vectorData);
    basisValuesHost  = BasisValues<Scalar,HostExecSpace>(basisValues);
    
    TEST_EQUALITY(0.0,           basisValuesHost(0,0,0) ); // zero first component
    TEST_EQUALITY(value * value, basisValuesHost(0,0,1) ); // second component is data * data
  }

  TEUCHOS_UNIT_TEST(HostCopy, Data)
  {
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using Scalar = double;
    const double value = 3.0;
    Data<Scalar,ExecutionSpace> data(value, Kokkos::Array<int,1>{1});
    
    using HostExecSpace = Kokkos::HostSpace::execution_space;
    Data<Scalar,HostExecSpace> dataHost(data);
    
    TEST_EQUALITY(value, dataHost(0));
    
    // now, test that the host-copy of an invalid data object is invalid, too
    data = Data<Scalar,ExecutionSpace>(); // empty/invalid Data
    dataHost = Data<Scalar,HostExecSpace>(data);

    TEST_EQUALITY(data.isValid(), dataHost.isValid());
  }

  TEUCHOS_UNIT_TEST(HostCopy, TensorData)
  {
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using Scalar = double;
    const double value = 3.0;
    Data<Scalar,ExecutionSpace> data(value, Kokkos::Array<int,1>{1});
    
    std::vector< Data<Scalar,ExecutionSpace> > tensorComponents {data, data};
    TensorData<Scalar,ExecutionSpace> tensorData(tensorComponents);
    
    using HostExecSpace = Kokkos::HostSpace::execution_space;
    TensorData<Scalar,HostExecSpace> tensorDataHost(tensorData);
    
    TEST_EQUALITY(value*value, tensorDataHost(0));

    // now, test that the host-copy of an invalid TensorData object is invalid, too
    tensorData = TensorData<Scalar,ExecutionSpace>(); // empty/invalid TensorData
    tensorDataHost = TensorData<Scalar,HostExecSpace>(tensorData);

    TEST_EQUALITY(tensorData.isValid(), tensorDataHost.isValid());
  }

  TEUCHOS_UNIT_TEST(HostCopy, TensorPoints)
  {
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using HostExecSpace = Kokkos::HostSpace::execution_space;
    using PointScalar = double;
    using WeightScalar = double;
    DefaultCubatureFactory cub_factory;
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >());
    auto cellTopoKey = cellTopo.getKey();
    const int quadratureDegree = 3;
    auto quadrature = cub_factory.create<ExecutionSpace, PointScalar, WeightScalar>(cellTopoKey, quadratureDegree);
    ordinal_type numRefPoints = quadrature->getNumPoints();
    const int spaceDim = cellTopo.getDimension();
    ViewType<PointScalar> points = ViewType<PointScalar>("quadrature points ref cell", numRefPoints, spaceDim);
    ViewType<WeightScalar> weights = ViewType<WeightScalar>("quadrature weights ref cell", numRefPoints);
    quadrature->getCubature(points, weights);
    
    TensorPoints<PointScalar,ExecutionSpace> tensorPoints;
    TensorData<WeightScalar,ExecutionSpace>  tensorWeights;
    
    using CubatureTensorType = CubatureTensor<ExecutionSpace,PointScalar,WeightScalar>;
    CubatureTensorType* tensorQuadrature = dynamic_cast<CubatureTensorType*>(quadrature.get());

    TEST_ASSERT(tensorQuadrature != NULL);
    
    if (tensorQuadrature)
    {
      tensorPoints  = tensorQuadrature->allocateCubaturePoints();
      tensorWeights = tensorQuadrature->allocateCubatureWeights();
      tensorQuadrature->getCubature(tensorPoints, tensorWeights);
    }
    
    // copy everything to host
    auto pointsHost = Kokkos::create_mirror_view_and_copy(HostExecSpace::memory_space(), points);
    
    // this copy[-like] constructor is the one that's actually under test:
    TensorPoints<PointScalar,HostExecSpace> tensorPointsHost(tensorPoints);
    
    for (int pointOrdinal=0; pointOrdinal<numRefPoints; pointOrdinal++)
    {
      for (int d=0; d<spaceDim; d++)
      {
        TEST_EQUALITY(pointsHost(pointOrdinal,d), tensorPointsHost(pointOrdinal,d));
      }
    }
  }

  TEUCHOS_UNIT_TEST(HostCopy, TransformedVectorData)
  {
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using Scalar = double;
    const double value = 3.0;
    Data<Scalar,ExecutionSpace> data(value, Kokkos::Array<int,2>{1,1}); // F,P (but just 1 for each, here)
    
    std::vector< Data<Scalar,ExecutionSpace> > tensorComponents {data, data};
    TensorData<Scalar,ExecutionSpace> tensorData(tensorComponents);

    // one family with vector entries of form (0, value*value)
    TensorData<Scalar,ExecutionSpace> emptyTensorData; // represents 0 in the vector
    std::vector< TensorData<Scalar,ExecutionSpace> > vectorComponents { emptyTensorData, tensorData };
    std::vector< std::vector< TensorData<Scalar,ExecutionSpace> > > families { vectorComponents };
    
    VectorData<Scalar,ExecutionSpace> vectorData(families);

    // set up a simple diagonal scaling
    const double scaling = 0.1;
    ViewType<Scalar> scalingView("scaling", 2);
    Kokkos::deep_copy(scalingView, scaling);
    const int rank = 4; // (C,P,D,D)
    const int spaceDim = 2;
    Kokkos::Array<int,7> extents {1,1,spaceDim,spaceDim,1,1,1};
    Kokkos::Array<DataVariationType,7> variationTypes {CONSTANT,CONSTANT,BLOCK_PLUS_DIAGONAL,BLOCK_PLUS_DIAGONAL,CONSTANT,CONSTANT,CONSTANT};
    const int blockPlusDiagonalLastNonDiagonal = -1; // only diagonal
    Data<Scalar,ExecutionSpace> transform(scalingView,rank,extents,variationTypes,blockPlusDiagonalLastNonDiagonal);
    
    TransformedVectorData<Scalar,ExecutionSpace> transformedVectorData(transform,vectorData);
    
    using HostExecSpace = Kokkos::HostSpace::execution_space;
    TransformedVectorData<Scalar,HostExecSpace> transformedVectorDataHost(transformedVectorData);

    TEST_EQUALITY(transformedVectorData.rank(), transformedVectorDataHost.rank());
    TEST_EQUALITY(transformedVectorData.extent_int(0), transformedVectorDataHost.extent_int(0));
    TEST_EQUALITY(transformedVectorData.extent_int(1), transformedVectorDataHost.extent_int(1));
    TEST_EQUALITY(transformedVectorData.extent_int(2), transformedVectorDataHost.extent_int(2));
    TEST_EQUALITY(transformedVectorData.extent_int(3), transformedVectorDataHost.extent_int(3));
    
    TEST_EQUALITY(0.0, transformedVectorDataHost(0,0,0,0) ); // zero first component
    const double tol = 1e-13;
    TEST_FLOATING_EQUALITY(value * value * scaling, transformedVectorDataHost(0,0,0,1), tol); // second component is data * data * scaling
  }

  TEUCHOS_UNIT_TEST(HostCopy, VectorData)
  {
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using Scalar = double;
    const double value = 3.0;
    Data<Scalar,ExecutionSpace> data(value, Kokkos::Array<int,2>{1,1}); // F,P (but just 1 for each, here)
    
    std::vector< Data<Scalar,ExecutionSpace> > tensorComponents {data, data};
    TensorData<Scalar,ExecutionSpace> tensorData(tensorComponents);
    
    TensorData<Scalar,ExecutionSpace> emptyTensorData; // represents 0 in the vector
    std::vector< TensorData<Scalar,ExecutionSpace> > vectorComponents { emptyTensorData, tensorData };
    std::vector< std::vector< TensorData<Scalar,ExecutionSpace> > > families { vectorComponents };
    
    VectorData<Scalar,ExecutionSpace> vectorData(families);
    
    using HostExecSpace = Kokkos::HostSpace::execution_space;
    VectorData<Scalar,HostExecSpace> vectorDataHost(vectorData);
    
    TEST_EQUALITY(0.0,           vectorDataHost(0,0,0) ); // zero first component
    TEST_EQUALITY(value * value, vectorDataHost(0,0,1) ); // second component is data * data
  }
} // namespace
