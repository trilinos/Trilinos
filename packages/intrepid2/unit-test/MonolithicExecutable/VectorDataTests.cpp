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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Tests gainst VectorData class.
    \author Created by Nate Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_Data.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_TestUtils.hpp"

namespace
{
  using namespace Intrepid2;

  //! tests that a manually-constructed VectorData object with tensor-product gradient values matches the one we get from a standard basis.
  template<int spaceDim>
  void testRefSpaceVectorValues(Teuchos::FancyOStream &out, bool &success)
  {
    using Scalar = double;
    using PointScalar = double;
    
    const double relTol = 1e-12;
    const double absTol = 1e-12;
    
    const int polyOrder = 1;
    const int meshWidth = 1;
    
    auto fs = Intrepid2::FUNCTION_SPACE_HGRAD;
    
    auto lineBasis = Intrepid2::getLineBasis< Intrepid2::NodalBasisFamily<> >(fs, polyOrder);
    
    int numFields_1D = lineBasis->getCardinality();
    
    int numFields = 1;
    int numHypercubes = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numHypercubes *= meshWidth;
      numFields     *= numFields_1D;
    }
      
    shards::CellTopology lineTopo = shards::getCellTopologyData< shards::Line<> >();
    shards::CellTopology cellTopo;
    if      (spaceDim == 1) cellTopo = shards::getCellTopologyData< shards::Line<>          >();
    else if (spaceDim == 2) cellTopo = shards::getCellTopologyData< shards::Quadrilateral<> >();
    else if (spaceDim == 3) cellTopo = shards::getCellTopologyData< shards::Hexahedron<>    >();
    
    using ExecSpaceType = Kokkos::DefaultExecutionSpace;
    auto lineCubature = Intrepid2::DefaultCubatureFactory::create<ExecSpaceType>(lineTopo,polyOrder*2);
    int numPoints_1D = lineCubature->getNumPoints();
    ScalarView<PointScalar,ExecSpaceType> lineCubaturePoints("line cubature points",numPoints_1D,1);
    ScalarView<double,ExecSpaceType> lineCubatureWeights("line cubature weights", numPoints_1D);
    
    lineCubature->getCubature(lineCubaturePoints, lineCubatureWeights);
    
    // Allocate some intermediate containers
    ScalarView<Scalar,ExecSpaceType> lineBasisValues    ("line basis values",      numFields_1D, numPoints_1D   );
    ScalarView<Scalar,ExecSpaceType> lineBasisGradValues("line basis grad values", numFields_1D, numPoints_1D, 1);
    
    // for now, we use 1D values to build up the 2D or 3D gradients
    // eventually, TensorBasis should offer a getValues() variant that returns tensor basis data
    lineBasis->getValues(lineBasisValues,     lineCubaturePoints, Intrepid2::OPERATOR_VALUE );
    lineBasis->getValues(lineBasisGradValues, lineCubaturePoints, Intrepid2::OPERATOR_GRAD  );
    
    // drop the trivial space dimension in line gradient values:
    Kokkos::resize(lineBasisGradValues, numFields_1D, numPoints_1D);
      
    Kokkos::Array<TensorData<Scalar,ExecSpaceType>, spaceDim> vectorComponents;
    
    for (int d=0; d<spaceDim; d++)
    {
      Kokkos::Array<Data<Scalar,ExecSpaceType>, spaceDim> gradComponent_d;
      // gradComponent_d stores vector component d of the gradient, expressed as the product of values corresponding to each coordinate dimension
      // The gradient operator is (dx,dy,dz) in 3D; that is, the derivative taken is in the coordinate dimension that matches d.
      // Therefore, the operator leaves the tensorial components in dimension d2≠d unaffected, and results in a 1D "gradient" being taken in the dimension for which d2=d.
      // Hence, the assignment below.
      for (int d2=0; d2<spaceDim; d2++)
      {
        if (d2 == d) gradComponent_d[d2] = Data<Scalar,ExecSpaceType>(lineBasisGradValues);
        else         gradComponent_d[d2] = Data<Scalar,ExecSpaceType>(lineBasisValues);
      }
      vectorComponents[d] = TensorData<Scalar,ExecSpaceType>(gradComponent_d);
    }
    VectorData<Scalar,ExecSpaceType> gradientValues(vectorComponents, false); // false: not axis-aligned
    
    int numPoints = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numPoints *= numPoints_1D;
    }
    
    auto basis = Intrepid2::getBasis< Intrepid2::NodalBasisFamily<> >(cellTopo, fs, polyOrder);
    
    // Allocate some intermediate containers
    ScalarView<Scalar,ExecSpaceType> basisValues    ("basis values", numFields, numPoints );
    ScalarView<Scalar,ExecSpaceType> basisGradValues("basis grad values", numFields, numPoints, spaceDim);

    using Kokkos::DefaultExecutionSpace;
    auto cubature = Intrepid2::DefaultCubatureFactory::create<DefaultExecutionSpace>(cellTopo,polyOrder*2);
    TEST_EQUALITY( numPoints, cubature->getNumPoints());
    ScalarView<PointScalar,ExecSpaceType> cubaturePoints("cubature points",numPoints,spaceDim);
    ScalarView<double,ExecSpaceType> cubatureWeights("cubature weights", numPoints);
    
    cubature->getCubature(cubaturePoints, cubatureWeights);
    
    basis->getValues(basisValues,     cubaturePoints, Intrepid2::OPERATOR_VALUE );
    basis->getValues(basisGradValues, cubaturePoints, Intrepid2::OPERATOR_GRAD  );
    
    testFloatingEquality3(basisGradValues, gradientValues, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( VectorData, RefSpaceVectorValues_1D )
  {
    const int spaceDim = 1;
    testRefSpaceVectorValues<spaceDim>(out, success);
  }

  TEUCHOS_UNIT_TEST( VectorData, RefSpaceVectorValues_2D )
  {
    const int spaceDim = 2;
    testRefSpaceVectorValues<spaceDim>(out, success);
  }

  TEUCHOS_UNIT_TEST( VectorData, RefSpaceVectorValues_3D )
  {
    const int spaceDim = 3;
    testRefSpaceVectorValues<spaceDim>(out, success);
  }
   
  TEUCHOS_UNIT_TEST( VectorData, ZeroFirstComponent )
  {
    using Scalar  = double;
    using ExecSpaceType = Kokkos::DefaultExecutionSpace;
    
    const int spaceDim = 2;
    
    const int numComponentPoints = 1;
    const int numComponentFields = 1;
    const int numFields          = numComponentFields * numComponentFields;
    const int numPoints          = numComponentPoints * numComponentPoints;
    
    ScalarView<Scalar,ExecSpaceType> fieldComponentDataView = getView<Scalar>("field component data", numComponentFields);
    auto fieldComponentDataViewHost = Kokkos::create_mirror_view(fieldComponentDataView);
    fieldComponentDataViewHost(0) = 1.0;
    
    Kokkos::deep_copy(fieldComponentDataView, fieldComponentDataViewHost);
    
    const int fieldComponentDataRank = 2;
    Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numComponentFields,numComponentPoints};
    Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,CONSTANT};
    Data<Scalar,ExecSpaceType> fieldComponentData(fieldComponentDataView,fieldComponentExtents,fieldComponentVariationTypes);
    
    TensorData<Scalar,ExecSpaceType> nonzeroTensorData(std::vector< Data<Scalar,ExecSpaceType> >{fieldComponentData,fieldComponentData});
    
    const int numFamilies = 1;
    Kokkos::Array<TensorData<Scalar,ExecSpaceType>, spaceDim > family {TensorData<Scalar,ExecSpaceType>(), nonzeroTensorData}; // empty first component
    Kokkos::Array< Kokkos::Array<TensorData<Scalar,ExecSpaceType>, spaceDim>, numFamilies> vectorComponents {family};
    
    VectorData<Scalar,ExecSpaceType> vectorData(vectorComponents);
    
    TEST_EQUALITY(numFields, vectorData.extent_int(0)); // (F,P,D)
    TEST_EQUALITY(numPoints, vectorData.extent_int(1)); // (F,P,D)
    TEST_EQUALITY( spaceDim, vectorData.extent_int(2)); // (F,P,D)
    
    // check that the first component is identically zero (indicated by invalidity)
    TEST_EQUALITY( false, vectorData.getComponent(0,0).isValid() ); // getComponent(familyOrdinal, componentOrdinal)
    TEST_EQUALITY(  true, vectorData.getComponent(0,1).isValid() ); // getComponent(familyOrdinal, componentOrdinal)
    
    // test values
    TEST_EQUALITY(                          0.0, vectorData(0,0,0)); // (F,P,D)
    TEST_EQUALITY(fieldComponentDataViewHost(0), vectorData(0,0,1)); // (F,P,D)
  }

} // anonymous namespace
