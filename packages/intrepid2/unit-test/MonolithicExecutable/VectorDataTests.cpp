// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    using PointScalar = double;
    using WeightScalar = double;
    using CubatureType   = Cubature<DeviceType,PointScalar,WeightScalar>;
    using PointViewType  = typename CubatureType::PointViewTypeAllocatable;
    using WeightViewType = typename CubatureType::WeightViewTypeAllocatable;
    
    const double relTol = 1e-12;
    const double absTol = 1e-12;
    
    const int polyOrder = 1;
    const int meshWidth = 1;
    
    auto fs = Intrepid2::FUNCTION_SPACE_HGRAD;
    
    auto lineBasis = Intrepid2::getLineBasis< Intrepid2::NodalBasisFamily<DeviceType> >(fs, polyOrder);
    
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
    
    auto lineCubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(lineTopo,polyOrder*2);
    int numPoints_1D = lineCubature->getNumPoints();
    PointViewType lineCubaturePoints("line cubature points",numPoints_1D,1);
    WeightViewType lineCubatureWeights("line cubature weights", numPoints_1D);
    
    lineCubature->getCubature(lineCubaturePoints, lineCubatureWeights);
    
    // Allocate some intermediate containers
    ScalarView<Scalar,DeviceType> lineBasisValues    ("line basis values",      numFields_1D, numPoints_1D   );
    ScalarView<Scalar,DeviceType> lineBasisGradValues("line basis grad values", numFields_1D, numPoints_1D, 1);
    
    // for now, we use 1D values to build up the 2D or 3D gradients
    // eventually, TensorBasis should offer a getValues() variant that returns tensor basis data
    lineBasis->getValues(lineBasisValues,     lineCubaturePoints, Intrepid2::OPERATOR_VALUE );
    lineBasis->getValues(lineBasisGradValues, lineCubaturePoints, Intrepid2::OPERATOR_GRAD  );
    
    // drop the trivial space dimension in line gradient values:
    Kokkos::resize(lineBasisGradValues, numFields_1D, numPoints_1D);
      
    Kokkos::Array<TensorData<Scalar,DeviceType>, spaceDim> vectorComponents;
    
    for (int d=0; d<spaceDim; d++)
    {
      Kokkos::Array<Data<Scalar,DeviceType>, spaceDim> gradComponent_d;
      // gradComponent_d stores vector component d of the gradient, expressed as the product of values corresponding to each coordinate dimension
      // The gradient operator is (dx,dy,dz) in 3D; that is, the derivative taken is in the coordinate dimension that matches d.
      // Therefore, the operator leaves the tensorial components in dimension d2≠d unaffected, and results in a 1D "gradient" being taken in the dimension for which d2=d.
      // Hence, the assignment below.
      for (int d2=0; d2<spaceDim; d2++)
      {
        if (d2 == d) gradComponent_d[d2] = Data<Scalar,DeviceType>(lineBasisGradValues);
        else         gradComponent_d[d2] = Data<Scalar,DeviceType>(lineBasisValues);
      }
      vectorComponents[d] = TensorData<Scalar,DeviceType>(gradComponent_d);
    }
    VectorData<Scalar,DeviceType> gradientValues(vectorComponents, false); // false: not axis-aligned
    
    int numPoints = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numPoints *= numPoints_1D;
    }
    
    auto basis = Intrepid2::getBasis< Intrepid2::NodalBasisFamily<DeviceType> >(cellTopo, fs, polyOrder);
    
    // Allocate some intermediate containers
    ScalarView<Scalar,DeviceType> basisValues    ("basis values", numFields, numPoints );
    ScalarView<Scalar,DeviceType> basisGradValues("basis grad values", numFields, numPoints, spaceDim);

    auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
    TEST_EQUALITY( numPoints, cubature->getNumPoints());
    PointViewType cubaturePoints("cubature points",numPoints,spaceDim);
    WeightViewType cubatureWeights("cubature weights", numPoints);
    
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
    using DeviceType = Kokkos::DefaultExecutionSpace;
    
    const int spaceDim = 2;
    
    const int numComponentPoints = 1;
    const int numComponentFields = 1;
    const int numFields          = numComponentFields * numComponentFields;
    const int numPoints          = numComponentPoints * numComponentPoints;
    
    ScalarView<Scalar,DeviceType> fieldComponentDataView = getView<Scalar,DeviceType>("field component data", numComponentFields);
    auto fieldComponentDataViewHost = Kokkos::create_mirror_view(fieldComponentDataView);
    fieldComponentDataViewHost(0) = 1.0;
    
    Kokkos::deep_copy(fieldComponentDataView, fieldComponentDataViewHost);
    
    const int fieldComponentDataRank = 2;
    Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numComponentFields,numComponentPoints};
    Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,CONSTANT};
    Data<Scalar,DeviceType> fieldComponentData(fieldComponentDataView,fieldComponentExtents,fieldComponentVariationTypes);
    
    TensorData<Scalar,DeviceType> nonzeroTensorData(std::vector< Data<Scalar,DeviceType> >{fieldComponentData,fieldComponentData});
    
    const int numFamilies = 1;
    Kokkos::Array<TensorData<Scalar,DeviceType>, spaceDim > family {TensorData<Scalar,DeviceType>(), nonzeroTensorData}; // empty first component
    Kokkos::Array< Kokkos::Array<TensorData<Scalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {family};
    
    VectorData<Scalar,DeviceType> vectorData(vectorComponents);
    
    TEST_EQUALITY(numFields, vectorData.extent_int(0)); // (F,P,D)
    TEST_EQUALITY(numPoints, vectorData.extent_int(1)); // (F,P,D)
    TEST_EQUALITY( spaceDim, vectorData.extent_int(2)); // (F,P,D)
    
    Kokkos::View<bool*,DeviceType> vectorDataBools("vectorDataBools", 2);
    Kokkos::View<double*,DeviceType> vectorDataValues("vectorDataValues", 2);

    using HostSpaceType = Kokkos::DefaultHostExecutionSpace;
    Kokkos::parallel_for(Kokkos::RangePolicy<typename DeviceType::execution_space>(0,2),
    KOKKOS_LAMBDA (const int &i) {
      vectorDataBools(i) = vectorData.getComponent(0,i).isValid();
      vectorDataValues(i) = vectorData(0,0,i);
    });
    auto vectorDataBoolsHost = Kokkos::create_mirror_view_and_copy(HostSpaceType(), vectorDataBools);
    auto vectorDataValuesHost = Kokkos::create_mirror_view_and_copy(HostSpaceType(), vectorDataValues);

    // check that the first component is identically zero (indicated by invalidity)
    TEST_EQUALITY( false, vectorDataBoolsHost(0) ); // getComponent(familyOrdinal, componentOrdinal)
    TEST_EQUALITY( true,  vectorDataBoolsHost(1) ); // getComponent(familyOrdinal, componentOrdinal)

    // test values
    TEST_EQUALITY(                          0.0, vectorDataValuesHost(0)); // (F,P,D)
    TEST_EQUALITY(fieldComponentDataViewHost(0), vectorDataValuesHost(1)); // (F,P,D)
  }

} // anonymous namespace
