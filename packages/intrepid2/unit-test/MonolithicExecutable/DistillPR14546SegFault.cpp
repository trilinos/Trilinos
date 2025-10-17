// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Tests against structured integration facilities - "synthetic" test cases (i.e., no geometry specified).
    \author Nathan V. Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_CellGeometryTestUtils.hpp>
#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_IntegrationTools.hpp>
#include <Intrepid2_Kernels.hpp>
#include <Intrepid2_NodalBasisFamily.hpp>
#include <Intrepid2_TensorArgumentIterator.hpp>
#include <Intrepid2_TestUtils.hpp>

#include <Intrepid2_CellGeometry.hpp>
#include "Intrepid2_Data.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

#include "StructuredIntegrationTests_TagDefs.hpp"
#include "StructuredIntegrationTests_Utils.hpp"

namespace
{
 using namespace Intrepid2;

template<class Scalar, typename DeviceType>
void callIntegrate(const TransformedBasisValues<Scalar,DeviceType> vectorDataLeft,
                   const TensorData<Scalar,DeviceType> cellMeasures, const TransformedBasisValues<Scalar,DeviceType> vectorDataRight,
                   Teuchos::FancyOStream &out, bool &success)
{
  const double relTol = 1e-12;
  const double absTol = 1e-12;
  
  using IntegrationTools = Intrepid2::IntegrationTools<DeviceType>;
  
  auto integralsBaseline  = IntegrationTools::allocateIntegralData(vectorDataLeft, cellMeasures, vectorDataRight);
  auto integralsIntegrate = IntegrationTools::allocateIntegralData(vectorDataLeft, cellMeasures, vectorDataRight);
  
  integrate_baseline(integralsBaseline, vectorDataLeft, cellMeasures, vectorDataRight);
  IntegrationTools::integrate(integralsIntegrate, vectorDataLeft, cellMeasures, vectorDataRight);
}


TEUCHOS_UNIT_TEST( PR14546, Distill14546SegFault )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 2;
  const int numComponentPoints = 1;
  const int numFields1 = 2;
  const int numFields2 = 3;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2;
  const int numFieldsPerFamilyRight = numFields2 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  fieldComponentDataViewHost1(1) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0/2.0;
  fieldComponentDataViewHost2(1) = 1.0/3.0;
  fieldComponentDataViewHost2(2) = 1.0/4.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
    
  TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2});
  TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData2,fieldComponentData2});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  auto identityMatrixViewHost = getHostCopy(identityMatrixView);
  
  for (int d1=0; d1<spaceDim; d1++)
  {
    for (int d2=0; d2<spaceDim; d2++)
    {
      identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
    }
  }
  Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
  
  const int numPoints = numComponentPoints * numComponentPoints;
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 2;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyLeft {TensorData<DataScalar,DeviceType>(),tensorDataLeft};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft, secondFamilyLeft};
  
  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(1));
  
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyRight {TensorData<DataScalar,DeviceType>(),tensorDataRight};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight, secondFamilyRight};
  
  VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(1));
  
  TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  callIntegrate(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

} // anonymous namespace
