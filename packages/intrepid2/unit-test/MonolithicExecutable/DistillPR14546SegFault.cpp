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
  using D  = Data<DataScalar,DeviceType>;
  using TD = TensorData<DataScalar,DeviceType>;
  using VD = VectorData<DataScalar,DeviceType>;
  
  const int spaceDim = 2;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2;
  const int numFieldsPerFamilyRight = numFields2 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0/2.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  D fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  D fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
    
  TD  tensorDataLeft(std::vector<D>{fieldComponentData1,fieldComponentData2});
  TD tensorDataRight(std::vector<D>{fieldComponentData2,fieldComponentData2});
  
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
  
  D explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 2;
  TD nullTD;
  Kokkos::Array<TD, spaceDim > firstFamilyLeft  {tensorDataLeft,nullTD};
  Kokkos::Array<TD, spaceDim > secondFamilyLeft {nullTD,tensorDataLeft};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft, secondFamilyLeft};
  
  VD vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(1));
  
  Kokkos::Array<TD, spaceDim > firstFamilyRight  {tensorDataRight,nullTD};
  Kokkos::Array<TD, spaceDim > secondFamilyRight {nullTD,tensorDataRight};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight, secondFamilyRight};
  
  VD vectorDataRight(vectorComponentsRight);
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(1));
  
  TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
  D constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TD constantCellMeasures(constantCellMeasuresData);
  
  callIntegrate(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

} // anonymous namespace
