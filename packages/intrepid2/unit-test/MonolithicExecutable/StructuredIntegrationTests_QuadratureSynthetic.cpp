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

// #pragma mark StructuredIntegration: QuadratureSynthetic_AxisAlignedPath_Case1
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_AxisAlignedPath_Case1 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  const int numPoints = numComponentPoints * numComponentPoints;
  const int numFields = 1;
  
  Data<DataScalar,DeviceType> unitData(1.0, Kokkos::Array<int,2>{numFields,numComponentPoints});
  TensorData<DataScalar,DeviceType> unitTensorData(std::vector< Data<DataScalar,DeviceType> >{unitData,unitData});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim);
  Kokkos::deep_copy(identityMatrixView, 1.0);
  
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, BLOCK_PLUS_DIAGONAL, BLOCK_PLUS_DIAGONAL};
  
  const int blockPlusDiagonalLastNonDiagonal = -1; // no non-diagonals; -1 is the default value, but I prefer to be explicit
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType, blockPlusDiagonalLastNonDiagonal);
  
  const int numFamilies = 2;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamily  {unitTensorData,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamily {TensorData<DataScalar,DeviceType>(),unitTensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {firstFamily, secondFamily};
  
  VectorData<DataScalar,DeviceType> unitVectorData(vectorComponents);
  
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,unitVectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresCellComponent(1.0, Kokkos::Array<int,1>{numCells});
  Data<DataScalar,DeviceType> constantCellMeasuresPointComponent(1.0, Kokkos::Array<int,1>{numComponentPoints});
  std::vector< Data<DataScalar,DeviceType> > cellMeasuresComponents { constantCellMeasuresCellComponent, constantCellMeasuresPointComponent, constantCellMeasuresPointComponent};
  
  const bool separateFirstComponent = true; // so that constantCellMeasures advertises shape (C,P), instead of folding cell dimension into the point tensor product...
  TensorData<DataScalar,DeviceType> constantCellMeasures(cellMeasuresComponents, separateFirstComponent);
  
  testIntegrateMatchesBaseline(transformedUnitVectorData, constantCellMeasures, transformedUnitVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case1
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case1 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 1;
  const int numPoints = 1;
  const int numFields = 1;
  
  Data<DataScalar,DeviceType> unitData(1.0, Kokkos::Array<int,2>{numFields,numPoints});
  TensorData<DataScalar,DeviceType> unitTensorData(std::vector< Data<DataScalar,DeviceType> >{unitData,unitData});
  
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
  
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 2;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamily  {unitTensorData,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamily {TensorData<DataScalar,DeviceType>(),unitTensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {firstFamily, secondFamily};
  
  VectorData<DataScalar,DeviceType> unitVectorData(vectorComponents);
  
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,unitVectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorData, constantCellMeasures, transformedUnitVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case2
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case2 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  const int numComponentFields = 1;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numComponentFields);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numComponentFields);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 2.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numComponentFields,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,CONSTANT};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  
  TensorData<DataScalar,DeviceType> tensorData(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2});
  
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
  
  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim> family  {TensorData<DataScalar,DeviceType>(), tensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {family};
  
  VectorData<DataScalar,DeviceType> vectorData(vectorComponents);
  
  TransformedBasisValues<DataScalar,DeviceType> transformedVectorData(explicitIdentityMatrix,vectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
//  printFunctor2(  fieldComponentData1, std::cout, "fieldComponentData1"); // (F,P)
//  printFunctor2(  fieldComponentData2, std::cout, "fieldComponentData2"); // (F,P)
//  printFunctor2(    nonzeroTensorData, std::cout, "nonzeroTensorData");   // (F,P)
//  printFunctor2(       zeroTensorData, std::cout, "zeroTensorData");      // (F,P)
//
//  printFunctor3(          vectorData, std::cout, "vectorData");             // (F,P,D)
//  printFunctor2( constantCellMeasures, std::cout, "constantCellMeasures");  // (C,P)
//  printFunctor4(transformedVectorData, std::cout, "transformedVectorData"); // (C,F,P,D)
  
  const int numFields = numComponentFields * numComponentFields;
  TEST_EQUALITY(numFields, vectorData.extent_int(0)); // (F,P,D)
  TEST_EQUALITY(numFields, transformedVectorData.extent_int(1)); // (C,F,P,D)
  
  testIntegrateMatchesBaseline(transformedVectorData, constantCellMeasures, transformedVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case3
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case3 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 2;
  const int numComponentPoints = 2;
  const int numFields = 2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields,numComponentPoints);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0,0) = 1.0;
  fieldComponentDataViewHost1(0,1) = 2.0;
  fieldComponentDataViewHost1(1,0) = 3.0;
  fieldComponentDataViewHost1(1,1) = 4.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields,numComponentPoints);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0,0) = 1.0/1.0;
  fieldComponentDataViewHost2(0,1) = 1.0/2.0;
  fieldComponentDataViewHost2(1,0) = 1.0/3.0;
  fieldComponentDataViewHost2(1,1) = 1.0/4.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
    
  TensorData<DataScalar,DeviceType> tensorData(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2});
  
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
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamily  {tensorData,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamily {TensorData<DataScalar,DeviceType>(),tensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {firstFamily, secondFamily};
  
  VectorData<DataScalar,DeviceType> vectorData(vectorComponents);
  
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,vectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorData, constantCellMeasures, transformedUnitVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case4
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case4 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 2;
  const int numComponentPoints = 1;
  const int numFields1 = 2;
  const int numFields2 = 1;
  const int numFieldsPerFamily = numFields1 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  fieldComponentDataViewHost1(1) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0/2.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
    
  TensorData<DataScalar,DeviceType> nonzeroTensorData(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2});
  
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
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamily  {nonzeroTensorData,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamily {TensorData<DataScalar,DeviceType>(),nonzeroTensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {firstFamily, secondFamily};
  
  VectorData<DataScalar,DeviceType> vectorData(vectorComponents);
  
  TEST_EQUALITY(numFieldsPerFamily, vectorData.numFieldsInFamily(0));
  TEST_EQUALITY(numFieldsPerFamily, vectorData.numFieldsInFamily(1));
  
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,vectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorData, constantCellMeasures, transformedUnitVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case5
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case5 )
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
  
  testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_AxisAlignedPath_Case6_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_AxisAlignedPath_Case6_3D )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  
  const int numFields1 = 1;
  const int numFields2 = 3;
  const int numFields3 = 2;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields3;
  const int numFieldsPerFamilyRight = numFields1 * numFields3 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
//  fieldComponentDataViewHost1(1) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0;
  fieldComponentDataViewHost2(1) = 1.0/3.0;
  fieldComponentDataViewHost2(2) = 1.0/4.0;
  
  auto fieldComponentDataView3 = getFixedRankView<DataScalar>("field component data 3", numFields3);
  auto fieldComponentDataViewHost3 = Kokkos::create_mirror_view(fieldComponentDataView3);
  fieldComponentDataViewHost3(0) = 1.0;
  fieldComponentDataViewHost3(1) = 0.25;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  Kokkos::deep_copy(fieldComponentDataView3, fieldComponentDataViewHost3);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields3;
  Data<DataScalar,DeviceType> fieldComponentData3(fieldComponentDataView3,fieldComponentExtents,fieldComponentVariationTypes);
  
   TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData3});
   TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData3,fieldComponentData2});
   
   const int numPoints = numComponentPoints * numComponentPoints * numComponentPoints;
   
   auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim);
   Kokkos::deep_copy(identityMatrixView, 1.0);
   
   Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
   Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, BLOCK_PLUS_DIAGONAL, BLOCK_PLUS_DIAGONAL};
   
   const int blockPlusDiagonalLastNonDiagonal = -1; // no non-diagonals; -1 is the default value, but I prefer to be explicit
   Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType, blockPlusDiagonalLastNonDiagonal);
  
  // confirm that the matrix is diagonal (this is required to follow the axis-aligned path):
  TEST_EQUALITY(true, explicitIdentityMatrix.isDiagonal());

  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,tensorDataLeft,tensorDataLeft};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};
  
  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,tensorDataRight,tensorDataRight};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
   
   VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
   
   TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
   Data<DataScalar,DeviceType>  constantCellMeasuresCellComponent(1.0, Kokkos::Array<int,1>{numCells});
   Data<DataScalar,DeviceType> constantCellMeasuresPointComponent(1.0, Kokkos::Array<int,1>{numComponentPoints});
   std::vector< Data<DataScalar,DeviceType> > cellMeasuresComponents { constantCellMeasuresCellComponent, constantCellMeasuresPointComponent, constantCellMeasuresPointComponent, constantCellMeasuresPointComponent};
   
   const bool separateFirstComponent = true; // so that constantCellMeasures advertises shape (C,P), instead of folding cell dimension into the point tensor product...
   TensorData<DataScalar,DeviceType> constantCellMeasures(cellMeasuresComponents, separateFirstComponent);
   
   testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case6_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case6_3D )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  
  const int numFields1 = 1;
  const int numFields2 = 3;
  const int numFields3 = 2;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields3;
  const int numFieldsPerFamilyRight = numFields1 * numFields3 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
//  fieldComponentDataViewHost1(1) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0/2.0;
  fieldComponentDataViewHost2(1) = 1.0/3.0;
  fieldComponentDataViewHost2(2) = 1.0/4.0;
  
  auto fieldComponentDataView3 = getFixedRankView<DataScalar>("field component data 3", numFields3);
  auto fieldComponentDataViewHost3 = Kokkos::create_mirror_view(fieldComponentDataView3);
  fieldComponentDataViewHost3(0) = 0.5;
  fieldComponentDataViewHost3(1) = 0.25;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  Kokkos::deep_copy(fieldComponentDataView3, fieldComponentDataViewHost3);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields3;
  Data<DataScalar,DeviceType> fieldComponentData3(fieldComponentDataView3,fieldComponentExtents,fieldComponentVariationTypes);
  
  TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData3});
  TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData3,fieldComponentData2});
   
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
  
  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};
  
  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
   
   VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
   
   TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
   Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
   TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
   
   testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case7_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case7_3D )
{
  // super-simple case for 3D: symmetric data, 1 point, 1 field.
  
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFields3 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields3;
  const int numFieldsPerFamilyRight = numFields1 * numFields3 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0;
  
  auto fieldComponentDataView3 = getFixedRankView<DataScalar>("field component data 3", numFields3);
  auto fieldComponentDataViewHost3 = Kokkos::create_mirror_view(fieldComponentDataView3);
  fieldComponentDataViewHost3(0) = 1.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  Kokkos::deep_copy(fieldComponentDataView3, fieldComponentDataViewHost3);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields3;
  Data<DataScalar,DeviceType> fieldComponentData3(fieldComponentDataView3,fieldComponentExtents,fieldComponentVariationTypes);
  
   TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData3});
   TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData3,fieldComponentData2});
   
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
     
   const int numFamilies = 1;
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};
  
   VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
   TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
   
   VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
   
   TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
   Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
   TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
   
   testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case8_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case8_3D )
{
  // like case 7, but multi-family
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFields3 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields3;
  const int numFieldsPerFamilyRight = numFields1 * numFields3 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0/2.0;
  
  auto fieldComponentDataView3 = getFixedRankView<DataScalar>("field component data 3", numFields3);
  auto fieldComponentDataViewHost3 = Kokkos::create_mirror_view(fieldComponentDataView3);
  fieldComponentDataViewHost3(0) = 0.5;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  Kokkos::deep_copy(fieldComponentDataView3, fieldComponentDataViewHost3);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields3;
  Data<DataScalar,DeviceType> fieldComponentData3(fieldComponentDataView3,fieldComponentExtents,fieldComponentVariationTypes);
  
   TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData3});
   TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData3,fieldComponentData2});
   
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
   
   const int numFamilies = 3;
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyLeft {TensorData<DataScalar,DeviceType>(),tensorDataLeft,TensorData<DataScalar,DeviceType>()};
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > thirdFamilyLeft  {TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>(),tensorDataLeft};
   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft, secondFamilyLeft, thirdFamilyLeft};

   VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
   TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
   TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(1));

   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyRight {TensorData<DataScalar,DeviceType>(),tensorDataRight,TensorData<DataScalar,DeviceType>()};
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > thirdFamilyRight  {TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>(),tensorDataRight};
   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight, secondFamilyRight, thirdFamilyRight};
   VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(1));
  
   TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
   Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
   TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
   
   testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case9_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case9_3D )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 3;
  
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields1;
  const int numFieldsPerFamilyRight = numFields1 * numFields2 * numFields1;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("imitation VALUE data", numFields1, numComponentPoints);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0,0) = 1.0;
  fieldComponentDataViewHost1(0,1) = 2.0;
  fieldComponentDataViewHost1(0,2) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("imitation GRAD data", numFields2, numComponentPoints);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0,0) = 1./1.0;
  fieldComponentDataViewHost2(0,1) = 1./2.0;
  fieldComponentDataViewHost2(0,2) = 1./3.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  
  TensorData<DataScalar,DeviceType>  tensorDataLeft1(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData2,fieldComponentData1,fieldComponentData1});
  TensorData<DataScalar,DeviceType>  tensorDataLeft2(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData1});
  TensorData<DataScalar,DeviceType>  tensorDataLeft3(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData1,fieldComponentData2});
  TensorData<DataScalar,DeviceType> tensorDataRight1 = tensorDataLeft1;
  TensorData<DataScalar,DeviceType> tensorDataRight2 = tensorDataLeft2;
  TensorData<DataScalar,DeviceType> tensorDataRight3 = tensorDataLeft3;
   
  const int numPoints = numComponentPoints * numComponentPoints * numComponentPoints;
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", numPoints, spaceDim, spaceDim);
  auto identityMatrixViewHost = getHostCopy(identityMatrixView);
  
  for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
  {
    for (int d1=0; d1<spaceDim; d1++)
    {
      for (int d2=0; d2<spaceDim; d2++)
      {
        identityMatrixViewHost(pointOrdinal,d1,d2) = (d1 == d2) ? 1.0 : 0.0;
      }
    }
  }
  Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
  
  
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, GENERAL, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft1,tensorDataLeft2,tensorDataLeft3};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};

  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));

  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight1,tensorDataRight2,tensorDataRight3};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
  VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
  
  TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case10_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case10_3D )
{
  // test with variable quadrature weights
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 2;
  
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields1;
  const int numFieldsPerFamilyRight = numFields1 * numFields2 * numFields1;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("imitation VALUE data", numFields1, numComponentPoints);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0,0) = 1;
  fieldComponentDataViewHost1(0,1) = 1;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("imitation GRAD data", numFields2, numComponentPoints);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0,0) = 1;
  fieldComponentDataViewHost2(0,1) = 1;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  
  TensorData<DataScalar,DeviceType>  tensorDataLeft1(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData2,fieldComponentData1,fieldComponentData1});
  TensorData<DataScalar,DeviceType>  tensorDataLeft2(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData1});
  TensorData<DataScalar,DeviceType>  tensorDataLeft3(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData1,fieldComponentData2});
  TensorData<DataScalar,DeviceType> tensorDataRight1 = tensorDataLeft1;
  TensorData<DataScalar,DeviceType> tensorDataRight2 = tensorDataLeft2;
  TensorData<DataScalar,DeviceType> tensorDataRight3 = tensorDataLeft3;
   
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
  
  const int numPoints = numComponentPoints * numComponentPoints * numComponentPoints;
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft1,tensorDataLeft2,tensorDataLeft3};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};

  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));

  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight1,tensorDataRight2,tensorDataRight3};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
  VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
  
  TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
  auto cellMeasures = getFixedRankView<DataScalar>("cellMeasures", numCells, numPoints);
  
  auto cellMeasuresHost = getHostCopy(cellMeasures);
  
  for (int cellOrdinal=0; cellOrdinal<numCells; cellOrdinal++)
  {
    for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
    {
      cellMeasuresHost(cellOrdinal,pointOrdinal) = (pointOrdinal == 0) ? 1 : 0; //1.0 / (pointOrdinal + 1.0);
    }
  }
  Kokkos::deep_copy(cellMeasures, cellMeasuresHost);
  
  Data<DataScalar,DeviceType> cellMeasuresData(cellMeasures, Kokkos::Array<int,2>{numCells,numPoints}, Kokkos::Array<DataVariationType,2>{GENERAL,GENERAL});
  TensorData<DataScalar,DeviceType> cellMeasuresTensorData(cellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, cellMeasuresTensorData, transformedUnitVectorDataRight, out, success);
}

} // anonymous namespace
