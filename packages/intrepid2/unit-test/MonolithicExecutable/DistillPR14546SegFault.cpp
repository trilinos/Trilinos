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

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_IntegrationTools.hpp>
#include <Intrepid2_TestUtils.hpp>

#include "Intrepid2_Data.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

namespace
{
 using namespace Intrepid2;

template<class Scalar, typename DeviceType>
void integrate_baseline(Data<Scalar,DeviceType> integrals, const TransformedBasisValues<Scalar,DeviceType> vectorDataLeft,
                      const TensorData<Scalar,DeviceType> cellMeasures, const TransformedBasisValues<Scalar,DeviceType> vectorDataRight)
{
const int spaceDim       = vectorDataLeft.spaceDim();

//  printFunctor4(vectorDataLeft, std::cout, "vectorDataLeft");
//  printFunctor2(cellMeasures, std::cout, "cellMeasures");

// integral data may have shape (C,F1,F2) or (if the variation type is CONSTANT in the cell dimension) shape (F1,F2)
const int integralViewRank = integrals.getUnderlyingViewRank();

}

TEUCHOS_UNIT_TEST( PR14546, Distill14546SegFault )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  using D  = Data<DataScalar,DeviceType>;
  using TD = TensorData<DataScalar,DeviceType>;
  using VD = VectorData<DataScalar,DeviceType>;
  using TBV = TransformedBasisValues<DataScalar,DeviceType>;
  
  const int spaceDim = 1;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2;
  const int numFieldsPerFamilyRight = numFields2 * numFields2;
  
  auto oneElementView = getFixedRankView<DataScalar>("oneElementView", 1);
  Kokkos::deep_copy(oneElementView, 1.0);
  
  Kokkos::Array<int,2> extents {1,1};
  Kokkos::Array<DataVariationType,2> variationTypes {GENERAL,GENERAL};
  D fieldComponentData(oneElementView,extents,variationTypes);

  TD  tensorData(std::vector<D>{fieldComponentData,fieldComponentData});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  Kokkos::deep_copy(identityMatrixView, 1.0);
  
  Kokkos::Array<int,4> transformExtents {1, 1, 1, 1};
  Kokkos::Array<DataVariationType,4> transformationVariationType {GENERAL, GENERAL, GENERAL, GENERAL};
  
  D explicitIdentityMatrix(identityMatrixView, transformExtents, transformationVariationType);
  
  const int numFamilies = 2;
  TD nullTD;
  Kokkos::Array<TD, spaceDim > firstFamilyLeft  {tensorData};
  Kokkos::Array<TD, spaceDim > secondFamilyLeft {tensorData};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft, secondFamilyLeft};
  
  VD vectorDataLeft(vectorComponentsLeft);
  
  Kokkos::Array<TD, spaceDim > firstFamilyRight  {tensorData};
  Kokkos::Array<TD, spaceDim > secondFamilyRight {tensorData};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight, secondFamilyRight};
  
  VD vectorDataRight(vectorComponentsRight);
  
  TBV  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TBV transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataLeft);
  
  D constantCellMeasuresData(1.0, Kokkos::Array<int,2>{1,1});
  TD constantCellMeasures(constantCellMeasuresData);
  
  // these assignments imitate a function call with arguments (tbvLeft, cellMeasures, tbvRight)
  const TBV      tbvLeft = transformedUnitVectorDataLeft;
  const TD  cellMeasures = constantCellMeasures;
  const TBV     tbvRight = transformedUnitVectorDataLeft;
  
  using IT = Intrepid2::IntegrationTools<DeviceType>;
  
  D integralsBaseline  = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  D integralsIntegrate = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  
  integrate_baseline(integralsBaseline, tbvLeft, cellMeasures, tbvRight);
  IT::integrate(integralsIntegrate, tbvLeft, cellMeasures, tbvRight);
}

TEUCHOS_UNIT_TEST( PR14546, AllocationIssue )
{
  // not sure if this is the same issue -- I think so, but in case it isn't, I'm creating a separate test
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  using D  = Data<DataScalar,DeviceType>;
  using TD = TensorData<DataScalar,DeviceType>;
  using VD = VectorData<DataScalar,DeviceType>;
  using TBV = TransformedBasisValues<DataScalar,DeviceType>;
  
  const int spaceDim = 1;
  
  auto oneElementView = getFixedRankView<DataScalar>("oneElementView", 1);
  Kokkos::deep_copy(oneElementView, 1.0);
  
  Kokkos::Array<int,2> extents {1,1};
  Kokkos::Array<DataVariationType,2> variationTypes {GENERAL,GENERAL};
  D fieldComponentData(oneElementView,extents,variationTypes);

  TD tensorData(std::vector<D>{fieldComponentData});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  Kokkos::deep_copy(identityMatrixView, 1.0);
  
  Kokkos::Array<int,4> transformExtents {1, 1, 1, 1};
  Kokkos::Array<DataVariationType,4> transformationVariationType {GENERAL, GENERAL, GENERAL, GENERAL};
  
  D explicitIdentityMatrix(identityMatrixView, transformExtents, transformationVariationType);
  {
    auto data = getMatchingViewWithLabel(identityMatrixView, "Data mat-mat result", 1, 1, 1, 1);
    std::cout << "data.size(): " << data.size() << std::endl;
    std::cout << "Got to line " << __LINE__ << std::endl;
  }
  
  const int numFamilies = 1;
  TD nullTD;
  Kokkos::Array<TD, spaceDim > firstFamilyLeft  {tensorData};
//  Kokkos::Array<TD, spaceDim > secondFamilyLeft {tensorData};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};//, secondFamilyLeft};
  
  VD vectorDataLeft(vectorComponentsLeft);
  
  Kokkos::Array<TD, spaceDim > firstFamilyRight  {tensorData};
  Kokkos::Array<TD, spaceDim > secondFamilyRight {tensorData};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight}; //, secondFamilyRight};
  
  VD vectorDataRight(vectorComponentsRight);
  
  TBV  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TBV transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataLeft);
  
  D constantCellMeasuresData(1.0, Kokkos::Array<int,2>{1,1});
  TD constantCellMeasures(constantCellMeasuresData);
  
  // these assignments imitate a function call with arguments (tbvLeft, cellMeasures, tbvRight)
  const TBV      tbvLeft = transformedUnitVectorDataLeft;
  const TD  cellMeasures = constantCellMeasures;
  const TBV     tbvRight = transformedUnitVectorDataLeft;
  
  using IT = Intrepid2::IntegrationTools<DeviceType>;
  
  auto integralsBaseline  = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  auto integralsIntegrate = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  
  {
    auto data = getMatchingViewWithLabel(identityMatrixView, "Data mat-mat result", 1, 1, 1, 1);
    std::cout << "data.size(): " << data.size() << std::endl;
    std::cout << "Got to line " << __LINE__ << std::endl;
  }
  
  integrate_baseline(integralsBaseline, tbvLeft, cellMeasures, tbvRight);
  
  {
    auto data = getMatchingViewWithLabel(identityMatrixView, "Data mat-mat result", 1, 1, 1, 1);
    std::cout << "data.size(): " << data.size() << std::endl;
    std::cout << "Got to line " << __LINE__ << std::endl;
  }
  
  IT::integrate(integralsIntegrate, tbvLeft, cellMeasures, tbvRight);
}

} // anonymous namespace
