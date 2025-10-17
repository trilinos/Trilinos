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

TEUCHOS_UNIT_TEST( PR14546, Distill14546SegFault )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  using D  = Data<DataScalar,DeviceType>;
  using TD = TensorData<DataScalar,DeviceType>;
  using VD = VectorData<DataScalar,DeviceType>;
  using TBV = TransformedBasisValues<DataScalar,DeviceType>;
  
  const int spaceDim = 2;
  
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
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", 1, 1);
  Kokkos::deep_copy(identityMatrixView, 1.0);
  
  Kokkos::Array<int,4> transformExtents {1, 1, 1, 1};
  Kokkos::Array<DataVariationType,4> transformationVariationType {GENERAL, GENERAL, GENERAL, GENERAL};
  
  D explicitIdentityMatrix(identityMatrixView, transformExtents, transformationVariationType);
  
  const int numFamilies = 2;
  TD nullTD;
  Kokkos::Array<TD, spaceDim > firstFamilyLeft  {tensorData,nullTD};
  Kokkos::Array<TD, spaceDim > secondFamilyLeft {nullTD,tensorData};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft, secondFamilyLeft};
  
  VD vectorDataLeft(vectorComponentsLeft);
  
  Kokkos::Array<TD, spaceDim > firstFamilyRight  {tensorData,nullTD};
  Kokkos::Array<TD, spaceDim > secondFamilyRight {nullTD,tensorData};
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
  
  auto integralsBaseline  = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  auto integralsIntegrate = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  
  integrate_baseline(integralsBaseline, tbvLeft, cellMeasures, tbvRight);
  IT::integrate(integralsIntegrate, tbvLeft, cellMeasures, tbvRight);
}

} // anonymous namespace
