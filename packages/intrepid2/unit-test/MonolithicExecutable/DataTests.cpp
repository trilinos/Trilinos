// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   DataTests.cpp
    \brief  Tests against Intrepid2::Data.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_Data.hpp"
#include "Intrepid2_DataTools.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

namespace
{
  using namespace Intrepid2;

/** \brief Data has facilities for in-place combinations of logical data.  Suppose you have two containers of logical shape (C,P), one of which is constant across cells, the other of which is constant across points.  To combine these (e.g., sum them together entrywise), you want a container that varies in both cells and points.  The test below exercises the facility for allocation of the combined container.
*/
  TEUCHOS_UNIT_TEST( Data, AllocateInPlaceCombinationResult )
  {
    // test allocateInPlaceCombinationResult()
    // Use two Data objects A and B, each with logical shape (5,9,15) -- (C,F,P), say.
    // with A having variation types of GENERAL, MODULAR, and CONSTANT,
    // and B having variation types of CONSTANT, CONSTANT, and GENERAL.
    // Result should have variation types of GENERAL, MODULAR, GENERAL.
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    
    const int rank        = 3;
    const int cellCount   = 5;
    const int fieldCount  = 9;
    const int pointCount  = 15;
    
    const int fieldCountA = 3; // A is modular in field dimension, with variation mod 3.
    auto AView = getView<Scalar,DeviceType>("A", cellCount, fieldCountA);
    auto BView = getView<Scalar,DeviceType>("B", pointCount);
    
    auto ABView = getView<Scalar,DeviceType>("A+B", cellCount, fieldCountA, pointCount);
    
    Kokkos::Array<int,rank> extents {cellCount, fieldCount, pointCount};
    Kokkos::Array<DataVariationType,rank> A_variation {GENERAL, MODULAR, CONSTANT};
    Kokkos::Array<DataVariationType,rank> B_variation {CONSTANT, CONSTANT, GENERAL};
    
    Data<Scalar,DeviceType> A(AView,extents,A_variation);
    Data<Scalar,DeviceType> B(BView,extents,B_variation);
    
    // expected variation for A+B:
    Kokkos::Array<DataVariationType,3> AB_variation {GENERAL, MODULAR, GENERAL};
    // expected Data object for A+B:
    Data<Scalar,DeviceType> AB_expected(ABView,extents,AB_variation);
    
    auto AB_actual = Data<Scalar,DeviceType>::allocateInPlaceCombinationResult(A, B);
    
    TEST_EQUALITY(AB_actual.rank(), AB_expected.rank());
    for (int d=0; d<rank; d++)
    {
      const auto actualVariationType   = AB_actual.getVariationTypes()[d];
      const auto expectedVariationType = AB_expected.getVariationTypes()[d];
      TEST_EQUALITY(actualVariationType, expectedVariationType);
      
      const auto actualVariationModulus   = AB_actual.getVariationModulus(d);
      const auto expectedVariationModulus = AB_expected.getVariationModulus(d);
      TEST_EQUALITY(actualVariationModulus, expectedVariationModulus);
      
      const auto actualExtent   = AB_actual.extent_int(d);
      const auto expectedExtent = AB_expected.extent_int(d);
      TEST_EQUALITY(actualExtent, expectedExtent);
    }
    
    TEST_EQUALITY(AB_actual.getUnderlyingViewRank(), AB_expected.getUnderlyingViewRank());
    const int dataRank = AB_expected.getUnderlyingViewRank();
    if (AB_actual.getUnderlyingViewRank() == dataRank)
    {
      for (int d=0; d<dataRank; d++)
      {
        const auto actualDataExtent   = AB_actual.getDataExtent(d);
        const auto expectedDataExtent = AB_expected.getDataExtent(d);
        TEST_EQUALITY(actualDataExtent, expectedDataExtent);
      }
    }
  }

  TEUCHOS_UNIT_TEST( Data, CombinedDimensionInfo )
  {
    // test free function, combinedDimensionInfo()
    
    DimensionInfo A_dimInfo;
    DimensionInfo B_dimInfo;
    DimensionInfo AB_dimInfo;
    
    A_dimInfo.logicalExtent = 15;
    B_dimInfo.logicalExtent = 15;
    AB_dimInfo.logicalExtent = 15;
    
    A_dimInfo.blockPlusDiagonalLastNonDiagonal = -1;
    B_dimInfo.blockPlusDiagonalLastNonDiagonal = -1;
    AB_dimInfo.blockPlusDiagonalLastNonDiagonal = -1;
    
    A_dimInfo.variationModulus = 15;
    B_dimInfo.variationModulus = 1;
    AB_dimInfo.variationModulus = 15;
    
    A_dimInfo.variationType = GENERAL;
    B_dimInfo.variationType = CONSTANT;
    AB_dimInfo.variationType = GENERAL;
    
    A_dimInfo.dataExtent  =  A_dimInfo.logicalExtent / ( A_dimInfo.logicalExtent /  A_dimInfo.variationModulus);
    B_dimInfo.dataExtent  =  B_dimInfo.logicalExtent / ( B_dimInfo.logicalExtent /  B_dimInfo.variationModulus);
    AB_dimInfo.dataExtent = AB_dimInfo.logicalExtent / (AB_dimInfo.logicalExtent / AB_dimInfo.variationModulus);
    
    // combinedDimensionInfo should commute, so let's test both directions:
    DimensionInfo AB_dimInfoActual_LR = combinedDimensionInfo(A_dimInfo, B_dimInfo);
    DimensionInfo AB_dimInfoActual_RL = combinedDimensionInfo(B_dimInfo, A_dimInfo);
    
    std::vector<DimensionInfo> actualCombinations {AB_dimInfoActual_LR, AB_dimInfoActual_RL};
    
    for (const auto & dimInfoActual : actualCombinations)
    {
      TEST_EQUALITY(dimInfoActual.logicalExtent, AB_dimInfo.logicalExtent);
      TEST_EQUALITY(dimInfoActual.dataExtent, AB_dimInfo.dataExtent);
      TEST_EQUALITY(dimInfoActual.variationType, AB_dimInfo.variationType);
      TEST_EQUALITY(dimInfoActual.variationModulus, AB_dimInfo.variationModulus);
    }
  }

// #pragma mark Data: EmptyDataMarkedAsInvalid
/** \brief When Data containers are constructed without arguments, the isValid() method should return false.  This test confirms that that is the case.
 */
  TEUCHOS_UNIT_TEST( Data, EmptyDataMarkedAsInvalid )
  {
    using DeviceType = DefaultTestDeviceType;
    // check the new default constructor for Data
    Data<double,DeviceType> emptyData;
    TEST_EQUALITY(emptyData.isValid(), false); // empty data container should return false from isValid()
  }

// #pragma mark Data: GetWritableEntry
/** \brief Data has a method intended primarily for internal use, getWritableEntry(), which returns an lvalue reference to a specified location.  This test uses that method to set a particular entry, and checks that the destination data object has the correct value.
*/
  TEUCHOS_UNIT_TEST( Data, GetWritableEntry )
  {
    using Scalar = double;
    using DeviceType = DefaultTestDeviceType;
    
    // no arithmetic, so comparisons should be exact
    double relTol = 0.0;
    double absTol = 0.0;
    
    const int numRows = 3;
    const int numCols = 3;
    const Scalar lastValueToSet = 18;
    auto zeroView = getView<Scalar,DeviceType>("GetWritableEntry view", numRows, numCols);
    Kokkos::deep_copy(zeroView,0);
    Data<double,DeviceType> data(zeroView);
    
    using ExecSpaceType = typename DeviceType::execution_space;
    Kokkos::RangePolicy<ExecSpaceType> policy(0,1); // trivial policy: 1 entry
    Kokkos::parallel_for("set lastVal", policy,
    KOKKOS_LAMBDA (const int &i) {
      auto & lastVal = data.getWritableEntry(numRows-1, numCols-1);
      lastVal = lastValueToSet;
    });
    
    auto expectedView = getView<Scalar,DeviceType>("GetWritableEntry expected view", numRows, numCols);
    auto expectedViewHost = Kokkos::create_mirror_view(expectedView);
    Kokkos::deep_copy(expectedViewHost,0);
    expectedViewHost(numRows-1,numCols-1) = lastValueToSet;
    
    Kokkos::deep_copy(expectedView, expectedViewHost);
    testFloatingEquality2(expectedView, data, relTol, absTol, out, success);
  }

/** \brief Data has facilities for in-place combinations of logical data.  Suppose you have two containers of logical shape (C,P), one of which is constant across cells, the other of which is constant across points.  To combine these (e.g., sum them together entrywise), you want a container that varies in both cells and points.  The test below exercises the facility for allocation of the combined container.
*/

  TEUCHOS_UNIT_TEST( Data, InPlaceSum )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    // Use two Data objects A and B, each with logical shape (5,9,15) -- (C,F,P), say.
    // with A having variation types of GENERAL, MODULAR, and CONSTANT,
    // and B having variation types of CONSTANT, CONSTANT, and GENERAL.
    // Result should have variation types of GENERAL, MODULAR, GENERAL.
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    
    const int rank        = 3;
    const int cellCount   = 5;
    const int fieldCount  = 9;
    const int pointCount  = 15;
    
    auto formula_A = [] (int cellOrdinal, int fieldOrdinal, int pointOrdinal) -> double
    {
      // varies modulus 3 in fieldOrdinal; constant pointwise
      return double(cellOrdinal) + double(fieldOrdinal % 3);
    };
    
    auto formula_B = [] (int cellOrdinal, int fieldOrdinal, int pointOrdinal) -> double
    {
      // constant in cell, field; varies pointwise
      return double(pointOrdinal);
    };
    
    auto sum = [] (const Scalar &a, const Scalar &b) -> Scalar
    {
      return a + b;
    };
    
    const int fieldCountA = 3; // A is modular in field dimension, with variation mod 3.
    auto AView = getView<Scalar,DeviceType>("A", cellCount, fieldCountA);
    auto BView = getView<Scalar,DeviceType>("B", pointCount);
    
    auto ABView = getView<Scalar,DeviceType>("A+B", cellCount, fieldCountA, pointCount);
    
    auto AViewHost  = Kokkos::create_mirror(AView);
    auto BViewHost  = Kokkos::create_mirror(BView);
    auto ABViewHost = Kokkos::create_mirror(ABView);
    for (int cellOrdinal=0; cellOrdinal<cellCount; cellOrdinal++)
    {
      for (int fieldOrdinal=0; fieldOrdinal<fieldCountA; fieldOrdinal++)
      {
        for (int pointOrdinal=0; pointOrdinal<pointCount; pointOrdinal++)
        {
          auto a = formula_A(cellOrdinal,fieldOrdinal,pointOrdinal);
          auto b = formula_B(cellOrdinal,fieldOrdinal,pointOrdinal);
          AViewHost (cellOrdinal,fieldOrdinal) = a;
          BViewHost (pointOrdinal) = b;
          ABViewHost(cellOrdinal,fieldOrdinal,pointOrdinal) = sum(a,b);
        }
      }
    }
    Kokkos::deep_copy( AView,  AViewHost);
    Kokkos::deep_copy( BView,  BViewHost);
    Kokkos::deep_copy(ABView, ABViewHost);
    
    Kokkos::Array<int,rank> extents {cellCount, fieldCount, pointCount};
    Kokkos::Array<DataVariationType,rank> A_variation {GENERAL, MODULAR, CONSTANT};
    Kokkos::Array<DataVariationType,rank> B_variation {CONSTANT, CONSTANT, GENERAL};
    
    Data<Scalar,DeviceType> A(AView,extents,A_variation);
    Data<Scalar,DeviceType> B(BView,extents,B_variation);
    
    // expected variation for A+B:
    Kokkos::Array<DataVariationType,3> AB_variation {GENERAL, MODULAR, GENERAL};
    // expected Data object for A+B:
    Data<Scalar,DeviceType> AB_expected(ABView,extents,AB_variation);
    
    auto AB_actual = Data<Scalar,DeviceType>::allocateInPlaceCombinationResult(A, B);
    
    AB_actual.storeInPlaceSum(A, B);
    
    // test AB_actual equals AB_expected.  (This will iterate over the logical extents.)
    testFloatingEquality3(AB_actual, AB_expected, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( Data, InPlaceProduct_Matrix )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    // CASE 1
    // Use two Data objects A and B, each with logical shape (10,3,2,2) -- (C,P,D,D), say.
    // with A having variation types of CONSTANT, CONSTANT, and BLOCK_DIAGONAL (with the final two dimensions being purely diagonal)
    // and B having variation types of CONSTANT for all entries
    // Result should have variation types of CONSTANT, CONSTANT, and BLOCK_DIAGONAL, with the diagonal matrix entries weighted by the value from B.
    // (This test is modeled on Jacobians for a uniform geometry being weighted by inverse determinants, as happens when computing inputs to getHDIVtransformVALUE().
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    
    const int rank        = 4;
    const int cellCount   = 10;
    const int pointCount  = 3;
    const int spaceDim    = 2;
    
    const double bValue = M_PI;
    
    auto AView = getView<Scalar,DeviceType>("A", spaceDim); // diagonal: one entry per dimension
    auto ABView = getView<Scalar,DeviceType>("A .* B", spaceDim); // should be same shape as A
    
    auto AViewHost  = Kokkos::create_mirror(AView);
    auto ABViewHost = Kokkos::create_mirror(ABView);
    for (int d=0; d<spaceDim; d++)
    {
      AViewHost(d)  = d + 1.;
      ABViewHost(d) = AViewHost(d) * bValue;
    }
    Kokkos::deep_copy( AView,  AViewHost);
    Kokkos::deep_copy(ABView, ABViewHost);
    
    Kokkos::Array<int,rank> extents {cellCount, pointCount, spaceDim, spaceDim};
    Kokkos::Array<DataVariationType,rank> A_variation {CONSTANT, CONSTANT, BLOCK_PLUS_DIAGONAL, BLOCK_PLUS_DIAGONAL};
    
    Data<Scalar,DeviceType> A(AView,extents,A_variation);
    Data<Scalar,DeviceType> B(bValue,extents);
    
    // expected variation for A+B:
    Kokkos::Array<DataVariationType,4> AB_variation  = A_variation;
    // expected Data object for A+B:
    Data<Scalar,DeviceType> AB_expected(ABView,extents,AB_variation);
    
    auto AB_actual = Data<Scalar,DeviceType>::allocateInPlaceCombinationResult(A, B);
    
    AB_actual.storeInPlaceProduct(A, B);
    
    // test AB_actual equals AB_expected.  (This will iterate over the logical extents.)
    testFloatingEquality4(AB_actual, AB_expected, relTol, absTol, out, success);
    
    // CASE 2
    // now do the same, but with A having a full matrix, rather than just a diagonal
    AView = getView<Scalar,DeviceType>("A", spaceDim, spaceDim);
    ABView = getView<Scalar,DeviceType>("A .* B", spaceDim, spaceDim); // should be same shape as A
    
    AViewHost  = Kokkos::create_mirror(AView);
    ABViewHost = Kokkos::create_mirror(ABView);
    for (int d0=0; d0<spaceDim; d0++)
    {
      for (int d1=0; d1<spaceDim; d1++)
      {
        AViewHost(d0,d1)  = d0 + d0 * d1 + 1.; // some arbitrary data
        ABViewHost(d0,d1) = AViewHost(d0,d1) * bValue;
      }
    }
    Kokkos::deep_copy( AView,  AViewHost);
    Kokkos::deep_copy(ABView, ABViewHost);
   
    A_variation[2] = GENERAL; // replaces BLOCK_PLUS_DIAGONAL
    A_variation[3] = GENERAL; // replaces BLOCK_PLUS_DIAGONAL
    
    A = Data<Scalar,DeviceType>(AView,extents,A_variation);
    
    // expected variation for A*B:
    AB_variation  = A_variation;
    // expected Data object for A*B:
    AB_expected = Data<Scalar,DeviceType>(ABView,extents,AB_variation);
    
    AB_actual = Data<Scalar,DeviceType>::allocateInPlaceCombinationResult(A, B);
    
    AB_actual.storeInPlaceProduct(A, B);
    
    // test AB_actual equals AB_expected.  (This will iterate over the logical extents.)
    testFloatingEquality4(AB_actual, AB_expected, relTol, absTol, out, success);
    
    // CASE 3
    // now try a case where B has logical shape (C,P), and we use DataTools to produce the result.  This is what gets invoked by TransformedBasisValues.
    Kokkos::Array<int,2> b_extents {cellCount, pointCount};
    B = Data<Scalar,DeviceType> (bValue,b_extents);
    
    AB_actual = DataTools::multiplyByCPWeights(A, B);
    // test AB_actual equals AB_expected.  (This will iterate over the logical extents.)
    testFloatingEquality4(AB_actual, AB_expected, relTol, absTol, out, success);
    
    // CASE 4
    // and now one with an actual (C,P) view for B (instead of CONSTANT data in (C,P))
    auto BView = getView<Scalar,DeviceType>("B", cellCount, pointCount);
    Kokkos::deep_copy(BView, bValue);
    Kokkos::Array<DataVariationType,2> B_variation {GENERAL, GENERAL};
    B = Data<Scalar,DeviceType> (BView,b_extents,B_variation);
    AB_actual = DataTools::multiplyByCPWeights(A, B);
    // test AB_actual equals AB_expected.  (This will iterate over the logical extents.)
    testFloatingEquality4(AB_actual, AB_expected, relTol, absTol, out, success);
    
    // CASE 5
    // and now one with an actual (C,P,D,D) view for A (instead of CONSTANT data in (C,P))
    AView = getView<Scalar,DeviceType>("A", cellCount, pointCount, spaceDim, spaceDim);
    ABView = getView<Scalar,DeviceType>("A .* B", cellCount, pointCount, spaceDim, spaceDim); // same shape as A
    
    AViewHost  = Kokkos::create_mirror(AView);
    ABViewHost = Kokkos::create_mirror(ABView);
    for (int cellOrdinal=0; cellOrdinal<cellCount; cellOrdinal++)
    {
      for (int pointOrdinal=0; pointOrdinal<pointCount; pointOrdinal++)
      {
        for (int d0=0; d0<spaceDim; d0++)
        {
          for (int d1=0; d1<spaceDim; d1++)
          {
            AViewHost(cellOrdinal,pointOrdinal,d0,d1)  = d0 + d0 * d1 + 1. + cellOrdinal * (pointOrdinal+1); // some arbitrary data
            ABViewHost(cellOrdinal,pointOrdinal,d0,d1) = AViewHost(cellOrdinal,pointOrdinal,d0,d1) * bValue;
          }
        }
      }
    }
    Kokkos::deep_copy( AView,  AViewHost);
    Kokkos::deep_copy(ABView, ABViewHost);
   
    A_variation = Kokkos::Array<DataVariationType,rank>{GENERAL,GENERAL,GENERAL,GENERAL};
    
    A = Data<Scalar,DeviceType>(AView,extents,A_variation);
    
    // expected variation for A*B:
    AB_variation  = A_variation;
    // expected Data object for A*B:
    AB_expected = Data<Scalar,DeviceType>(ABView,extents,AB_variation);
    
    AB_actual = DataTools::multiplyByCPWeights(A, B);
    
    // test AB_actual equals AB_expected.  (This will iterate over the logical extents.)
    testFloatingEquality4(AB_actual, AB_expected, relTol, absTol, out, success);
  }

// #pragma mark Data: MatVec
/** \brief Data provides matrix-vector multiplication support.  This method checks correctness of the computed mat-vec for a particular case involving a 2x2 matrix and a 2x1 vector.
*/
  TEUCHOS_UNIT_TEST( Data, MatVec )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    const int spaceDim = 2;
    auto matrixView = getView<Scalar,DeviceType>("full matrix", spaceDim, spaceDim);
    auto matrixViewHost = Kokkos::create_mirror(matrixView);
    matrixViewHost(0,0) =  1.0;  matrixViewHost(0,1) =  2.0;
    matrixViewHost(1,0) = -1.0;  matrixViewHost(1,1) =  3.0;
    Kokkos::deep_copy(matrixView, matrixViewHost);
    
    auto vecView = getView<Scalar,DeviceType>("vector", spaceDim);
    auto vecViewHost = Kokkos::create_mirror(vecView);
    
    vecViewHost(0) = 1.0;
    vecViewHost(1) = 2.0;
    Kokkos::deep_copy(vecView, vecViewHost);
    
    auto expectedResultView = getView<Scalar,DeviceType>("result vector", spaceDim);
    auto expectedResultViewHost = Kokkos::create_mirror(expectedResultView);
    
    expectedResultViewHost(0) = matrixViewHost(0,0) * vecViewHost(0) + matrixViewHost(0,1) * vecViewHost(1);
    expectedResultViewHost(1) = matrixViewHost(1,0) * vecViewHost(0) + matrixViewHost(1,1) * vecViewHost(1);
    
    Kokkos::deep_copy(expectedResultView, expectedResultViewHost);
    
    Data<Scalar,DeviceType> matData(matrixView);
    Data<Scalar,DeviceType> vecData(vecView);
    auto actualResultData = Data<Scalar,DeviceType>::allocateMatVecResult(matData, vecData);
    actualResultData.storeMatVec(matData, vecData);
    
    testFloatingEquality1(expectedResultView, actualResultData.getUnderlyingView1(), relTol, absTol, out, success);
  }

// #pragma mark Data: MatVec_CPDD_transpose
/** \brief Data provides matrix-vector multiplication support.  This method checks correctness of the computed mat-vec for a particular case involving a 2x2 matrix stored in a rank 4 container such as will arise for Jacobians; here, the matrix is stored in a rank-4 View, and the Data object is a thin wrapper around it (i.e., the DataVariationType for each of the dimensions in data is GENERAL).
*/
  TEUCHOS_UNIT_TEST( Data, MatVec_CPDD_transpose )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    const int spaceDim = 2;
    const int cellCount = 1;
    const int pointCount = 1;
    auto matrixView = getView<Scalar,DeviceType>("matrix", cellCount, pointCount, spaceDim, spaceDim);
    auto matrixViewHost = Kokkos::create_mirror(matrixView);
    matrixViewHost(0,0,0,0) =  1.0;  matrixViewHost(0,0,0,1) =  2.0;
    matrixViewHost(0,0,1,0) = -1.0;  matrixViewHost(0,0,1,1) =  3.0;
    Kokkos::deep_copy(matrixView, matrixViewHost);
    
    auto vectorView = getView<Scalar,DeviceType>("vector", cellCount, pointCount, spaceDim);
    auto vectorViewHost = Kokkos::create_mirror(vectorView);
    vectorViewHost(0,0,0) =  1.0;
    vectorViewHost(0,0,1) = -1.0;
    Kokkos::deep_copy(vectorView, vectorViewHost);
    
    auto expectedResultView = getView<Scalar,DeviceType>("result vector", cellCount, pointCount, spaceDim);
    auto expectedResultViewHost = Kokkos::create_mirror(expectedResultView);
    
    std::vector<bool> transposeChoices {false, true};
    
    for (auto transpose : transposeChoices)
    {
      const int cellOrdinal  = 0;
      const int pointOrdinal = 0;
      for (int i=0; i<spaceDim; i++)
      {
        Scalar result_i =  0;
        for (int j=0; j<spaceDim; j++)
        {
          const auto & mat_ij = transpose ? matrixViewHost(cellOrdinal,pointOrdinal,j,i) : matrixViewHost(cellOrdinal,pointOrdinal,i,j);
          result_i += mat_ij * vectorViewHost(cellOrdinal,pointOrdinal,j);
        }
        expectedResultViewHost(cellOrdinal,pointOrdinal,i) = result_i;
      }
      Kokkos::deep_copy(expectedResultView, expectedResultViewHost);
      
      Data<Scalar,DeviceType> A_data(matrixView);
      Data<Scalar,DeviceType> x_data(vectorView);
      auto actualResultData = Data<Scalar,DeviceType>::allocateMatVecResult(A_data, x_data, transpose);
      
      TEST_EQUALITY(         3, actualResultData.rank());
      TEST_EQUALITY( cellCount, actualResultData.extent_int(0));
      TEST_EQUALITY(pointCount, actualResultData.extent_int(1));
      TEST_EQUALITY(  spaceDim, actualResultData.extent_int(2));
      
      actualResultData.storeMatVec(A_data, x_data, transpose);
      
      testFloatingEquality3(expectedResultView, actualResultData, relTol, absTol, out, success);
      
      printView(actualResultData.getUnderlyingView3(), out);
    }
    
    // now, check that u' A v = v' A' u for arbitrary vectors u,v
    
    // set up a second vector (v)
    auto vector2View = getView<Scalar,DeviceType>("vector2", cellCount, pointCount, spaceDim);
    auto vector2ViewHost = Kokkos::create_mirror(vector2View);
    vector2ViewHost(0,0,0) =  3.0;
    vector2ViewHost(0,0,1) =  2.0;
    Kokkos::deep_copy(vector2View, vector2ViewHost);
    
    Data<Scalar,DeviceType> u_data(vectorView);
    Data<Scalar,DeviceType> A_data(matrixView);
    Data<Scalar,DeviceType> v_data(vector2View);
    
    auto AvResultData = Data<Scalar,DeviceType>::allocateMatVecResult(A_data, v_data, false);
    AvResultData.storeMatVec(A_data, v_data, false);
    
    auto upAvResultData = Data<Scalar,DeviceType>::allocateDotProductResult(u_data, AvResultData);
    upAvResultData.storeDotProduct(u_data, AvResultData);
      
    auto ApuResultData = Data<Scalar,DeviceType>::allocateMatVecResult(A_data, u_data, true);
    ApuResultData.storeMatVec(A_data, u_data, true);
    
    auto vpAuResultData = Data<Scalar,DeviceType>::allocateDotProductResult(v_data, ApuResultData);
    vpAuResultData.storeDotProduct(v_data, ApuResultData);
    
    testFloatingEquality2(upAvResultData, vpAuResultData, relTol, absTol, out, success);
    printView(upAvResultData.getUnderlyingView2(), out);
    printView(vpAuResultData.getUnderlyingView2(), out);
  }

// #pragma mark Data: MatMat
/** \brief Data provides matrix-matrix multiplication support.  This method checks correctness of the computed mat-mat for a particular case involving two 2x2 matrices; here, the matrices are each stored in a rank-3 View, and the Data object is a thin wrapper around these (i.e., the DataVariationType for each of the dimensions in data is GENERAL).
*/
  TEUCHOS_UNIT_TEST( Data, MatMat )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    const int spaceDim = 2;
    const int cellCount = 1;
    auto leftMatrixView = getView<Scalar,DeviceType>("left matrix", cellCount, spaceDim, spaceDim);
    auto leftMatrixViewHost = Kokkos::create_mirror(leftMatrixView);
    leftMatrixViewHost(0,0,0) =  1.0;  leftMatrixViewHost(0,0,1) =  2.0;
    leftMatrixViewHost(0,1,0) = -1.0;  leftMatrixViewHost(0,1,1) =  3.0;
    Kokkos::deep_copy(leftMatrixView, leftMatrixViewHost);
    
    auto rightMatrixView = getView<Scalar,DeviceType>("right matrix", cellCount, spaceDim, spaceDim);
    auto rightMatrixViewHost = Kokkos::create_mirror(rightMatrixView);
    rightMatrixViewHost(0,0,0) =  1.0;  rightMatrixViewHost(0,0,1) =  2.0;
    rightMatrixViewHost(0,1,0) = -1.0;  rightMatrixViewHost(0,1,1) =  3.0;
    Kokkos::deep_copy(rightMatrixView, rightMatrixViewHost);
    
    auto expectedResultView = getView<Scalar,DeviceType>("result matrix", cellCount, spaceDim, spaceDim);
    auto expectedResultViewHost = Kokkos::create_mirror(expectedResultView);
    
    const int cellOrdinal = 0;
    for (int i=0; i<spaceDim; i++)
    {
      for (int j=0; j<spaceDim; j++)
      {
        Scalar result =  0;
        for (int k=0; k<spaceDim; k++)
        {
          const auto & left  =  leftMatrixViewHost(cellOrdinal,i,k);
          const auto & right = rightMatrixViewHost(cellOrdinal,k,j);
          result += left * right;
        }
        expectedResultViewHost(cellOrdinal,i,j) = result;
      }
    }
    Kokkos::deep_copy(expectedResultView, expectedResultViewHost);
    
    // TODO: add tests for other transpose possibilities
    const bool transposeA = false;
    const bool transposeB = false;
    
    Data<Scalar,DeviceType> A_data(leftMatrixView);
    Data<Scalar,DeviceType> B_data(rightMatrixView);
    auto actualResultData = Data<Scalar,DeviceType>::allocateMatMatResult(transposeA, A_data, transposeB, B_data);
    
    TEST_EQUALITY(       3,  actualResultData.rank());
    TEST_EQUALITY(cellCount, actualResultData.extent_int(0));
    TEST_EQUALITY(spaceDim,  actualResultData.extent_int(1));
    TEST_EQUALITY(spaceDim,  actualResultData.extent_int(2));
    
    actualResultData.storeMatMat(transposeA, A_data, transposeB, B_data);
    
    testFloatingEquality3(expectedResultView, actualResultData, relTol, absTol, out, success);
    
    printView(actualResultData.getUnderlyingView3(), out);
  }

/** \brief Data provides matrix-matrix multiplication support.  This method checks correctness of the computed mat-mat for a case arising from taking the outer product of two vectors.
*/
  TEUCHOS_UNIT_TEST( Data, MatMatOuterProduct )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    const int spaceDim = 2;
    const int cellCount = 1;
    const int pointCount = 1;
    auto leftVectorView = getView<Scalar,DeviceType>("left vector", cellCount, pointCount, spaceDim);
    auto leftVectorViewHost = Kokkos::create_mirror(leftVectorView);
    leftVectorViewHost(0,0,0) = 1.0;
    leftVectorViewHost(0,0,1) = 0.5;
    Kokkos::deep_copy(leftVectorView, leftVectorViewHost);
    
    Data<Scalar,DeviceType> leftVector(leftVectorView);
    
    auto rightVectorView = getView<Scalar,DeviceType>("right vector", cellCount, pointCount, spaceDim);
    auto rightVectorViewHost = Kokkos::create_mirror(rightVectorView);
    rightVectorViewHost(0,0,0) = 0.5;
    rightVectorViewHost(0,0,1) = 1.0;
    Kokkos::deep_copy(rightVectorView, rightVectorViewHost);
    Data<Scalar,DeviceType> rightVector(rightVectorView);
    
    // re-cast leftVector as a rank 4 (C,P,D,1) object -- a D x 1 matrix at each (C,P).
    const int newRank   = 4;
    auto extents        = leftVector.getExtents();
    auto variationTypes = leftVector.getVariationTypes();
    auto leftMatrix = leftVector.shallowCopy(newRank, extents, variationTypes);
    
    // re-cast rightVector as a rank 4 (C,P,1,D) object -- a 1 x D matrix at each (C,P)
    extents           = rightVector.getExtents();
    extents[3]        = extents[2];
    extents[2]        = 1;
    variationTypes    = rightVector.getVariationTypes();
    variationTypes[3] = variationTypes[2];
    variationTypes[2] = CONSTANT;
    auto rightMatrix  = rightVector.shallowCopy(newRank, extents, variationTypes);
    
    auto expectedResultView = getView<Scalar,DeviceType>("result matrix", cellCount, pointCount, spaceDim, spaceDim);
    auto expectedResultViewHost = Kokkos::create_mirror(expectedResultView);
    
    const int cellOrdinal = 0;
    for (int i=0; i<spaceDim; i++)
    {
      for (int j=0; j<spaceDim; j++)
      {
        const auto & left  =  leftVectorViewHost(cellOrdinal,0,i);
        const auto & right = rightVectorViewHost(cellOrdinal,0,j);
        Scalar result = left * right;
        expectedResultViewHost(cellOrdinal,0,i,j) = result;
      }
    }
    Kokkos::deep_copy(expectedResultView, expectedResultViewHost);
    
    const bool transposeA = false;
    const bool transposeB = false;
    
    auto actualResultData = Data<Scalar,DeviceType>::allocateMatMatResult(transposeA, leftMatrix, transposeB, rightMatrix);
    
    TEST_EQUALITY(         4, actualResultData.rank());
    TEST_EQUALITY( cellCount, actualResultData.extent_int(0));
    TEST_EQUALITY(pointCount, actualResultData.extent_int(1));
    TEST_EQUALITY(  spaceDim, actualResultData.extent_int(2));
    TEST_EQUALITY(  spaceDim, actualResultData.extent_int(3));
    
    actualResultData.storeMatMat(transposeA, leftMatrix, transposeB, rightMatrix);
    
    testFloatingEquality4(expectedResultView, actualResultData, relTol, absTol, out, success);
    
    printView(actualResultData.getUnderlyingView(), out);
  }

// #pragma mark Data: MatMatExplicitIdentity_PDD
/** \brief Data provides matrix-matrix multiplication support.  This method checks correctness of the computed mat-mat for several cases involving 3x3 identity matrices.  Here, the logical dimensions (C,P,D,D) differ from the stored dimensions of (P,D,D).  We test each possible transpose combination.
*/
TEUCHOS_UNIT_TEST( Data, MatMatExplicitIdentity_PDD ) // (P,D,D) underlying; notionally (C,P,D,D)
{
  // tests that identity * identity = identity
  using DeviceType = DefaultTestDeviceType;
  
  double relTol = 1e-13;
  double absTol = 1e-13;
  
  using Scalar = double;
  const int spaceDim = 3;
  const int cellCount = 1;
  const int numPoints = 4;
  
  auto identityMatrixView = getFixedRankView<Scalar>("identity matrix", numPoints, spaceDim, spaceDim);
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
  
  Kokkos::Array<int,4> transformationExtents {cellCount, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, GENERAL, GENERAL, GENERAL};
  
  Data<Scalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  using std::vector;
  using std::pair;
  vector<pair<bool,bool>> transposeABChoices {
    pair<bool,bool>{false,false},
    pair<bool,bool>{false,true},
    pair<bool,bool>{true,false},
    pair<bool,bool>{true,true}
  };
    
  for (const auto &transposeAB : transposeABChoices)
  {
    const bool transposeA = transposeAB.first;
    const bool transposeB = transposeAB.second;
    
    using std::string;
    const string transposeAString = transposeA ? "true" : "false";
    const string transposeBString = transposeB ? "true" : "false";
  
    out << "*** Testing transposeA = " << transposeAString << ", transposeB = " << transposeBString << " ***\n";
    auto actualResultData = Data<Scalar,DeviceType>::allocateMatMatResult(transposeA, explicitIdentityMatrix, transposeB, explicitIdentityMatrix);
    
    TEST_EQUALITY(       4,  actualResultData.rank());
    TEST_EQUALITY(cellCount, actualResultData.extent_int(0));
    TEST_EQUALITY(numPoints, actualResultData.extent_int(1));
    TEST_EQUALITY(spaceDim,  actualResultData.extent_int(2));
    TEST_EQUALITY(spaceDim,  actualResultData.extent_int(3));
    
    actualResultData.storeMatMat(transposeA, explicitIdentityMatrix, transposeB, explicitIdentityMatrix);

    testFloatingEquality4(explicitIdentityMatrix, actualResultData, relTol, absTol, out, success, "expected", "actual");
    
    printView(actualResultData.getUnderlyingView3(), out);
  }
}

  // #pragma mark Data: MatMatBlockPlusDiagonal
/** \brief Data provides matrix-matrix multiplication support.  This method checks correctness of the computed mat-mat for a case involving one 3x3 matrix that has a 2x2 upper left block, and diagonal entry in the (3,3) position, and one 3x3 matrix that is entirely diagonal.  Here, the logical dimensions (C,D,D) match the stored dimensions.
*/
  TEUCHOS_UNIT_TEST( Data, MatMatBlockPlusDiagonal )
  {
    using DeviceType = DefaultTestDeviceType;
    
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    using Scalar = double;
    const int spaceDim = 3;
    const int leftLastNonDiagonal = 1;
    const int cellCount = 1;
    auto leftMatrixView = getView<Scalar,DeviceType>("left matrix", cellCount, (leftLastNonDiagonal+1) * (leftLastNonDiagonal+1) + (spaceDim - leftLastNonDiagonal - 1));
    auto leftMatrixViewHost = Kokkos::create_mirror(leftMatrixView);
    int entryIndex = 0;
    // Block:
    leftMatrixViewHost(0,entryIndex++) =  1.0;  leftMatrixViewHost(0,entryIndex++) =  2.0;
    leftMatrixViewHost(0,entryIndex++) = -1.0;  leftMatrixViewHost(0,entryIndex++) =  3.0;
    // Diagonal:
    leftMatrixViewHost(0,entryIndex++) = 3.0;
    Kokkos::deep_copy(leftMatrixView, leftMatrixViewHost);
    
    const int rightLastNonDiagonal = 0;
    auto rightMatrixView = getView<Scalar,DeviceType>("right matrix", cellCount, (rightLastNonDiagonal+1) * (rightLastNonDiagonal+1) + (spaceDim - rightLastNonDiagonal - 1));
    auto rightMatrixViewHost = Kokkos::create_mirror(rightMatrixView);
    entryIndex = 0;
    // Diagonal:
    rightMatrixViewHost(0,entryIndex++) =  1.0;
                                                rightMatrixViewHost(0,entryIndex++) =  3.0;
                                                                                            rightMatrixViewHost(0,entryIndex++) =  2.0;
    Kokkos::deep_copy(rightMatrixView, rightMatrixViewHost);
    
    // TODO: add tests for other transpose possibilities
    const bool transposeA = false;
    const bool transposeB = false;
    
    const int rank = 3;
    Kokkos::Array<int,7> extents {cellCount,spaceDim,spaceDim,0,0,0,0};
    Kokkos::Array<DataVariationType,7> variationTypes {GENERAL,BLOCK_PLUS_DIAGONAL,BLOCK_PLUS_DIAGONAL,CONSTANT,CONSTANT,CONSTANT,CONSTANT};
    
    Data<Scalar,DeviceType> A( leftMatrixView, rank, extents, variationTypes, leftLastNonDiagonal);
    Data<Scalar,DeviceType> B(rightMatrixView, rank, extents, variationTypes, rightLastNonDiagonal);
    
    auto expectedResultView = getView<Scalar,DeviceType>("result matrix", cellCount, spaceDim, spaceDim);
    
    const int cellOrdinal = 0;
    auto policy = Kokkos::MDRangePolicy<typename DeviceType::execution_space,Kokkos::Rank<2>>({0,0},{spaceDim,spaceDim});
    Kokkos::parallel_for("compute first-order simplex Jacobians", policy,
    KOKKOS_LAMBDA(const int &i, const int &j)
    {
      Scalar result =  0;
      for (int k=0; k<spaceDim; k++)
      {
        Scalar left  = transposeA ? A(cellOrdinal,k,i) : A(cellOrdinal,i,k);
        Scalar right = transposeB ? B(cellOrdinal,j,k) : B(cellOrdinal,k,j);
        result += left * right;
      }
      expectedResultView(cellOrdinal, i, j) = result;
    });
        
    auto actualResultData = Data<Scalar,DeviceType>::allocateMatMatResult(transposeA, A, transposeB, B);
    
    TEST_EQUALITY(       3,  actualResultData.rank());
    TEST_EQUALITY(cellCount, actualResultData.extent_int(0));
    TEST_EQUALITY(spaceDim,  actualResultData.extent_int(1));
    TEST_EQUALITY(spaceDim,  actualResultData.extent_int(2));
    
    const int resultLastNonDiagonal = std::max(leftLastNonDiagonal,rightLastNonDiagonal);
    TEST_EQUALITY(resultLastNonDiagonal, actualResultData.blockPlusDiagonalLastNonDiagonal());
    
    actualResultData.storeMatMat(transposeA, A, transposeB, B);
    
    testFloatingEquality3(expectedResultView, actualResultData, relTol, absTol, out, success);
    
    printView(actualResultData.getUnderlyingView2(), out);
  }

// #pragma mark Data: VecDotProduct
/** \brief Data provides vector dot product multiplication support.  This method checks correctness of the computed dot product for a particular case involving 2x1 vectors.
*/
  TEUCHOS_UNIT_TEST( Data, VecDotProduct )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    const int numCells = 1;
    const int spaceDim = 2;
    
    auto vec1View = getView<Scalar,DeviceType>("vector", numCells, spaceDim);
    auto vec1ViewHost = Kokkos::create_mirror(vec1View);
    
    vec1ViewHost(0,0) = 1.0;
    vec1ViewHost(0,1) = 2.0;
    Kokkos::deep_copy(vec1View, vec1ViewHost);
    
    auto vec2View = getView<Scalar,DeviceType>("vector", numCells, spaceDim);
    auto vec2ViewHost = Kokkos::create_mirror(vec1View);
    
    vec2ViewHost(0,0) = 3.0;
    vec2ViewHost(0,1) = 2.0;
    Kokkos::deep_copy(vec2View, vec2ViewHost);
    
    auto expectedResultView = getView<Scalar,DeviceType>("result",numCells);
    auto expectedResultViewHost = Kokkos::create_mirror(expectedResultView);
    
    expectedResultViewHost(0) = vec1ViewHost(0,0) * vec2ViewHost(0,0) + vec1ViewHost(0,1) * vec2ViewHost(0,1);
    
    Kokkos::deep_copy(expectedResultView, expectedResultViewHost);
    
    Data<Scalar,DeviceType> vec1Data(vec1View);
    Data<Scalar,DeviceType> vec2Data(vec2View);
    auto actualResultData = Data<Scalar,DeviceType>::allocateDotProductResult(vec1Data, vec2Data);
    actualResultData.storeDotProduct(vec1Data, vec2Data);
    
    testFloatingEquality1(expectedResultView, actualResultData.getUnderlyingView1(), relTol, absTol, out, success);
  }
  
  // test statically that Data supports all 7 rank operators
  static_assert(supports_rank<Data<double,DefaultTestDeviceType>,1>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double,DefaultTestDeviceType>,2>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double,DefaultTestDeviceType>,3>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double,DefaultTestDeviceType>,4>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double,DefaultTestDeviceType>,5>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double,DefaultTestDeviceType>,6>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double,DefaultTestDeviceType>,7>::value, "Data is expected to support up to rank 7");
} // anonymous namespace
