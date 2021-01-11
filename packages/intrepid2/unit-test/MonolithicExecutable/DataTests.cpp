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

/** \file   DataTests.cpp
    \brief  Tests against Intrepid2::Data.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_Data.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

namespace
{
  using namespace Intrepid2;

// #pragma mark Data: EmptyDataMarkedAsInvalid
/** \brief When Data containers are constructed without arguments, the isValid() method should return false.  This test confirms that that is the case.
 */
  TEUCHOS_UNIT_TEST( Data, EmptyDataMarkedAsInvalid )
  {
    // check the new default constructor for Data
    Data<double> emptyData;
    TEST_EQUALITY(emptyData.isValid(), false); // empty data container should return false from isValid()
  }

// #pragma mark Data: GetWritableEntry
/** \brief Data has a method intended primarily for internal use, getWritableEntry(), which returns an lvalue reference to a specified location.  This test uses that method to set a particular entry, and checks that the destination data object has the correct value.
*/
  TEUCHOS_UNIT_TEST( Data, GetWritableEntry )
  {
    using Scalar = double;
    using ExecSpaceType = Kokkos::DefaultExecutionSpace;
    
    // no arithmetic, so comparisons should be exact
    double relTol = 0.0;
    double absTol = 0.0;
    
    const int numRows = 3;
    const int numCols = 3;
    const Scalar lastValueToSet = 18;
    auto zeroView = getView<Scalar>("GetWritableEntry view", numRows, numCols);
    Kokkos::deep_copy(zeroView,0);
    Data<double> data(zeroView);
    
    Kokkos::RangePolicy<ExecSpaceType> policy(0,1); // trivial policy: 1 entry
    Kokkos::parallel_for("set lastVal", policy,
    KOKKOS_LAMBDA (const int &i) {
      auto & lastVal = data.getWritableEntry(numRows-1, numCols-1, 0, 0, 0, 0, 0);
      lastVal = lastValueToSet;
    });
    
    auto expectedView = getView<Scalar>("GetWritableEntry expected view", numRows, numCols);
    auto expectedViewHost = Kokkos::create_mirror_view(expectedView);
    Kokkos::deep_copy(expectedViewHost,0);
    expectedViewHost(numRows-1,numCols-1) = lastValueToSet;
    
    Kokkos::deep_copy(expectedView, expectedViewHost);
    testFloatingEquality2(expectedView, data, relTol, absTol, out, success);
  }

// #pragma mark Data: MatVec
/** \brief Data provides matrix-vector multiplication support.  This method checks correctness of the computed mat-vec for a particular case involving a 2x2 matrix and a 2x1 vector.
*/
  TEUCHOS_UNIT_TEST( Data, MatVec )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    using Scalar = double;
    const int spaceDim = 2;
    auto matrixView = getView<Scalar>("full matrix", spaceDim, spaceDim);
    auto matrixViewHost = Kokkos::create_mirror(matrixView);
    matrixViewHost(0,0) =  1.0;  matrixViewHost(0,1) =  2.0;
    matrixViewHost(1,0) = -1.0;  matrixViewHost(1,1) =  3.0;
    Kokkos::deep_copy(matrixView, matrixViewHost);
    
    auto vecView = getView<Scalar>("vector", spaceDim);
    auto vecViewHost = Kokkos::create_mirror(vecView);
    
    vecViewHost(0) = 1.0;
    vecViewHost(1) = 2.0;
    Kokkos::deep_copy(vecView, vecViewHost);
    
    auto expectedResultView = getView<Scalar>("result vector", spaceDim);
    auto expectedResultViewHost = Kokkos::create_mirror(expectedResultView);
    
    expectedResultViewHost(0) = matrixViewHost(0,0) * vecViewHost(0) + matrixViewHost(0,1) * vecViewHost(1);
    expectedResultViewHost(1) = matrixViewHost(1,0) * vecViewHost(0) + matrixViewHost(1,1) * vecViewHost(1);
    
    Kokkos::deep_copy(expectedResultView, expectedResultViewHost);
    
    Data<Scalar> matData(matrixView);
    Data<Scalar> vecData(vecView);
    auto actualResultData = Data<Scalar>::allocateMatVecResult(matData, vecData);
    actualResultData.storeMatVec(matData, vecData);
    
    testFloatingEquality1(expectedResultView, actualResultData.getUnderlyingView1(), relTol, absTol, out, success);
  }

// #pragma mark Data: MatMat
/** \brief Data provides matrix-matrix multiplication support.  This method checks correctness of the computed mat-mat for a particular case involving two 2x2 matrices; here, the matrices are each stored in a rank-3 View, and the Data object is a thin wrapper around these (i.e., the DataVariationType for each of the dimensions in data is GENERAL).
*/
  TEUCHOS_UNIT_TEST( Data, MatMat )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    using Scalar = double;
    const int spaceDim = 2;
    const int cellCount = 1;
    auto leftMatrixView = getView<Scalar>("left matrix", cellCount, spaceDim, spaceDim);
    auto leftMatrixViewHost = Kokkos::create_mirror(leftMatrixView);
    leftMatrixViewHost(0,0,0) =  1.0;  leftMatrixViewHost(0,0,1) =  2.0;
    leftMatrixViewHost(0,1,0) = -1.0;  leftMatrixViewHost(0,1,1) =  3.0;
    Kokkos::deep_copy(leftMatrixView, leftMatrixViewHost);
    
    auto rightMatrixView = getView<Scalar>("right matrix", cellCount, spaceDim, spaceDim);
    auto rightMatrixViewHost = Kokkos::create_mirror(rightMatrixView);
    rightMatrixViewHost(0,0,0) =  1.0;  rightMatrixViewHost(0,0,1) =  2.0;
    rightMatrixViewHost(0,1,0) = -1.0;  rightMatrixViewHost(0,1,1) =  3.0;
    Kokkos::deep_copy(rightMatrixView, rightMatrixViewHost);
    
    auto expectedResultView = getView<Scalar>("result matrix", cellCount, spaceDim, spaceDim);
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
    
    Data<Scalar> A_data(leftMatrixView);
    Data<Scalar> B_data(rightMatrixView);
    auto actualResultData = Data<Scalar>::allocateMatMatResult(transposeA, A_data, transposeB, B_data);
    
    TEST_EQUALITY(       3,  actualResultData.rank());
    TEST_EQUALITY(cellCount, actualResultData.extent_int(0));
    TEST_EQUALITY(spaceDim,  actualResultData.extent_int(1));
    TEST_EQUALITY(spaceDim,  actualResultData.extent_int(2));
    
    actualResultData.storeMatMat(transposeA, A_data, transposeB, B_data);
    
    testFloatingEquality3(expectedResultView, actualResultData, relTol, absTol, out, success);
    
    printView(actualResultData.getUnderlyingView3(), out);
  }

// #pragma mark Data: MatMatExplicitIdentity_PDD
/** \brief Data provides matrix-matrix multiplication support.  This method checks correctness of the computed mat-mat for several cases involving 3x3 identity matrices.  Here, the notional dimensions (C,P,D,D) differ from the stored dimensions of (P,D,D).  We test each possible transpose combination.
*/
TEUCHOS_UNIT_TEST( Data, MatMatExplicitIdentity_PDD ) // (P,D,D) underlying; notionally (C,P,D,D)
{
  // tests that identity * identity = identity
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
  
  Data<Scalar> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
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
    auto actualResultData = Data<Scalar>::allocateMatMatResult(transposeA, explicitIdentityMatrix, transposeB, explicitIdentityMatrix);
    
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
/** \brief Data provides matrix-matrix multiplication support.  This method checks correctness of the computed mat-mat for a case involving one 3x3 matrix that has a 2x2 upper left block, and diagonal entry in the (3,3) position, and one 3x3 matrix that is entirely diagonal.  Here, the notional dimensions (C,D,D) match the stored dimensions.
*/
  TEUCHOS_UNIT_TEST( Data, MatMatBlockPlusDiagonal )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;
    
    using Scalar = double;
    const int spaceDim = 3;
    const int leftLastNonDiagonal = 1;
    const int cellCount = 1;
    auto leftMatrixView = getView<Scalar>("left matrix", cellCount, (leftLastNonDiagonal+1) * (leftLastNonDiagonal+1) + (spaceDim - leftLastNonDiagonal - 1));
    auto leftMatrixViewHost = Kokkos::create_mirror(leftMatrixView);
    int entryIndex = 0;
    // Block:
    leftMatrixViewHost(0,entryIndex++) =  1.0;  leftMatrixViewHost(0,entryIndex++) =  2.0;
    leftMatrixViewHost(0,entryIndex++) = -1.0;  leftMatrixViewHost(0,entryIndex++) =  3.0;
    // Diagonal:
    leftMatrixViewHost(0,entryIndex++) = 3.0;
    Kokkos::deep_copy(leftMatrixView, leftMatrixViewHost);
    
    const int rightLastNonDiagonal = 0;
    auto rightMatrixView = getView<Scalar>("right matrix", cellCount, (rightLastNonDiagonal+1) * (rightLastNonDiagonal+1) + (spaceDim - rightLastNonDiagonal - 1));
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
    
    Data<Scalar> A( leftMatrixView, rank, extents, variationTypes, leftLastNonDiagonal);
    Data<Scalar> B(rightMatrixView, rank, extents, variationTypes, rightLastNonDiagonal);
    
    auto expectedResultView = getView<Scalar>("result matrix", cellCount, spaceDim, spaceDim);
    auto expectedResultViewHost = Kokkos::create_mirror(expectedResultView);
    
    const int cellOrdinal = 0;
    for (int i=0; i<spaceDim; i++)
    {
      for (int j=0; j<spaceDim; j++)
      {
        Scalar result =  0;
        for (int k=0; k<spaceDim; k++)
        {
          Scalar left  = transposeA ? A(cellOrdinal,k,i) : A(cellOrdinal,i,k);
          Scalar right = transposeB ? B(cellOrdinal,j,k) : B(cellOrdinal,k,j);
          result += left * right;
        }
        expectedResultViewHost(cellOrdinal, i, j) = result;
      }
    }
    
    Kokkos::deep_copy(expectedResultView, expectedResultViewHost);
    
    auto actualResultData = Data<Scalar>::allocateMatMatResult(transposeA, A, transposeB, B);
    
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
  
  // test statically that Data supports all 7 rank operators
  static_assert(supports_rank<Data<double>,1>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double>,2>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double>,3>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double>,4>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double>,5>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double>,6>::value, "Data is expected to support up to rank 7");
  static_assert(supports_rank<Data<double>,7>::value, "Data is expected to support up to rank 7");
} // anonymous namespace
