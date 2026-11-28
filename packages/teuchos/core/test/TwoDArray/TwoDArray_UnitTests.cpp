// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TwoDArray.hpp"



namespace Teuchos{

TwoDArray<int> getSimpleTestTwoDArray(){
  TwoDArray<int> simpleArray(2,2);
  simpleArray(0,0) =1;
  simpleArray(0,1) =2;
  simpleArray(1,0) =3;
  simpleArray(1,1) =4;
  return simpleArray;
}


/**
 * Test all the validator dependencies.
 */
TEUCHOS_UNIT_TEST(Teuchos_TwoDArrays, simpleTest){
  TwoDArray<int> simpleArray(3,2);
  simpleArray[0][0] =1;
  simpleArray[0][1] =2;
  simpleArray[1][0] =3;
  simpleArray[1][1] =4;
  simpleArray[2][0] =5;
  simpleArray[2][1] =6;

  TEST_EQUALITY_CONST(simpleArray[0][0],1)
  TEST_EQUALITY_CONST(simpleArray[0][1],2)
  TEST_EQUALITY_CONST(simpleArray[1][0],3)
  TEST_EQUALITY_CONST(simpleArray[1][1],4)
  TEST_EQUALITY_CONST(simpleArray[2][0],5)
  TEST_EQUALITY_CONST(simpleArray[2][1],6)

  TwoDArray<int> const& simpleArrayConst = simpleArray;

  TEST_EQUALITY_CONST(simpleArrayConst[0][0],1)
  TEST_EQUALITY_CONST(simpleArrayConst[0][1],2)
  TEST_EQUALITY_CONST(simpleArrayConst[1][0],3)
  TEST_EQUALITY_CONST(simpleArrayConst[1][1],4)
  TEST_EQUALITY_CONST(simpleArrayConst[2][0],5)
  TEST_EQUALITY_CONST(simpleArrayConst[2][1],6)

  TEST_EQUALITY_CONST(simpleArray(0,0),1)
  TEST_EQUALITY_CONST(simpleArray(0,1),2)
  TEST_EQUALITY_CONST(simpleArray(1,0),3)
  TEST_EQUALITY_CONST(simpleArray(1,1),4)
  TEST_EQUALITY_CONST(simpleArray(2,0),5)
  TEST_EQUALITY_CONST(simpleArray(2,1),6)

  TEST_EQUALITY_CONST(simpleArray.getNumRows(), 3)
  TEST_EQUALITY_CONST(simpleArray.getNumCols(), 2)
  Array<int> oneDArray = tuple<int>(1,2,3,4,5,6);	
  TEST_COMPARE_ARRAYS(oneDArray, simpleArray.getDataArray())

}

TEUCHOS_UNIT_TEST(Teuchos_TwoDArrays, stringFunctions){
  TwoDArray<int> simpleArray = getSimpleTestTwoDArray();
  std::string stringRep = TwoDArray<int>::toString(simpleArray);
  TwoDArray<int> convertedArray = TwoDArray<int>::fromString(stringRep);
  TEST_EQUALITY(simpleArray, convertedArray)

  std::string badStringRep = "4x4:{1.0,1.0}";
  TEST_THROW(TwoDArray<int>::fromString(badStringRep),
    InvalidArrayStringRepresentation)
}

TEUCHOS_UNIT_TEST(Teuchos_TwoDArrays, fromStringRejectsMalformedDimensions){
  /*
   * Regression coverage for parsing dimensions from untrusted serialized
   * TwoDArray strings.  These cases must be rejected before constructing a
   * TwoDArray with dimensions that are inconsistent with the backing Array.
   */

  // Product wraps to 4 on 64-bit size_t:
  //   (2^62 + 1) * 4 == 2^64 + 4 == 4 mod 2^64
  const std::string wrappedProductString =
    "4611686018427387905x4:{1,2,3,4}";
  TEST_THROW(TwoDArray<int>::fromString(wrappedProductString),
    InvalidArrayStringRepresentation)

  // Signed dimension tokens must not be accepted as unsigned size_t values.
  const std::string negativeRowString =
    "-1x4:{1,2,3,4}";
  TEST_THROW(TwoDArray<int>::fromString(negativeRowString),
    InvalidArrayStringRepresentation)

  const std::string negativeColString =
    "1x-4:{1,2,3,4}";
  TEST_THROW(TwoDArray<int>::fromString(negativeColString),
    InvalidArrayStringRepresentation)

  // Missing dimensions delimiter.
  const std::string missingDimensionsDelimiterString =
    "2:{1,2}";
  TEST_THROW(TwoDArray<int>::fromString(missingDimensionsDelimiterString),
    InvalidArrayStringRepresentation)

  // Missing metadata separator between dimensions and payload.
  const std::string missingMetadataSeparatorString =
    "2x2";
  TEST_THROW(TwoDArray<int>::fromString(missingMetadataSeparatorString),
    InvalidArrayStringRepresentation)

  // Valid control case.
  const TwoDArray<int> validArray =
    TwoDArray<int>::fromString("2x2:{1,2,3,4}");
  TEST_EQUALITY_CONST(validArray.getNumRows(), 2)
  TEST_EQUALITY_CONST(validArray.getNumCols(), 2)
  TEST_EQUALITY_CONST(validArray(0,0), 1)
  TEST_EQUALITY_CONST(validArray(0,1), 2)
  TEST_EQUALITY_CONST(validArray(1,0), 3)
  TEST_EQUALITY_CONST(validArray(1,1), 4)
}

TEUCHOS_UNIT_TEST(Teuchos_TwoDArrays, emptyTest){
  TwoDArray<int> emptyArray;
  TEST_EQUALITY_CONST(emptyArray.getNumRows(), 0)
  TEST_EQUALITY_CONST(emptyArray.getNumCols(), 0)
  TEST_EQUALITY_CONST(emptyArray.getDataArray().size(), 0)
  TEST_ASSERT(emptyArray.isEmpty());
}

TEUCHOS_UNIT_TEST(Teuchos_TwoDArrays, streamTests){
  TwoDArray<int> simpleArray = getSimpleTestTwoDArray();
  std::stringstream ss;
  ss << simpleArray;
  TwoDArray<int> readArray;
  std::istringstream instream(ss.str());
  instream >> readArray;
  TEST_EQUALITY(simpleArray, readArray);
}

TEUCHOS_UNIT_TEST(Teuchos_TwoDArray, clearTest){
  TwoDArray<int> simpleArray = getSimpleTestTwoDArray();

  simpleArray.clear();
  TEST_ASSERT(simpleArray.isEmpty());
}

TEUCHOS_UNIT_TEST(Teuchos_TwoDArray, resizeTest){
  TwoDArray<int> simpleArray = getSimpleTestTwoDArray();

  simpleArray.resizeRows(4);
  TEST_EQUALITY_CONST(simpleArray.getNumRows(), 4);
  TEST_EQUALITY_CONST(simpleArray.getNumCols(), 2);
  TEST_EQUALITY_CONST(simpleArray(3,1), 0);
  TEST_EQUALITY_CONST(simpleArray(1,1), 4);

  simpleArray.resizeRows(2);
  TEST_EQUALITY_CONST(simpleArray.getNumRows(), 2);
  TEST_EQUALITY_CONST(simpleArray.getNumCols(), 2);
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW(simpleArray(3,1), RangeError);
#endif
  TEST_EQUALITY_CONST(simpleArray(1,1), 4);

  simpleArray.resizeCols(4);
  TEST_EQUALITY_CONST(simpleArray.getNumCols(), 4);
  TEST_EQUALITY_CONST(simpleArray.getNumRows(), 2);
  TEST_EQUALITY_CONST(simpleArray(1,3), 0);
  TEST_EQUALITY_CONST(simpleArray(1,1), 4);

  simpleArray.resizeCols(2);
  TEST_EQUALITY_CONST(simpleArray.getNumCols(), 2);
  TEST_EQUALITY_CONST(simpleArray.getNumRows(), 2);
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW(simpleArray(1,3), RangeError);
#endif
  TEST_EQUALITY_CONST(simpleArray(1,1), 4);

}

TEUCHOS_UNIT_TEST(Teuchos_TwoDArray, symmetryTest){
  TwoDArray<int> simpleArray = getSimpleTestTwoDArray();
  TEST_ASSERT(!simpleArray.isSymmetrical());
  simpleArray.setSymmetrical(true);
  TEST_ASSERT(simpleArray.isSymmetrical());

}

TEUCHOS_UNIT_TEST(Teuchos_TwoDArray, symmetrySerialization){
  TwoDArray<int> simpleArray = getSimpleTestTwoDArray();
  simpleArray.setSymmetrical(true);
  std::string arrayString = TwoDArray<int>::toString(simpleArray);
  TwoDArray<int> readIn = TwoDArray<int>::fromString(arrayString);
  TEST_ASSERT(readIn.isSymmetrical());
}


} //namespace Teuchos

