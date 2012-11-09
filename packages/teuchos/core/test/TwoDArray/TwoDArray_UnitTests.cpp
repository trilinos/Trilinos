// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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

