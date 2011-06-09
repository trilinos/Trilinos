// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TwoDArray.hpp"



namespace Teuchos{

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
  TwoDArray<int> simpleArray(2,2);
  simpleArray[0][0] =1;
  simpleArray[0][1] =2;
  simpleArray[1][0] =3;
  simpleArray[1][1] =4;
  std::string stringRep = TwoDArray<int>::toString(simpleArray);
  TwoDArray<int> convertedArray = TwoDArray<int>::fromString(stringRep);
  TEST_EQUALITY(simpleArray, convertedArray)
}

TEUCHOS_UNIT_TEST(Teuchos_TwoDArrays, emptyTest){
  TwoDArray<int> emptyArray;
  TEST_EQUALITY_CONST(emptyArray.getNumRows(), 0)
  TEST_EQUALITY_CONST(emptyArray.getNumCols(), 0)
  TEST_EQUALITY_CONST(emptyArray.getDataArray().size(), 0)

}

} //namespace Teuchos

