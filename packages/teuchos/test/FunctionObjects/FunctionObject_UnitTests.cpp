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

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardFunctionObjects.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos{

/**
 * Tests for subtraction functions
 */
TEUCHOS_UNIT_TEST(Teuchos_Functions, SubtractionTests){
  SubtractionFunction<int> intTester(10);
  TEST_ASSERT(intTester.runFunction(10) == 0);
  TEST_ASSERT(intTester.runFunction(-10) == -20);

  SubtractionFunction<double> doubleTester(5.5);
  TEST_ASSERT(doubleTester.runFunction(10) == 4.5);
}

/**
 * Tests for addition functions
 */
TEUCHOS_UNIT_TEST(Teuchos_Functions, AdditionTests){
  AdditionFunction<int> intTester(10);
  TEST_ASSERT(intTester.runFunction(10) == 20);
  TEST_ASSERT(intTester.runFunction(-10) == 0);

  AdditionFunction<double> doubleTester(5.5);
  TEST_ASSERT(doubleTester.runFunction(10) == 15.5);
}
  
/**
 * Tests for multiplication functions
 */
TEUCHOS_UNIT_TEST(Teuchos_Functions, MultiplicationTests){
  MultiplicationFunction<int> intTester(10);
  TEST_ASSERT(intTester.runFunction(10) == 100);
  TEST_ASSERT(intTester.runFunction(-10) == -100);

  MultiplicationFunction<double> doubleTester(5.5);
  TEST_ASSERT(doubleTester.runFunction(10) == 55);
}
  
/**
 * Tests for division functions
 */
TEUCHOS_UNIT_TEST(Teuchos_Functions, DivisionTests){
  DivisionFunction<int> intTester(10);
  TEST_ASSERT(intTester.runFunction(100) == 10);
  TEST_ASSERT(intTester.runFunction(-10) == -1);

  DivisionFunction<double> doubleTester(2);
  TEST_ASSERT(doubleTester.runFunction(7.5) == 3.75);
}
  


} //namespace Teuchos

