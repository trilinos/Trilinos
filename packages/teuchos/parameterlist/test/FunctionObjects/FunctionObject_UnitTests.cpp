// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

