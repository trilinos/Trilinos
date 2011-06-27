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

