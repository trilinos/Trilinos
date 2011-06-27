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
#include "Teuchos_FunctionObjectXMLConverterDB.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos{

/**
 * Serialization Tests for subtraction functions
 */
TEUCHOS_UNIT_TEST(Teuchos_Functions, SubtractionTests){
  RCP<SubtractionFunction<int> > intTester = rcp(
    new SubtractionFunction<int>(10));

  XMLObject subFuncXML = FunctionObjectXMLConverterDB::convertFunctionObject(
    intTester);

  std::string type = subFuncXML.getRequired(
    FunctionObjectXMLConverter::getTypeAttributeName());
  TEST_ASSERT(type == intTester->getTypeAttributeValue() );
  int operand = subFuncXML.getRequired<int>(
    SimpleFunctionXMLConverter<int>::getOperandAttributeName());
  TEST_ASSERT(operand == intTester->getModifiyingOperand());

  RCP<FunctionObject> readIn = 
    FunctionObjectXMLConverterDB::convertXML(subFuncXML);
  RCP<SubtractionFunction<int> > readInCasted = 
    rcp_dynamic_cast<SubtractionFunction<int> >(readIn);
  TEST_ASSERT(readInCasted.get() != NULL);
  TEST_ASSERT(
    readInCasted->getModifiyingOperand()
    ==
    intTester->getModifiyingOperand());
}

TEUCHOS_UNIT_TEST(Teuchos_Functions, AdditionTests){
  RCP<AdditionFunction<int> > intTester = rcp(
    new AdditionFunction<int>(10));

  XMLObject addFuncXML = FunctionObjectXMLConverterDB::convertFunctionObject(
    intTester);

  std::string type = addFuncXML.getRequired(
    FunctionObjectXMLConverter::getTypeAttributeName());
  TEST_ASSERT(type == intTester->getTypeAttributeValue() );
  int operand = addFuncXML.getRequired<int>(
    SimpleFunctionXMLConverter<int>::getOperandAttributeName());
  TEST_ASSERT(operand == intTester->getModifiyingOperand());

  RCP<FunctionObject> readIn = 
    FunctionObjectXMLConverterDB::convertXML(addFuncXML);
  RCP<AdditionFunction<int> > readInCasted = 
    rcp_dynamic_cast<AdditionFunction<int> >(readIn);
  TEST_ASSERT(readInCasted.get() != NULL);
  TEST_ASSERT(
    readInCasted->getModifiyingOperand()
    ==
    intTester->getModifiyingOperand());
}

TEUCHOS_UNIT_TEST(Teuchos_Functions, MultiplicationTests){
  RCP<MultiplicationFunction<int> > intTester = rcp(
    new MultiplicationFunction<int>(10));

  XMLObject multiFuncXML = FunctionObjectXMLConverterDB::convertFunctionObject(
    intTester);

  std::string type = multiFuncXML.getRequired(
    FunctionObjectXMLConverter::getTypeAttributeName());
  TEST_ASSERT(type == intTester->getTypeAttributeValue() );
  int operand = multiFuncXML.getRequired<int>(
    SimpleFunctionXMLConverter<int>::getOperandAttributeName());
  TEST_ASSERT(operand == intTester->getModifiyingOperand());

  RCP<FunctionObject> readIn = 
    FunctionObjectXMLConverterDB::convertXML(multiFuncXML);
  RCP<MultiplicationFunction<int> > readInCasted = 
    rcp_dynamic_cast<MultiplicationFunction<int> >(readIn);
  TEST_ASSERT(readInCasted.get() != NULL);
  TEST_ASSERT(
    readInCasted->getModifiyingOperand()
    ==
    intTester->getModifiyingOperand());
}

TEUCHOS_UNIT_TEST(Teuchos_Functions, DivisionTests){
  RCP<DivisionFunction<int> > intTester = rcp(
    new DivisionFunction<int>(10));

  XMLObject divisFuncXML = FunctionObjectXMLConverterDB::convertFunctionObject(
    intTester);

  std::string type = divisFuncXML.getRequired(
    FunctionObjectXMLConverter::getTypeAttributeName());
  TEST_ASSERT(type == intTester->getTypeAttributeValue() );
  int operand = divisFuncXML.getRequired<int>(
    SimpleFunctionXMLConverter<int>::getOperandAttributeName());
  TEST_ASSERT(operand == intTester->getModifiyingOperand());

  RCP<FunctionObject> readIn = 
    FunctionObjectXMLConverterDB::convertXML(divisFuncXML);
  RCP<DivisionFunction<int> > readInCasted = 
    rcp_dynamic_cast<DivisionFunction<int> >(readIn);
  TEST_ASSERT(readInCasted.get() != NULL);
  TEST_ASSERT(
    readInCasted->getModifiyingOperand()
    ==
    intTester->getModifiyingOperand());
}



} //namespace Teuchos

