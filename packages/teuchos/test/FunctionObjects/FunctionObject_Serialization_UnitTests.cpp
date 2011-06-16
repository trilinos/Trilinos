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

