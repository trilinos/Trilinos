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

