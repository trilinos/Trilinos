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
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardConditions.hpp"
#include "Teuchos_ConditionXMLConverterDB.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
namespace Teuchos{


/**
 * Test all the conditions
 */
TEUCHOS_UNIT_TEST(Teuchos_Conditions, testStringCondition){
  ConditionXMLConverterDB::printKnownConverters(out);
  std::string paramName1 = "string param";
  std::string paramName2 = "string param2";
  std::string dependent1Name = "dependent1";
  std::string dependent2Name = "dependent2";
  std::string paramValue = "cheese";
  StringCondition::ValueList conditionVal1 = tuple<std::string>("steve");
  StringCondition::ValueList conditionVal2 = 
    tuple<std::string>("steve", "blah", "your face");
  ParameterList testList("Condition Test List");
  testList.set(paramName1, paramValue); 
  testList.set(paramName2, paramValue);
  testList.set(dependent1Name, paramValue);
  testList.set(dependent2Name, paramValue);
  RCP<StringCondition> simpleStringCon = 
    rcp(new StringCondition(testList.getEntryRCP(paramName1), conditionVal1));
  RCP<StringCondition> complexStringCon = 
    rcp(new StringCondition(
      testList.getEntryRCP(paramName2), conditionVal2, false));
  
  RCP<ConditionVisualDependency> simpleConDep = 
    rcp(new ConditionVisualDependency(
      simpleStringCon, 
      testList.getEntryRCP(dependent1Name)));

  RCP<ConditionVisualDependency> complexConDep = 
    rcp(new ConditionVisualDependency(
      complexStringCon, 
      testList.getEntryRCP(dependent2Name)));
 
  RCP<DependencySheet> depSheet1 = rcp(new DependencySheet);
  depSheet1->addDependency(simpleConDep);
  depSheet1->addDependency(complexConDep);

  RCP<DependencySheet> depSheetIn = rcp(new DependencySheet);
  RCP<ParameterList> readinList = 
    writeThenReadPL(testList, depSheet1, depSheetIn);
  
  RCP<ParameterEntry> readInDependee1 = readinList->getEntryRCP(paramName1);
  RCP<ParameterEntry> readInDependee2 = readinList->getEntryRCP(paramName2);

  RCP<ConditionVisualDependency> simpleReadInDep = 
    rcp_dynamic_cast<ConditionVisualDependency>(
      *(depSheetIn->getDependenciesForParameter(readInDependee1)->begin()));
  TEST_EQUALITY(
    simpleReadInDep->getCondition()->getTypeAttributeValue(),
    DummyObjectGetter<StringCondition>::getDummyObject()->getTypeAttributeValue());
  RCP<const StringCondition> simpleReadInCon = 
    rcp_dynamic_cast<const StringCondition>(simpleReadInDep->getCondition());


  RCP<ConditionVisualDependency> complexReadInDep = 
    rcp_dynamic_cast<ConditionVisualDependency>(
      *(depSheetIn->getDependenciesForParameter(readInDependee2)->begin()));
  TEST_EQUALITY(
    complexReadInDep->getCondition()->getTypeAttributeValue(),
    DummyObjectGetter<StringCondition>::getDummyObject()->getTypeAttributeValue());
  RCP<const StringCondition> complexReadInCon = 
    rcp_dynamic_cast<const StringCondition>(complexReadInDep->getCondition());
 
    
  TEST_COMPARE_ARRAYS(
    simpleReadInCon->getValueList(), simpleStringCon->getValueList());
  TEST_COMPARE_ARRAYS(
    complexReadInCon->getValueList(), complexStringCon->getValueList());

  TEST_EQUALITY(
    simpleReadInCon->getWhenParamEqualsValue(), 
    simpleStringCon->getWhenParamEqualsValue());
  TEST_EQUALITY(
    complexReadInCon->getWhenParamEqualsValue(), 
    complexStringCon->getWhenParamEqualsValue());
}

TEUCHOS_UNIT_TEST(Teuchos_Conditions, testBoolCondition){
  ConditionXMLConverterDB::printKnownConverters(out);
  std::string paramName1 = "bool param";
  std::string dependent1Name = "dependent1";
  bool paramValue = true;
  std::string dependentValue = "hi there!";
  ParameterList testList("Condition Test List");
  testList.set(paramName1, paramValue); 
  testList.set(dependent1Name, dependentValue);
  RCP<BoolCondition> boolCon = 
    rcp(new BoolCondition(testList.getEntryRCP(paramName1)));
  
  RCP<ConditionVisualDependency> boolConDep = 
    rcp(new ConditionVisualDependency(
      boolCon, 
      testList.getEntryRCP(dependent1Name)));

  RCP<DependencySheet> depSheet1 = rcp(new DependencySheet);
  depSheet1->addDependency(boolConDep);

  RCP<DependencySheet> depSheetIn = rcp(new DependencySheet);
  RCP<ParameterList> readinList = 
    writeThenReadPL(testList, depSheet1, depSheetIn);
  
  RCP<ParameterEntry> readInDependee1 = readinList->getEntryRCP(paramName1);

  RCP<ConditionVisualDependency> simpleReadInDep = 
    rcp_dynamic_cast<ConditionVisualDependency>(
      *(depSheetIn->getDependenciesForParameter(readInDependee1)->begin()));
  TEST_EQUALITY(
    simpleReadInDep->getCondition()->getTypeAttributeValue(),
    DummyObjectGetter<BoolCondition>::getDummyObject()->getTypeAttributeValue());
  RCP<const BoolCondition> simpleReadInCon = 
    rcp_dynamic_cast<const BoolCondition>(simpleReadInDep->getCondition());

  TEST_EQUALITY(
    simpleReadInCon->getWhenParamEqualsValue(), 
    boolCon->getWhenParamEqualsValue());
}


} // namespace Teuchos
