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
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_StandardConditionXMLConverters.hpp"
#include "Teuchos_XMLConditionExceptions.hpp"

#include "Teuchos_XMLParameterListTestHelpers.hpp"


namespace Teuchos{


/**
 * Test all the conditions
 */
TEUCHOS_UNIT_TEST(Teuchos_Conditions, StringConditionSerialization){
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
      testList.getEntryRCP(paramName2), conditionVal2));
  
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

  writeParameterListToXmlOStream(testList, out, depSheet1);

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
    rcp_dynamic_cast<const StringCondition>(simpleReadInDep->getCondition(), true);
  TEST_ASSERT(nonnull(simpleReadInCon));


  RCP<ConditionVisualDependency> complexReadInDep = 
    rcp_dynamic_cast<ConditionVisualDependency>(
      *(depSheetIn->getDependenciesForParameter(readInDependee2)->begin()));
  TEST_EQUALITY(
    complexReadInDep->getCondition()->getTypeAttributeValue(),
    DummyObjectGetter<StringCondition>::getDummyObject()->getTypeAttributeValue());
  RCP<const StringCondition> complexReadInCon = 
    rcp_dynamic_cast<const StringCondition>(complexReadInDep->getCondition(), true);
  TEST_ASSERT(nonnull(complexReadInCon));
 
    
  TEST_COMPARE_ARRAYS(
    simpleReadInCon->getValueList(), simpleStringCon->getValueList());
  TEST_COMPARE_ARRAYS(
    complexReadInCon->getValueList(), complexStringCon->getValueList());

}

TEUCHOS_UNIT_TEST(Teuchos_Conditions, BoolConditionSerialization){
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
    rcp_dynamic_cast<const BoolCondition>(simpleReadInDep->getCondition(), true);
  TEST_ASSERT(nonnull(simpleReadInCon));

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Conditions, NumberConditionSerialization, T){
  ConditionXMLConverterDB::printKnownConverters(out);
  std::string paramName1 = "T param";
  std::string paramName2 = "T param 2";
  std::string dependent1Name = "dependent1";
  std::string dependent2Name = "dependent2";
  T paramValue = ScalarTraits< T >::one();
  T ten = 10 * ScalarTraits< T >::one();
  std::string dependentValue = "hi there!";
  ParameterList testList("Condition Test List");
  testList.set(paramName1, paramValue); 
  testList.set(paramName2, paramValue); 
  testList.set(dependent1Name, dependentValue);
  testList.set(dependent2Name, dependentValue);

  RCP<NumberCondition< T > > numberCon = 
    rcp(new NumberCondition< T >(testList.getEntryRCP(paramName1)));

  RCP<SubtractionFunction< T > > funcTester = 
    rcp(new SubtractionFunction< T >(ten));

  RCP<NumberCondition< T > > numberFuncCon = 
    rcp(new NumberCondition< T >(testList.getEntryRCP(paramName2), funcTester));
  
  RCP<ConditionVisualDependency> numberConDep = 
    rcp(new ConditionVisualDependency(
      numberCon, 
      testList.getEntryRCP(dependent1Name)));

  RCP<ConditionVisualDependency> funcNumberConDep = 
    rcp(new ConditionVisualDependency(
      numberFuncCon, 
      testList.getEntryRCP(dependent2Name)));

  RCP<DependencySheet> depSheet1 = rcp(new DependencySheet);
  depSheet1->addDependency(numberConDep);
  depSheet1->addDependency(funcNumberConDep);

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
    DummyObjectGetter<NumberCondition< T > >::getDummyObject()->getTypeAttributeValue());
  RCP<const NumberCondition< T > > simpleReadInCon = 
    rcp_dynamic_cast<const NumberCondition< T > >(simpleReadInDep->getCondition(), true);
  TEST_ASSERT(nonnull(simpleReadInCon));


  RCP<ConditionVisualDependency> funcReadInDep = 
    rcp_dynamic_cast<ConditionVisualDependency>(
      *(depSheetIn->getDependenciesForParameter(readInDependee2)->begin()));
  TEST_ASSERT(funcReadInDep != null);

  RCP<const NumberCondition< T > > funcReadInCon = 
    rcp_dynamic_cast<const NumberCondition< T > >(funcReadInDep->getCondition());

  TEST_ASSERT(funcReadInCon != null);

  RCP<const SubtractionFunction< T > > funcReadInFunc = 
    rcp_dynamic_cast<const SubtractionFunction< T > >(
      funcReadInCon->getFunctionObject());
  TEST_ASSERT(funcReadInFunc != null);
  TEST_EQUALITY(
    funcReadInFunc->getModifiyingOperand(),
    funcTester->getModifiyingOperand());


}

#define NUMBER_PARAM_TYPE_TEST( T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Conditions, NumberConditionSerialization, T )

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long ulong;

NUMBER_PARAM_TYPE_TEST(int)
NUMBER_PARAM_TYPE_TEST(uint)
NUMBER_PARAM_TYPE_TEST(short)
NUMBER_PARAM_TYPE_TEST(ushort)
NUMBER_PARAM_TYPE_TEST(long)
NUMBER_PARAM_TYPE_TEST(ulong)
NUMBER_PARAM_TYPE_TEST(float)
NUMBER_PARAM_TYPE_TEST(double)
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long int llint;
typedef unsigned long long int ullint;
NUMBER_PARAM_TYPE_TEST(llint)
NUMBER_PARAM_TYPE_TEST(ullint)
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Conditions, BoolLogicConditionSerialization, BinCondition){
  ConditionXMLConverterDB::printKnownConverters(out);
  std::string paramName1 = "bool param1";
  std::string paramName2 = "bool param2";
  std::string dependent1Name = "dependent1";
  bool paramValue1 = true;
  bool paramValue2 = false;
  std::string dependentValue = "hi there!";
  ParameterList testList("Condition Test List");
  testList.set(paramName1, paramValue1); 
  testList.set(paramName2, paramValue2); 
  testList.set(dependent1Name, dependentValue);
  RCP<BoolCondition> boolCon1 = 
    rcp(new BoolCondition(testList.getEntryRCP(paramName1)));
  RCP<BoolCondition> boolCon2 = 
    rcp(new BoolCondition(testList.getEntryRCP(paramName1)));

  Condition::ConstConditionList conList = 
    tuple<RCP<const Condition> >(boolCon1, boolCon2);

  RCP< BinCondition > binCon = rcp(new BinCondition (conList));
  
  RCP<ConditionVisualDependency> binConDep = 
    rcp(new ConditionVisualDependency(
      binCon, 
      testList.getEntryRCP(dependent1Name)));

  RCP<DependencySheet> depSheet1 = rcp(new DependencySheet);
  depSheet1->addDependency(binConDep);

  RCP<DependencySheet> depSheetIn = rcp(new DependencySheet);
  RCP<ParameterList> readinList = 
    writeThenReadPL(testList, depSheet1, depSheetIn);
  
  RCP<ParameterEntry> readInDependee1 = readinList->getEntryRCP(paramName1);
  RCP<ParameterEntry> readInDependee2 = readinList->getEntryRCP(paramName2);

  RCP<ConditionVisualDependency> readInDep1 = 
    rcp_dynamic_cast<ConditionVisualDependency>(
      *(depSheetIn->getDependenciesForParameter(readInDependee1)->begin()));
  RCP<ConditionVisualDependency> readInDep2 = 
    rcp_dynamic_cast<ConditionVisualDependency>(
      *(depSheetIn->getDependenciesForParameter(readInDependee1)->begin()));
  TEST_EQUALITY(readInDep1.get(), readInDep1.get());
  TEST_EQUALITY(
    readInDep1->getCondition()->getTypeAttributeValue(),
    DummyObjectGetter< BinCondition >::getDummyObject()->getTypeAttributeValue());
  RCP<const BinCondition > readInCon = 
    rcp_dynamic_cast<const BinCondition >(readInDep1->getCondition(), true);
  TEST_ASSERT(nonnull(readInCon));

  Condition::ConstConditionList readInConList = readInCon->getConditions();
  TEST_ASSERT(readInConList.size() ==2);

}

#define BIN_CON_TEST( BinCondition ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Conditions, BoolLogicConditionSerialization, BinCondition)

BIN_CON_TEST(AndCondition)
BIN_CON_TEST(OrCondition)
BIN_CON_TEST(EqualsCondition)

TEUCHOS_UNIT_TEST(Teuchos_Conditions, NotConditionSerialization){
  ConditionXMLConverterDB::printKnownConverters(out);
  std::string paramName1 = "bool param1";
  std::string dependent1Name = "dependent1";
  bool paramValue1 = true;
  std::string dependentValue = "hi there!";
  ParameterList testList("Condition Test List");
  testList.set(paramName1, paramValue1); 
  testList.set(dependent1Name, dependentValue);
  RCP<BoolCondition> boolCon1 = 
    rcp(new BoolCondition(testList.getEntryRCP(paramName1)));


  RCP<NotCondition> notCon = rcp(new NotCondition(boolCon1));
  
  RCP<ConditionVisualDependency> notConDep = 
    rcp(new ConditionVisualDependency(
      notCon, 
      testList.getEntryRCP(dependent1Name)));

  RCP<DependencySheet> depSheet1 = rcp(new DependencySheet);
  depSheet1->addDependency(notConDep);

  RCP<DependencySheet> depSheetIn = rcp(new DependencySheet);
  RCP<ParameterList> readinList = 
    writeThenReadPL(testList, depSheet1, depSheetIn);
  
  RCP<ParameterEntry> readInDependee1 = readinList->getEntryRCP(paramName1);

  RCP<ConditionVisualDependency> readInDep1 = 
    rcp_dynamic_cast<ConditionVisualDependency>(
      *(depSheetIn->getDependenciesForParameter(readInDependee1)->begin()));
  TEST_EQUALITY(
    readInDep1->getCondition()->getTypeAttributeValue(),
    DummyObjectGetter<NotCondition>::getDummyObject()->getTypeAttributeValue());
  RCP<const NotCondition> readInCon = 
    rcp_dynamic_cast<const NotCondition>(readInDep1->getCondition(), true);
  TEST_ASSERT(nonnull(readInCon));
}

TEUCHOS_UNIT_TEST(Teuchos_Conditions, ConditionSerializationExceptions){
  ConditionXMLConverterDB::printKnownConverters(out);
  RCP<DependencySheet> depSheet = rcp(new DependencySheet);


  TEST_THROW(RCP<ParameterList> missingParameterList = 
    getParametersFromXmlFile(
      "MissingParameterEntryDefinition.xml", depSheet),
    MissingParameterEntryDefinitionException);
    
  RCP<ParameterEntry> notInListParam = rcp(new ParameterEntry(3.0));
  RCP<NumberCondition<double> > doubleCon =
    rcp(new NumberCondition<double>(notInListParam));

  NumberConditionConverter<double> doubleConConverter;
  XMLParameterListWriter::EntryIDsMap emptyMap;
  XMLObject toWriteTo;
  TEST_THROW(doubleConConverter.fromConditiontoXML(doubleCon, emptyMap),
    MissingParameterEntryDefinitionException);

  TEST_THROW(RCP<ParameterList> missingValuesList = 
    getParametersFromXmlFile(
      "MissingValuesTag.xml", depSheet),
    MissingValuesTagException);


}


} // namespace Teuchos
