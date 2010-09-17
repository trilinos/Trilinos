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
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_DependencySheet.hpp"
#include "Teuchos_StandardConditions.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DependencyXMLConverterDB.hpp"
#include "Teuchos_StandardDependencyXMLConverters.hpp"
#include "Teuchos_ParameterList.cpp"


namespace Teuchos{


typedef unsigned short int ushort;
typedef unsigned int uint;
typedef unsigned long int ulong;
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long int llint;
typedef unsigned long long int ullint;
#endif


TEUCHOS_UNIT_TEST(Teuchos_Dependencies, stringVisualDepTest){
  std::string dependee1 = "string param";
  std::string dependee2 = "string param2";
  std::string dependent1 = "dependent param1";
  std::string dependent2 = "dependent param2";
  ParameterList myDepList("String Visual Dep List");
  RCP<DependencySheet> myDepSheet = rcp(new DependencySheet);
  myDepList.set(dependee1, "val1");
  myDepList.set(dependee2, "val2");
  myDepList.set(dependent1, 1.0);
  myDepList.set(dependent2, 1.0);

  StringVisualDependency::ValueList valList1 = tuple<std::string>("val1");

  RCP<StringVisualDependency> basicStringVisDep = rcp(
    new StringVisualDependency(
      myDepList.getEntryRCP(dependee1),
      myDepList.getEntryRCP(dependent1),
      valList1));

  Dependency::ParameterEntryList dependentList;
  dependentList.insert(myDepList.getEntryRCP(dependent1));
  dependentList.insert(myDepList.getEntryRCP(dependent2));
  StringVisualDependency::ValueList valList2 = 
    tuple<std::string>("val1", "val2");

  RCP<StringVisualDependency> complexStringVisDep = rcp(
    new StringVisualDependency(
      myDepList.getEntryRCP(dependee2),
      dependentList,
      valList2,
      false));

  myDepSheet->addDependency(basicStringVisDep);
  myDepSheet->addDependency(complexStringVisDep);

  RCP<DependencySheet> readInDepSheet = rcp(new DependencySheet);

  XMLParameterListWriter plWriter;
  XMLObject xmlOut = plWriter.toXML(myDepList, myDepSheet);
  out << xmlOut.toString();

  RCP<ParameterList> readInList = 
    writeThenReadPL(myDepList, myDepSheet, readInDepSheet); 

  RCP<ParameterEntry> readinDependee1 = readInList->getEntryRCP(dependee1);
  RCP<ParameterEntry> readinDependent1 = readInList->getEntryRCP(dependent1);
  RCP<ParameterEntry> readinDependee2 = readInList->getEntryRCP(dependee2);
  RCP<ParameterEntry> readinDependent2 = readInList->getEntryRCP(dependent2);
  
  RCP<Dependency> readinDep1 =
    *(readInDepSheet->getDependenciesForParameter(readinDependee1)->begin());

  RCP<Dependency> readinDep2 =
    *(readInDepSheet->getDependenciesForParameter(readinDependee2)->begin());

  std::string stringVisXMLTag = 
    DummyObjectGetter<StringVisualDependency>::getDummyObject()->getTypeAttributeValue();

  TEST_ASSERT(readinDep1->getTypeAttributeValue() == stringVisXMLTag);
  TEST_ASSERT(readinDep2->getTypeAttributeValue() == stringVisXMLTag);

  TEST_ASSERT(readinDep1->getFirstDependee().get() == readinDependee1.get());
  TEST_ASSERT(readinDep1->getDependents().size() == 1);
  TEST_ASSERT((*readinDep1->getDependents().begin()).get() == readinDependent1.get());

  TEST_ASSERT(readinDep2->getFirstDependee().get() == readinDependee2.get());
  TEST_ASSERT(readinDep2->getDependents().size() == 2);
  TEST_ASSERT(
    readinDep2->getDependents().find(readinDependent1) 
    !=
    readinDep2->getDependents().end()
  );
  TEST_ASSERT(
    readinDep2->getDependents().find(readinDependent2)
    !=
    readinDep2->getDependents().end()
  );
    
  RCP<StringVisualDependency> castedDep1 =
    rcp_dynamic_cast<StringVisualDependency>(readinDep1);
  RCP<StringVisualDependency> castedDep2 =
    rcp_dynamic_cast<StringVisualDependency>(readinDep2);

  TEST_COMPARE_ARRAYS(
    castedDep1->getValues(), basicStringVisDep->getValues());
  TEST_COMPARE_ARRAYS(
    castedDep2->getValues(), complexStringVisDep->getValues());

  TEST_EQUALITY(castedDep1->getShowIf(), basicStringVisDep->getShowIf());
  TEST_EQUALITY(castedDep2->getShowIf(), complexStringVisDep->getShowIf());
  

}

TEUCHOS_UNIT_TEST(Teuchos_Dependencies, boolVisualDepTest){
  std::string dependee1 = "bool param";
  std::string dependee2 = "bool param2";
  std::string dependent1 = "dependent param1";
  std::string dependent2 = "dependent param2";
  ParameterList myDepList("Bool Visual Dep List");
  RCP<DependencySheet> myDepSheet = rcp(new DependencySheet);
  myDepList.set(dependee1, true);
  myDepList.set(dependee2, true);
  myDepList.set(dependent1, 1.0);
  myDepList.set(dependent2, 1.0);

  RCP<BoolVisualDependency> trueBoolVisDep = rcp(
    new BoolVisualDependency(
      myDepList.getEntryRCP(dependee1),
      myDepList.getEntryRCP(dependent1)));

  Dependency::ParameterEntryList dependentList;
  dependentList.insert(myDepList.getEntryRCP(dependent1));
  dependentList.insert(myDepList.getEntryRCP(dependent2));

  RCP<BoolVisualDependency> falseBoolVisDep = rcp(
    new BoolVisualDependency(
      myDepList.getEntryRCP(dependee2),
      dependentList,
      false));

  myDepSheet->addDependency(trueBoolVisDep);
  myDepSheet->addDependency(falseBoolVisDep);

  RCP<DependencySheet> readInDepSheet = rcp(new DependencySheet);

  XMLParameterListWriter plWriter;
  XMLObject xmlOut = plWriter.toXML(myDepList, myDepSheet);
  out << xmlOut.toString();

  RCP<ParameterList> readInList = 
    writeThenReadPL(myDepList, myDepSheet, readInDepSheet); 

  RCP<ParameterEntry> readinDependee1 = readInList->getEntryRCP(dependee1);
  RCP<ParameterEntry> readinDependent1 = readInList->getEntryRCP(dependent1);
  RCP<ParameterEntry> readinDependee2 = readInList->getEntryRCP(dependee2);
  RCP<ParameterEntry> readinDependent2 = readInList->getEntryRCP(dependent2);
  
  RCP<Dependency> readinDep1 =
    *(readInDepSheet->getDependenciesForParameter(readinDependee1)->begin());

  RCP<Dependency> readinDep2 =
    *(readInDepSheet->getDependenciesForParameter(readinDependee2)->begin());

  std::string boolVisXMLTag = 
    DummyObjectGetter<BoolVisualDependency>::getDummyObject()->getTypeAttributeValue();

  TEST_ASSERT(readinDep1->getTypeAttributeValue() == boolVisXMLTag);
  TEST_ASSERT(readinDep2->getTypeAttributeValue() == boolVisXMLTag);

  TEST_ASSERT(readinDep1->getFirstDependee().get() == readinDependee1.get());
  TEST_ASSERT(readinDep1->getDependents().size() == 1);
  TEST_ASSERT((*readinDep1->getDependents().begin()).get() == readinDependent1.get());

  TEST_ASSERT(readinDep2->getFirstDependee().get() == readinDependee2.get());
  TEST_ASSERT(readinDep2->getDependents().size() == 2);
  TEST_ASSERT(
    readinDep2->getDependents().find(readinDependent1) 
    !=
    readinDep2->getDependents().end()
  );
  TEST_ASSERT(
    readinDep2->getDependents().find(readinDependent2)
    !=
    readinDep2->getDependents().end()
  );
    
  RCP<BoolVisualDependency> castedDep1 =
    rcp_dynamic_cast<BoolVisualDependency>(readinDep1);
  RCP<BoolVisualDependency> castedDep2 =
    rcp_dynamic_cast<BoolVisualDependency>(readinDep2);

  TEST_EQUALITY(castedDep1->getShowIf(), trueBoolVisDep->getShowIf());
  TEST_EQUALITY(castedDep2->getShowIf(), falseBoolVisDep->getShowIf());
  

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Dependencies, numberVisualDepTest, T){
  std::string dependee1 = "num param";
  std::string dependee2 = "num param2";
  std::string dependent1 = "dependent param1";
  std::string dependent2 = "dependent param2";
  ParameterList myDepList("Number Visual Dep List");
  RCP<DependencySheet> myDepSheet = rcp(new DependencySheet);
  myDepList.set(dependee1, ScalarTraits<T>::one());
  myDepList.set(dependee2, ScalarTraits<T>::one());
  myDepList.set(dependent1, true);
  myDepList.set(dependent2, "vale");

  RCP<NumberVisualDependency< T > > simpleNumVisDep = rcp(
    new NumberVisualDependency< T >(
      myDepList.getEntryRCP(dependee1),
      myDepList.getEntryRCP(dependent1)));

  Dependency::ParameterEntryList dependentList;
  dependentList.insert(myDepList.getEntryRCP(dependent1));
  dependentList.insert(myDepList.getEntryRCP(dependent2));

  RCP<NumberVisualDependency< T > > complexNumVisDep = rcp(
    new NumberVisualDependency< T >(
      myDepList.getEntryRCP(dependee2),
      dependentList));

  myDepSheet->addDependency(simpleNumVisDep);
  myDepSheet->addDependency(complexNumVisDep);

  RCP<DependencySheet> readInDepSheet = rcp(new DependencySheet);

  XMLParameterListWriter plWriter;
  XMLObject xmlOut = plWriter.toXML(myDepList, myDepSheet);
  out << xmlOut.toString();

  RCP<ParameterList> readInList = 
    writeThenReadPL(myDepList, myDepSheet, readInDepSheet); 

  RCP<ParameterEntry> readinDependee1 = readInList->getEntryRCP(dependee1);
  RCP<ParameterEntry> readinDependent1 = readInList->getEntryRCP(dependent1);
  RCP<ParameterEntry> readinDependee2 = readInList->getEntryRCP(dependee2);
  RCP<ParameterEntry> readinDependent2 = readInList->getEntryRCP(dependent2);
  
  RCP<Dependency> readinDep1 =
    *(readInDepSheet->getDependenciesForParameter(readinDependee1)->begin());

  RCP<Dependency> readinDep2 =
    *(readInDepSheet->getDependenciesForParameter(readinDependee2)->begin());

  std::string numVisXMLTag = 
    DummyObjectGetter<NumberVisualDependency<T> >::getDummyObject()->getTypeAttributeValue();

  TEST_ASSERT(readinDep1->getTypeAttributeValue() == numVisXMLTag);
  TEST_ASSERT(readinDep2->getTypeAttributeValue() == numVisXMLTag);

  TEST_ASSERT(readinDep1->getFirstDependee().get() == readinDependee1.get());
  TEST_ASSERT(readinDep1->getDependents().size() == 1);
  TEST_ASSERT((*readinDep1->getDependents().begin()).get() == readinDependent1.get());

  TEST_ASSERT(readinDep2->getFirstDependee().get() == readinDependee2.get());
  TEST_ASSERT(readinDep2->getDependents().size() == 2);
  TEST_ASSERT(
    readinDep2->getDependents().find(readinDependent1) 
    !=
    readinDep2->getDependents().end()
  );
  TEST_ASSERT(
    readinDep2->getDependents().find(readinDependent2)
    !=
    readinDep2->getDependents().end()
  );
    
  RCP<NumberVisualDependency<T> > castedDep1 =
    rcp_dynamic_cast<NumberVisualDependency<T> >(readinDep1);
  RCP<NumberVisualDependency<T> > castedDep2 =
    rcp_dynamic_cast<NumberVisualDependency<T> >(readinDep2);

  TEST_EQUALITY(castedDep1->getShowIf(), simpleNumVisDep->getShowIf());
  TEST_EQUALITY(castedDep2->getShowIf(), complexNumVisDep->getShowIf());
  

}

#define NUMBER_VIS_TEST(T) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Dependencies, numberVisualDepTest, T)

NUMBER_VIS_TEST(int)
NUMBER_VIS_TEST(uint)
NUMBER_VIS_TEST(short)
NUMBER_VIS_TEST(ushort)
NUMBER_VIS_TEST(long)
NUMBER_VIS_TEST(ulong)
NUMBER_VIS_TEST(float)
NUMBER_VIS_TEST(double)
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
NUMBER_VIS_TEST(llint)
NUMBER_VIS_TEST(ullint)
#endif

} //namespace Teuchos

