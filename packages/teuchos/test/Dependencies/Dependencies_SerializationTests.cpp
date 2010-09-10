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


TEUCHOS_UNIT_TEST(Teuchos_Dependencies, stringVisualDepTest){
  std::string dependee1 = "string param";
  /*std::string dependee2 = "string param2";
  std::string dependent1 = "dependent param1";
  std::string dependent2 = "dependent param2";*/
  ParameterList myDepList("String Visual Dep List");
  //RCP<DependencySheet> myDepSheet = rcp(new DependencySheet);
  myDepList.set(dependee1, "val1");
  /*myDepList.set(dependee2, "val2");
  myDepList.set(dependent1, 1.0);
  myDepList.set(dependent2, 1.0);

  RCP<StringVisualDependency> basicStringVisDep = rcp(
    new StringVisualDependency(
      myDepList.getEntryRCP(dependee1),
      myDepList.getEntryRCP(dependent1),
      "val1"));

  Dependency::ParameterEntryList dependentList;
  dependentList.insert(myDepList.getEntryRCP(dependent1));
  dependentList.insert(myDepList.getEntryRCP(dependent2));
  StringVisualDependency::ValueList valList = 
    tuple<std::string>("val1", "val2");

  RCP<StringVisualDependency> complexStringVisDep = rcp(
    new StringVisualDependency(
      myDepList.getEntryRCP(dependee2),
      dependentList,
      valList));

  myDepSheet->addDependency(basicStringVisDep);
  myDepSheet->addDependency(complexStringVisDep);

  RCP<DependencySheet> readInDepSheet = rcp(new DependencySheet);

*/
 /* RCP<ParameterList> readInList = 
    writeThenReadPL(myDepList, myDepSheet, readInDepSheet); */
  RCP<ParameterList> readInList = 
    writeThenReadPL(myDepList);
 

  //TEST_ASSERT(
    //readInDepSheet->hasDependents(readInList->getEntryRCP(dependee1)));

}


} //namespace Teuchos

