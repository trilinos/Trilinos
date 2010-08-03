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

#include "Teuchos_LocalTestingHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardConditions.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardConditions.hpp"

double doubleTesterFunc(double argument){
	return argument-100.0;
}

/**
 * Test all the conditions
 */
int testConditions(Teuchos::FancyOStream &out){
	bool success = true;
	//Settin up initial list
	Teuchos::RCP<Teuchos::ParameterList> testingList = Teuchos::rcp(new Teuchos::ParameterList("Condition Testing List"));

	/*
	 * Testing for string condition
	 */
	Teuchos::Array<std::string> validValues(Teuchos::tuple<std::string>("mountain dew", "pepsi", "coke", "fanta"));
	Teuchos::RCP<Teuchos::StringValidator> stringVali1 = Teuchos::rcp(new Teuchos::StringValidator(validValues));

	testingList->set("string param", "fanta", "parameter for testing string conditions", stringVali1);

	Teuchos::StringCondition::ValueList conValues1(Teuchos::tuple<std::string>("pepsi", "coke"));
	Teuchos::RCP<Teuchos::StringCondition> stringCon1 = rcp( new Teuchos::StringCondition("string param", testingList, conValues1));
	TEST_ASSERT(!stringCon1->isConditionTrue());
	testingList->set("string param", "coke");
	TEST_ASSERT(stringCon1->isConditionTrue());
	Teuchos::RCP<Teuchos::StringCondition> stringCon2 = Teuchos::rcp( new Teuchos::StringCondition("string param", testingList, conValues1, false));
	testingList->set("string param", "fanta");
	TEST_ASSERT(stringCon2->isConditionTrue());
	testingList->set("string param", "coke");
	TEST_ASSERT(!stringCon2->isConditionTrue());

	/*
	 * Testing for number condition
	 */
	testingList->set("double param", 5.0, "parameter for testing number conditions");

	Teuchos::RCP<Teuchos::NumberCondition<double> > numberCon1 = Teuchos::rcp( new Teuchos::NumberCondition<double>("double param", testingList, true));
	TEST_ASSERT(numberCon1->isConditionTrue());
	testingList->set("double param", -1.0);
	TEST_ASSERT(!numberCon1->isConditionTrue());

	Teuchos::RCP<Teuchos::NumberCondition<double> > numberCon2 = Teuchos::rcp( new Teuchos::NumberCondition<double>("double param", testingList, doubleTesterFunc, false));
	TEST_ASSERT(numberCon2->isConditionTrue());
	testingList->set("double param", 101.0);
	TEST_ASSERT(!numberCon2->isConditionTrue());

	/*
	 * Testing bool conditions
	 */
	testingList->set("bool param", true, "parameter for testing bool conditions");

	Teuchos::RCP<Teuchos::BoolCondition> boolCon1 = Teuchos::rcp( new Teuchos::BoolCondition("bool param", testingList));
	TEST_ASSERT(boolCon1->isConditionTrue());
	testingList->set("bool param", false);
	TEST_ASSERT(!boolCon1->isConditionTrue());

	Teuchos::RCP<Teuchos::BoolCondition> boolCon2 = Teuchos::rcp( new Teuchos::BoolCondition("bool param", testingList, false));
	TEST_ASSERT(boolCon2->isConditionTrue());
	testingList->set("bool param", true);
	TEST_ASSERT(!boolCon2->isConditionTrue());

	/*
	 * Test Not condition
	 */
	Teuchos::RCP<Teuchos::NotCondition> notCon1 = Teuchos::rcp(new Teuchos::NotCondition(numberCon1));
	TEST_ASSERT(!notCon1->isConditionTrue());
	testingList->set("double param", -1.0);
	TEST_ASSERT(notCon1->isConditionTrue());

	/*
	 * Test And condition
	 */
	Teuchos::Condition::ConditionList conList1(Teuchos::tuple<Teuchos::RCP<Teuchos::Condition> >(stringCon1, boolCon1));
	Teuchos::RCP<Teuchos::AndCondition> andCon1 = Teuchos::rcp(new Teuchos::AndCondition(conList1));
	TEST_ASSERT(andCon1->isConditionTrue());
	Teuchos::Condition::ConditionList conList2(Teuchos::tuple<Teuchos::RCP<Teuchos::Condition> >(stringCon1, boolCon2));
	Teuchos::RCP<Teuchos::AndCondition> andCon2 = Teuchos::rcp(new Teuchos::AndCondition(conList2));
	TEST_ASSERT(!andCon2->isConditionTrue());
	Teuchos::Condition::ConditionList conList3(Teuchos::tuple<Teuchos::RCP<Teuchos::Condition> >(stringCon2, boolCon2));
	Teuchos::RCP<Teuchos::AndCondition> andCon3 = Teuchos::rcp(new Teuchos::AndCondition(conList3));
	TEST_ASSERT(!andCon3->isConditionTrue());

	/*
	 * Testing or condition
	 */
	Teuchos::RCP<Teuchos::OrCondition> orCon1 = Teuchos::rcp(new Teuchos::OrCondition(conList1));
	TEST_ASSERT(orCon1->isConditionTrue());
	Teuchos::RCP<Teuchos::OrCondition> orCon2 = Teuchos::rcp(new Teuchos::OrCondition(conList2));
	TEST_ASSERT(orCon2->isConditionTrue());
	Teuchos::RCP<Teuchos::OrCondition> orCon3 = Teuchos::rcp(new Teuchos::OrCondition(conList3));
	TEST_ASSERT(!orCon3->isConditionTrue());

	/*
	 * Testing equal condition
	 */
	Teuchos::RCP<Teuchos::EqualsCondition> equalsCon1 = Teuchos::rcp(new Teuchos::EqualsCondition(conList1));
	TEST_ASSERT(equalsCon1->isConditionTrue());
	Teuchos::RCP<Teuchos::EqualsCondition> equalsCon2 = Teuchos::rcp(new Teuchos::EqualsCondition(conList2));
	TEST_ASSERT(!equalsCon2->isConditionTrue());
	Teuchos::RCP<Teuchos::EqualsCondition> equalsCon3 = Teuchos::rcp(new Teuchos::EqualsCondition(conList3));
	TEST_ASSERT(equalsCon3->isConditionTrue());

	return (success ? 0:1);
}

//Test getters and setters
int testConditionGetterAndSetters(Teuchos::FancyOStream &out){
	bool success = true;
	//Settin up initial list
	Teuchos::RCP<Teuchos::ParameterList> testingList = Teuchos::rcp(new Teuchos::ParameterList("Condition Testing List"));

	Teuchos::Array<std::string> validValues(Teuchos::tuple<std::string>("mountain dew", "pepsi", "coke", "fanta"));
	Teuchos::RCP<Teuchos::StringValidator> stringVali1 = Teuchos::rcp(new Teuchos::StringValidator(validValues));

	testingList->set("string param", "fanta", "parameter for testing string conditions", stringVali1);

	Teuchos::StringCondition::ValueList conValues1(Teuchos::tuple<std::string>("pepsi", "coke"));
	Teuchos::RCP<Teuchos::StringCondition> stringCon1 = Teuchos::rcp( new Teuchos::StringCondition("string param", testingList, conValues1));
	TEST_ASSERT(stringCon1->getType() == Teuchos::Condition::ParamCon);
	Teuchos::Dependency::ParameterParentMap stringParameters = stringCon1->getAllParameters();
	TEST_ASSERT(stringParameters.size() == 1);
	TEST_ASSERT(stringParameters.find("string param") != stringParameters.end());

	/*
	 * Testing for number condition
	 */
	testingList->set("double param", 5.0, "parameter for testing number conditions");

	Teuchos::RCP<Teuchos::NumberCondition<double> > numberCon1 = Teuchos::rcp( new Teuchos::NumberCondition<double>("double param", testingList, true));
	TEST_ASSERT(numberCon1->getType() == Teuchos::Condition::ParamCon);
	Teuchos::Dependency::ParameterParentMap numberParameters = numberCon1->getAllParameters();
	TEST_ASSERT(numberParameters.size() == 1);
	TEST_ASSERT(numberParameters.find("double param") != numberParameters.end());

	/*
	 * Testing bool conditions
	 */
	testingList->set("bool param", true, "parameter for testing bool conditions");

	Teuchos::RCP<Teuchos::BoolCondition> boolCon1 = Teuchos::rcp( new Teuchos::BoolCondition("bool param", testingList));
	TEST_ASSERT(boolCon1->getType() == Teuchos::Condition::ParamCon);
	Teuchos::Dependency::ParameterParentMap boolParameters = boolCon1->getAllParameters();
	TEST_ASSERT(boolParameters.size() == 1);
	TEST_ASSERT(boolParameters.find("bool param") != boolParameters.end());

	/*
	 * Test Not condition
	 */
	Teuchos::RCP<Teuchos::NotCondition> notCon1 = Teuchos::rcp(new Teuchos::NotCondition(numberCon1));
	TEST_ASSERT(notCon1->getType() == Teuchos::Condition::NotCon);
	Teuchos::Dependency::ParameterParentMap notParameters = notCon1->getAllParameters();
	TEST_ASSERT(notParameters.size() == 1);
	TEST_ASSERT(notParameters.find("double param") != notParameters.end());

	/*
	 * Test And condition
	 */
	Teuchos::Condition::ConditionList conList1(Teuchos::tuple<Teuchos::RCP<Teuchos::Condition> >(stringCon1, boolCon1));
	Teuchos::RCP<Teuchos::AndCondition> andCon1 = Teuchos::rcp(new Teuchos::AndCondition(conList1));
	TEST_ASSERT(andCon1->getType() == Teuchos::Condition::BinLogicCon);
	Teuchos::Dependency::ParameterParentMap andParameters = andCon1->getAllParameters();
	TEST_ASSERT(andParameters.size() == 2);
	TEST_ASSERT(andParameters.find("string param") != andParameters.end());
	TEST_ASSERT(andParameters.find("bool param") != andParameters.end());

	/*
	 * Testing or condition
	 */
	Teuchos::RCP<Teuchos::OrCondition> orCon1 = Teuchos::rcp(new Teuchos::OrCondition(conList1));
	TEST_ASSERT(orCon1->getType() == Teuchos::Condition::BinLogicCon);
	Teuchos::Dependency::ParameterParentMap orParameters = orCon1->getAllParameters();
	TEST_ASSERT(orParameters.size() == 2);
	TEST_ASSERT(orParameters.find("string param") != orParameters.end());
	TEST_ASSERT(orParameters.find("bool param") != orParameters.end());

	/*
	 * Testing Equsl condition
	 */
	Teuchos::Condition::ConditionList conList2(Teuchos::tuple<Teuchos::RCP<Teuchos::Condition> >(numberCon1, boolCon1));
	Teuchos::RCP<Teuchos::EqualsCondition> equalsCon1 = Teuchos::rcp(new Teuchos::EqualsCondition(conList2));
	TEST_ASSERT(equalsCon1->getType() == Teuchos::Condition::BinLogicCon);
	Teuchos::Dependency::ParameterParentMap equalsParameters = equalsCon1->getAllParameters();
	TEST_ASSERT(equalsParameters.size() == 2);
	TEST_ASSERT(equalsParameters.find("double param") != equalsParameters.end());
	TEST_ASSERT(equalsParameters.find("bool param") != equalsParameters.end());

	/*
	 * Testing BinaryLogicCondition add
	 */
	equalsCon1->addCondition(orCon1);
	Teuchos::Dependency::ParameterParentMap equalsParameters2 = equalsCon1->getAllParameters();
	TEST_ASSERT(equalsParameters2.size() == 3);
	TEST_ASSERT(equalsParameters2.find("double param") != equalsParameters2.end());
	TEST_ASSERT(equalsParameters2.find("bool param") != equalsParameters2.end());
	TEST_ASSERT(equalsParameters2.find("string param") != equalsParameters2.end());



	return (success ? 0:1);
}

//Test that exceptions get thrown when they should.
int testConditionException(Teuchos::FancyOStream &out){
	bool success = true;
	//Settin up initial list
	Teuchos::RCP<Teuchos::ParameterList> testingList = Teuchos::rcp(new Teuchos::ParameterList("Condition Testing List"));
	testingList->set("double param",1.0);
	testingList->set("string param", "awesome");
	Teuchos::RCP<Teuchos::ParameterList> testingList2 = Teuchos::rcp(new Teuchos::ParameterList("Condition Testing List"));
	testingList2->set("bool param", true);

	TEST_THROW(Teuchos::BoolCondition boolCon1("bool param", testingList), Teuchos::InvalidConditionException);
	TEST_THROW(Teuchos::StringCondition stringCon1("double param", testingList, "coke"), Teuchos::InvalidConditionException);
	TEST_THROW(Teuchos::NumberCondition<double> doubleCon1("string param", testingList, true), Teuchos::InvalidConditionException);
	TEST_THROW(Teuchos::BoolCondition boolCon1("double param", testingList), Teuchos::InvalidConditionException);
	Teuchos::Condition::ConditionList conList1;
	TEST_THROW(Teuchos::AndCondition andCon1(conList1), Teuchos::InvalidConditionException);
	Teuchos::RCP<Teuchos::Condition> con1;
	TEST_THROW(Teuchos::NotCondition notCon1(con1), Teuchos::InvalidConditionException);

	return (success ? 0:1);
}

int main(int argc, char* argv[]){
	bool success = true;
	Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
	if(testConditions(*out) == 1){
		success = false;
	}

	if(testConditionGetterAndSetters(*out) == 1){
		success = false;
	}

	return (success ? 0:1);
}


