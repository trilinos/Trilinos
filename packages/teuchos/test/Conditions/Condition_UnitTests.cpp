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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardConditions.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardConditions.hpp"

namespace Teuchos{

/**
 * Test all the conditions
 */
TEUCHOS_UNIT_TEST(Teuchos_Conditions, testConditions){
	//Settin up initial list
	RCP<ParameterList> testingList = rcp(new ParameterList("Condition Testing List"));

	/*
	 * Testing for string condition
	 */
	Array<std::string> validValues(tuple<std::string>("mountain dew", "pepsi", "coke", "fanta"));
	RCP<StringValidator> stringVali1 = rcp(new StringValidator(validValues));

	testingList->set("string param", "fanta", "parameter for testing string conditions", stringVali1);

	StringCondition::ValueList conValues1(tuple<std::string>("pepsi", "coke"));
	RCP<StringCondition> stringCon1 = rcp( new StringCondition(testingList->getEntryRCP("string param"), conValues1));
	TEST_ASSERT(!stringCon1->isConditionTrue());
	testingList->set("string param", "coke");
	TEST_ASSERT(stringCon1->isConditionTrue());

	/*
	 * Testing for number condition
	 */
	testingList->set("double param", 5.0, "parameter for testing number conditions");

	RCP<NumberCondition<double> > numberCon1 = 
    rcp( new NumberCondition<double>(testingList->getEntryRCP("double param")));
	TEST_ASSERT(numberCon1->isConditionTrue());
	testingList->set("double param", -1.0);
	TEST_ASSERT(!numberCon1->isConditionTrue());
RCP<SubtractionFunction<double> > doubleTesterFunc = rcp( new SubtractionFunction<double>(100));
	RCP<NumberCondition<double> > numberCon2 = 
    rcp( new NumberCondition<double>(testingList->getEntryRCP("double param"), doubleTesterFunc));
	TEST_ASSERT(!numberCon2->isConditionTrue());
	testingList->set("double param", 101.0);
	TEST_ASSERT(numberCon2->isConditionTrue());

	/*
	 * Testing bool conditions
	 */
	testingList->set("bool param", true, "parameter for testing bool conditions");

	RCP<BoolCondition> boolCon1 = rcp( new BoolCondition(testingList->getEntryRCP("bool param")));
	TEST_ASSERT(boolCon1->isConditionTrue());
	testingList->set("bool param", false);
	TEST_ASSERT(!boolCon1->isConditionTrue());

	/*
	 * Test Not condition
	 */
	RCP<NotCondition> notCon1 = rcp(new NotCondition(numberCon1));
	TEST_ASSERT(!notCon1->isConditionTrue());
	testingList->set("double param", -1.0);
	TEST_ASSERT(notCon1->isConditionTrue());

	/*
	 * Test And condition
	 */
	Condition::ConstConditionList conList1(tuple<RCP<const Condition> >(stringCon1, boolCon1));
	RCP<AndCondition> andCon1 = rcp(new AndCondition(conList1));
	TEST_ASSERT(!andCon1->isConditionTrue());
	testingList->set("bool param", true);
	TEST_ASSERT(andCon1->isConditionTrue());

	/*
	 * Testing or condition
	 */
	testingList->set("bool param", false);
	RCP<OrCondition> orCon1 = rcp(new OrCondition(conList1));
	TEST_ASSERT(orCon1->isConditionTrue());
	testingList->set("string param", "fanta");

	/*
	 * Testing equal condition
	 */
	RCP<EqualsCondition> equalsCon1 = rcp(new EqualsCondition(conList1));
	TEST_ASSERT(equalsCon1->isConditionTrue());
	testingList->set("bool param", true);
	TEST_ASSERT(!equalsCon1->isConditionTrue());
}

//Test getters and setters
TEUCHOS_UNIT_TEST(Teuchos_Conditions, testConditionGetterAndSetters){
	//Settin up initial list
	RCP<ParameterList> testingList = rcp(new ParameterList("Condition Testing List"));

	Array<std::string> validValues(tuple<std::string>("mountain dew", "pepsi", "coke", "fanta"));
	RCP<StringValidator> stringVali1 = rcp(new StringValidator(validValues));

	testingList->set("string param", "fanta", "parameter for testing string conditions", stringVali1);

	StringCondition::ValueList conValues1(tuple<std::string>("pepsi", "coke"));
	RCP<StringCondition> stringCon1 = rcp( new StringCondition(testingList->getEntryRCP("string param"), conValues1));
	Dependency::ConstParameterEntryList stringParameters = stringCon1->getAllParameters();
	TEST_ASSERT(stringParameters.size() == 1);
	TEST_ASSERT(stringParameters.find(testingList->getEntryRCP("string param")) != stringParameters.end());

	/*
	 * Testing for number condition
	 */
	testingList->set("double param", 5.0, "parameter for testing number conditions");

	RCP<NumberCondition<double> > numberCon1 = rcp( new NumberCondition<double>(testingList->getEntryRCP("double param")));
	Dependency::ConstParameterEntryList numberParameters = numberCon1->getAllParameters();
	TEST_ASSERT(numberParameters.size() == 1);
	TEST_ASSERT(numberParameters.find(testingList->getEntryRCP("double param")) != numberParameters.end());

	/*
	 * Testing bool conditions
	 */
	testingList->set("bool param", true, "parameter for testing bool conditions");

	RCP<BoolCondition> boolCon1 = rcp( new BoolCondition(testingList->getEntryRCP("bool param")));
	Dependency::ConstParameterEntryList boolParameters = boolCon1->getAllParameters();
	TEST_ASSERT(boolParameters.size() == 1);
	TEST_ASSERT(boolParameters.find(testingList->getEntryRCP("bool param")) != boolParameters.end());

	/*
	 * Test Not condition
	 */
	RCP<NotCondition> notCon1 = rcp(new NotCondition(numberCon1));
	Dependency::ConstParameterEntryList notParameters = notCon1->getAllParameters();
	TEST_ASSERT(notParameters.size() == 1);
	TEST_ASSERT(notParameters.find(testingList->getEntryRCP("double param")) != notParameters.end());

	/*
	 * Test And condition
	 */
	Condition::ConstConditionList conList1(tuple<RCP<const Condition> >(stringCon1, boolCon1));
	RCP<AndCondition> andCon1 = rcp(new AndCondition(conList1));
	Dependency::ConstParameterEntryList andParameters = andCon1->getAllParameters();
	TEST_ASSERT(andParameters.size() == 2);
	TEST_ASSERT(andParameters.find(testingList->getEntryRCP("string param")) != andParameters.end());
	TEST_ASSERT(andParameters.find(testingList->getEntryRCP("bool param")) != andParameters.end());

	/*
	 * Testing or condition
	 */
	RCP<OrCondition> orCon1 = rcp(new OrCondition(conList1));
	Dependency::ConstParameterEntryList orParameters = orCon1->getAllParameters();
	TEST_ASSERT(orParameters.size() == 2);
	TEST_ASSERT(orParameters.find(testingList->getEntryRCP("string param")) != orParameters.end());
	TEST_ASSERT(orParameters.find(testingList->getEntryRCP("bool param")) != orParameters.end());

	/*
	 * Testing Equsl condition
	 */
	Condition::ConstConditionList conList2(tuple<RCP<const Condition> >(numberCon1, boolCon1));
	RCP<EqualsCondition> equalsCon1 = rcp(new EqualsCondition(conList2));
	Dependency::ConstParameterEntryList equalsParameters = equalsCon1->getAllParameters();
	TEST_ASSERT(equalsParameters.size() == 2);
	TEST_ASSERT(equalsParameters.find(testingList->getEntryRCP("double param")) != equalsParameters.end());
	TEST_ASSERT(equalsParameters.find(testingList->getEntryRCP("bool param")) != equalsParameters.end());

	/*
	 * Testing BoolLogicCondition add
	 */
	equalsCon1->addCondition(orCon1);
	Dependency::ConstParameterEntryList equalsParameters2 = equalsCon1->getAllParameters();
	TEST_ASSERT(equalsParameters2.size() == 3);
	TEST_ASSERT(equalsParameters2.find(testingList->getEntryRCP("string param")) != equalsParameters2.end());
	TEST_ASSERT(equalsParameters2.find(testingList->getEntryRCP("double param")) != equalsParameters2.end());
	TEST_ASSERT(equalsParameters2.find(testingList->getEntryRCP("bool param")) != equalsParameters2.end());

}

//Test that exceptions get thrown when they should.
TEUCHOS_UNIT_TEST(Teuchos_Conditions, testConditionException){
	//Settin up initial list
	RCP<ParameterList> testingList = rcp(new ParameterList("Condition Testing List"));
	testingList->set("double param",1.0);
	testingList->set("string param", "awesome");
	RCP<ParameterList> testingList2 = rcp(new ParameterList("Condition Testing List"));
	testingList2->set("bool param", true);

	TEST_THROW(BoolCondition boolCon1(testingList->getEntryRCP("bool param")), InvalidConditionException);
	TEST_THROW(StringCondition stringCon1(testingList->getEntryRCP("double param"), "coke"), InvalidConditionException);
	TEST_THROW(BoolCondition boolCon1(testingList->getEntryRCP("double param")), InvalidConditionException);
	Condition::ConstConditionList conList1;
	TEST_THROW(AndCondition andCon1(conList1), InvalidConditionException);
	RCP<const Condition> con1;
	TEST_THROW(NotCondition notCon1(con1), InvalidConditionException);
}

} //namespace Teuchos
