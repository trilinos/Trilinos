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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_DependencySheet.hpp"
#include "Teuchos_StandardConditions.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos{

/**
 * Test all the validator dependencies.
 */
TEUCHOS_UNIT_TEST(Teuchos_Dependencies, testValiDeps){
	RCP<ParameterList> My_deplist = rcp(new ParameterList);
	RCP<DependencySheet> depSheet1 = rcp(new DependencySheet);

	/*
	 * Testing StringValidatorDependency
	 */
 	RCP<StringToIntegralParameterEntryValidator<int> >
   	stringFoodTypeValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
		tuple<std::string>( "Cheese", "Soda", "Chips" )
		,"Food Type"
		)
	);

	RCP<StringToIntegralParameterEntryValidator<int> >
    cheeseValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
   			tuple<std::string>( "Swiss", "American", "Super Awesome Cheese" )
			,"Food Selector"
			)
	);

	RCP<StringToIntegralParameterEntryValidator<int> >
	sodaValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( "Pepsi", "Coke", "Kurtis Cola", "Bad Cola" )
			,"Food Selector"
			)
		);

	RCP<StringToIntegralParameterEntryValidator<int> >
	chipsValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( "Lays", "Doritos", "Kurtis Super Awesome Brand" )
			,"Food Selector"
		)
	);

	StringValidatorDependency::ValueToValidatorMap testValidatorMap1;
	testValidatorMap1["Cheese"] = cheeseValidator;
	testValidatorMap1["Soda"] = sodaValidator;
	testValidatorMap1["Chips"] = chipsValidator;

	ParameterList stringValiDepList = My_deplist->sublist(
    "String Validator Dependency", false, 
    "String Validator Dependency testing list.");
	stringValiDepList.set(
    "Food Selector", "Swiss", "select the food you want", cheeseValidator);
	stringValiDepList.set(
    "Food Type", 
    "Cheese", 
    "String Validator Dependency Tester", 
    stringFoodTypeValidator);

	RCP<StringValidatorDependency> 
	stringValiDep = rcp(
		new StringValidatorDependency(
			stringValiDepList.getEntryRCP("Food Type"),
			stringValiDepList.getEntryRCP("Food Selector"),
			testValidatorMap1, 
			cheeseValidator
		)
	);

	depSheet1->addDependency(stringValiDep);
	
	TEST_NOTHROW(stringValiDepList.validateParameters(stringValiDepList));
	TEST_ASSERT(depSheet1->hasDependents(
    stringValiDepList.getEntryRCP("Food Type")));
	RCP<const DependencySheet::DepSet> stringValiDepSet = 
    depSheet1->getDependenciesForParameter(
      stringValiDepList.getEntryRCP("Food Type"));
	TEST_ASSERT(stringValiDepSet->size() == 1);
	stringValiDepList.set("Food Type","Soda");
	stringValiDep->evaluate();
	TEST_ASSERT(stringValiDepList.getEntry("Food Selector").validator()
    ==
    sodaValidator);
	TEST_THROW(stringValiDepList.validateParameters(stringValiDepList), 
    Exceptions::InvalidParameterValue);
	stringValiDepList.set("Food Selector", "Pepsi");
	TEST_NOTHROW(stringValiDepList.validateParameters(stringValiDepList));
 

	/*
	 * Tesing some different aspects of the StringValidatorDependency
	 */
	ParameterList 
	stringValiDepList2 = My_deplist->sublist(
		"String Validator Dependency (other validators)",
		false,
		"String validator testing"
	);

	RCP<StringToIntegralParameterEntryValidator<int> >
	stringRangeValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
		tuple<std::string>( "1-10", "10-33", "50-60" ),
		"Range selector"
		)
	);

	RCP<EnhancedNumberValidator<int> > range110Vali = 
	rcp(new EnhancedNumberValidator<int>(1,10));
	RCP<EnhancedNumberValidator<int> > range1033Vali = 
	rcp(new EnhancedNumberValidator<int>(10,33));
	RCP<EnhancedNumberValidator<int> > range5060Vali = 
	rcp(new EnhancedNumberValidator<int>(50,60));

	stringValiDepList2.set("Range selector", "1-10", 
    "selects the range to validate", stringRangeValidator);

	StringValidatorDependency::ValueToValidatorMap rangeValidatorMap1;
	rangeValidatorMap1["1-10"] = range110Vali;
	rangeValidatorMap1["10-33"] = range1033Vali;
	rangeValidatorMap1["50-60"] = range5060Vali;
	stringValiDepList2.set(
    "RangeValue", 3, "the value of the range", range110Vali);

	RCP<StringValidatorDependency> 
	stringValiDep2 = RCP<StringValidatorDependency>(
		new StringValidatorDependency(
			stringValiDepList2.getEntryRCP("Range selector"),
			stringValiDepList2.getEntryRCP("RangeValue"),
			rangeValidatorMap1, 
      range110Vali
		)
	);

	depSheet1->addDependency(stringValiDep2);

	TEST_NOTHROW(stringValiDepList2.validateParameters(stringValiDepList2));
	TEST_ASSERT(depSheet1->hasDependents(
    stringValiDepList2.getEntryRCP("Range selector")));
	RCP<const DependencySheet::DepSet> stringValiDepSet2 = 
    depSheet1->getDependenciesForParameter(
      stringValiDepList2.getEntryRCP("Range selector"));
	TEST_ASSERT(stringValiDepSet2->size() == 1);
	stringValiDepList2.set("Range selector","50-60");
	stringValiDep2->evaluate();
	TEST_ASSERT(stringValiDepList2.getEntry("RangeValue").validator() 
    ==
    range5060Vali);
	TEST_THROW(stringValiDepList2.validateParameters(stringValiDepList2), 
    Exceptions::InvalidParameterValue);
	stringValiDepList2.set("RangeValue", 55);
	TEST_NOTHROW(stringValiDepList2.validateParameters(stringValiDepList2));

	/*
	 * Testing the BoolValidatorDependency.
	 */
	ParameterList
	boolValidatorDepList = My_deplist->sublist(
		"Bool Validator Dependency List", 
		false,
		"Bool Validator Dependency testing list."
	);

	boolValidatorDepList.set("Use Validator?", 
    true, "truns the validator on and off");
	RCP<EnhancedNumberValidator<int> > basicVali = 
    rcp(new EnhancedNumberValidator<int>(1,10));
	RCP<EnhancedNumberValidator<int> > basicVali2 = 
    rcp(new EnhancedNumberValidator<int>());
	boolValidatorDepList.set("do I have a validator?", 
    4, "does it have a validator?", basicVali);

	RCP<BoolValidatorDependency> 
	boolValiDep = RCP<BoolValidatorDependency>(
		new BoolValidatorDependency(
			boolValidatorDepList.getEntryRCP("Use Validator?"),
			boolValidatorDepList.getEntryRCP("do I have a validator?"),
			basicVali, 
			basicVali2
		)
	);

	depSheet1->addDependency(boolValiDep);

	TEST_ASSERT(depSheet1->hasDependents(
    boolValidatorDepList.getEntryRCP("Use Validator?")));
	TEST_ASSERT(
    boolValidatorDepList.getEntry("do I have a validator?").validator()
    == 
    basicVali);
	TEST_NOTHROW(
    boolValidatorDepList.validateParameters(boolValidatorDepList));
	RCP<const DependencySheet::DepSet> boolValiDepSet = 
    depSheet1->getDependenciesForParameter(boolValidatorDepList.getEntryRCP(
      "Use Validator?"));
	TEST_ASSERT(boolValiDepSet->size() == 1);
	boolValidatorDepList.set("Use Validator?",false);
	boolValiDep->evaluate();
	TEST_ASSERT(
    boolValidatorDepList.getEntry("do I have a validator?").validator() 
    ==
    basicVali2);


	/*
	 * Testing the RangeValidatorDependency
	 */
	RCP<StringToIntegralParameterEntryValidator<int> >
	lowTempCheeseValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( "PepperJack", "Swiss", "American" ),
			"Cheese to Fondue"
		)
	);

	RCP<StringToIntegralParameterEntryValidator<int> >
	highTempCheeseValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( 
        "Munster", "Provalone", "Kurtis Super Awesome Cheese"),
			"Cheese to Fondue"
		)
	);

	RCP<StringToIntegralParameterEntryValidator<int> >
	defaultCheeseValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( 
        "Other cheese", "other cheese 1", "other cheese 3"),
			"Cheese to Fondue"
		)
	);

	ParameterList& 
	rangeValidatorDepList = My_deplist->sublist(
		"Range Validator Dependency List",
		false,
		"Range Validator Dependency testing list.\nWorking June 27th 2009"
	);
	rangeValidatorDepList.set(
    "Temperature",101.0, "The temperature of the fondue");
	rangeValidatorDepList.set(
    "Cheese to Fondue", "Swiss", 
    "The cheese we'll be using in our fondue pot.", lowTempCheeseValidator);
	RangeValidatorDependency<double>::RangeToValidatorMap tempranges;
	tempranges[std::pair<double,double>(100,200)] = lowTempCheeseValidator;
	tempranges[std::pair<double,double>(200,300)] = highTempCheeseValidator;
	RCP<RangeValidatorDependency<double> > 
	cheeseTempDep = RCP<RangeValidatorDependency<double> >(
		new RangeValidatorDependency<double>(
			rangeValidatorDepList.getEntryRCP("Temperature"),
			rangeValidatorDepList.getEntryRCP("Cheese to Fondue"),
			tempranges,
      defaultCheeseValidator
		)
	);
	depSheet1->addDependency(cheeseTempDep);

	TEST_ASSERT(depSheet1->hasDependents(
    rangeValidatorDepList.getEntryRCP("Temperature")));
	RCP<const DependencySheet::DepSet> rangeValiDepSet = 
    depSheet1->getDependenciesForParameter(
      rangeValidatorDepList.getEntryRCP("Temperature"));
	TEST_ASSERT(rangeValiDepSet->size() == 1);
	rangeValidatorDepList.set("Temperature",250.0);
	cheeseTempDep->evaluate();
	TEST_ASSERT(
    rangeValidatorDepList.getEntry("Cheese to Fondue").validator()
    == 
    highTempCheeseValidator);
	TEST_THROW(
    rangeValidatorDepList.validateParameters(rangeValidatorDepList), 
    Exceptions::InvalidParameterValue);
	rangeValidatorDepList.set("Cheese to Fondue", "Provalone");
	TEST_NOTHROW(
    rangeValidatorDepList.validateParameters(rangeValidatorDepList));
  rangeValidatorDepList.set("Temperature", 50.0);
  cheeseTempDep->evaluate();
  TEST_ASSERT(
    rangeValidatorDepList.getEntry("Cheese to Fondue").validator()
    ==
    defaultCheeseValidator
  );

}
  
/**
 * Testing all the visual dependencies
 */
TEUCHOS_UNIT_TEST(Teuchos_Dependencies, testVisualDeps){
	RCP<ParameterList> My_deplist = RCP<ParameterList>(new ParameterList);
	RCP<DependencySheet> depSheet1 = 
    RCP<DependencySheet>(new DependencySheet);
  /*
   * Two Simple NumberVisualDependency test
   */

	ParameterList
	simpleNumDepTestList = My_deplist->sublist(
		"NumberVisual Dependency List (double)", 
		false, 
		"Number visual Dependency testing list"
	);
		
	simpleNumDepTestList.set("Temperature",101.0);
	simpleNumDepTestList.set("Cheese to Fondue", "Swiss", "The cheese to fondue");
	simpleNumDepTestList.set("reverse param", "hello");

	RCP<NumberVisualDependency<double> > simpleNumDep = 
	RCP<NumberVisualDependency<double> >(
		new NumberVisualDependency<double>(
			simpleNumDepTestList.getEntryRCP("Temperature"),
			simpleNumDepTestList.getEntryRCP("Cheese to Fondue"),
      true
		)
	);
	RCP<NumberVisualDependency<double> > reverseNumDep = 
	RCP<NumberVisualDependency<double> >(
		new NumberVisualDependency<double>(
			simpleNumDepTestList.getEntryRCP("Temperature"),
			simpleNumDepTestList.getEntryRCP("reverse param"),
      false
		)
	);
  depSheet1->addDependency(simpleNumDep);
  depSheet1->addDependency(reverseNumDep);
	simpleNumDep->evaluate();
	reverseNumDep->evaluate();
	TEST_ASSERT(simpleNumDep->isDependentVisible());
	TEST_ASSERT(!reverseNumDep->isDependentVisible());
	simpleNumDepTestList.set("Temperature",-1.0);
	simpleNumDep->evaluate();
	reverseNumDep->evaluate();
	TEST_ASSERT(!simpleNumDep->isDependentVisible());
	TEST_ASSERT(reverseNumDep->isDependentVisible());


	/*
	 * complex Testing the NumberVisualDependency
	 */
	ParameterList
	doubleVisualDepList = My_deplist->sublist(
		"NumberVisual Dependency List (double)", 
		false, 
		"Number visual Dependency testing list"
	);
		
	doubleVisualDepList.set(
    "Temperature",101.0, "The temperature of the fondue");
	doubleVisualDepList.set(
    "Cheese to Fondue", "Swiss", "The cheese to fondue");
  doubleVisualDepList.set("reverse param", "hello");
  RCP<SubtractionFunction<double> > fondueFunc = rcp(new
    SubtractionFunction<double>(100));

	RCP<NumberVisualDependency<double> > fondueDep = 
	RCP<NumberVisualDependency<double> >(
		new NumberVisualDependency<double>(
			doubleVisualDepList.getEntryRCP("Temperature"),
			doubleVisualDepList.getEntryRCP("Cheese to Fondue"),
      true,
			fondueFunc
		)
	);
	RCP<NumberVisualDependency<double> > reverseFondueDep = 
	RCP<NumberVisualDependency<double> >(
		new NumberVisualDependency<double>(
			doubleVisualDepList.getEntryRCP("Temperature"),
			doubleVisualDepList.getEntryRCP("reverse param"),
      false,
			fondueFunc
		)
	);
	depSheet1->addDependency(fondueDep);
	depSheet1->addDependency(reverseFondueDep);
	fondueDep->evaluate();
  reverseFondueDep->evaluate();
	TEST_ASSERT(fondueDep->isDependentVisible());
	TEST_ASSERT(!reverseFondueDep->isDependentVisible());
	doubleVisualDepList.set("Temperature",99.0);
	fondueDep->evaluate();
  reverseFondueDep->evaluate();
	TEST_ASSERT(!fondueDep->isDependentVisible());
	TEST_ASSERT(reverseFondueDep->isDependentVisible());

	/*
	 * Testing the BoolVisualDependency
	 */
	ParameterList&
	boolVisDepList = My_deplist->sublist(
		"Bool Visual Dependency List", 
		false,
		"Bool Visual Dependency testing list."
	);
	boolVisDepList.set(
    "ShowPrecs", true, "Whether or not to should the Preciondtioner list");
	ParameterList
	Prec_List0 = boolVisDepList.sublist(
    "Preconditioner",false,"Sublist that defines the preconditioner.");
	Prec_List0.set("Type", "ILU", "The tpye of preconditioner to use");
	RCP<EnhancedNumberValidator<double> > droptolValidator = 
    rcp(new EnhancedNumberValidator<double>(0,10,1e-3));
	Prec_List0.set(
    "Drop Tolerance", 1e-3,
    "The tolerance below which entries from the "
    "factorization are left out of the factors.", droptolValidator);
	RCP<BoolVisualDependency> 
	precDep1 = RCP<BoolVisualDependency>(
		new BoolVisualDependency(
			boolVisDepList.getEntryRCP("ShowPrecs"),
			boolVisDepList.getEntryRCP("Preconditioner"),
			true
		)
	);
	depSheet1->addDependency(precDep1);
	precDep1->evaluate();
	TEST_ASSERT(precDep1->isDependentVisible());
	boolVisDepList.set("ShowPrecs", false);
	precDep1->evaluate();
	TEST_ASSERT(!precDep1->isDependentVisible());



	/*
	 * Testing the StringVisualDepenency
	 */
	ParameterList&
    stringVisDepList = My_deplist->sublist(
		"String Visual Dependency List",
		false,
		"String Visual Dependency testing list."
	);
	RCP<StringToIntegralParameterEntryValidator<int> >
	favCheeseValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( "Swiss", "American", "Cheder" ),
			"Favorite Cheese"
		)
	);
   
	stringVisDepList.set("Favorite Cheese", 
    "American", "Your favorite type of cheese", favCheeseValidator);
	RCP<EnhancedNumberValidator<int> > 
	swissValidator = rcp(new EnhancedNumberValidator<int>(0,10));
	stringVisDepList.set("Swiss rating", 0, 
    "How you rate swiss on a scale of 1 to 10", swissValidator);
	RCP<StringVisualDependency> 
	swissDep1 = RCP<StringVisualDependency>(
		new StringVisualDependency(
			stringVisDepList.getEntryRCP("Favorite Cheese"),
			stringVisDepList.getEntryRCP("Swiss rating"),
			"Swiss", 
			true
		)
	);
	depSheet1->addDependency(swissDep1);
	swissDep1->evaluate();
	TEST_ASSERT(!swissDep1->isDependentVisible());
	stringVisDepList.set("Favorite Cheese", "Swiss");
	swissDep1->evaluate();
	TEST_ASSERT(swissDep1->isDependentVisible());

	/*
	 * String Visual Tester with multiple values
	 */
	ParameterList multiStringVisDepList = My_deplist->sublist(
		"Multi String Visual Dependency List",
		false
	);
	RCP<StringToIntegralParameterEntryValidator<int> >
	favCheeseValidator2 = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( "Provalone", "Swiss", "American", "Cheder" ),
			"Favorite Cheese"
		)
	);
   
	multiStringVisDepList.set(
    "Favorite Cheese", "American", 
    "Your favorite type of cheese", favCheeseValidator2);
	multiStringVisDepList.set("Swiss rating", 0, 
    "How you rate swiss on a scale of 1 to 10", swissValidator);
	RCP<StringVisualDependency> 
	swissDep2 = RCP<StringVisualDependency>(
		new StringVisualDependency(
			multiStringVisDepList.getEntryRCP("Favorite Cheese"),
			multiStringVisDepList.getEntryRCP("Swiss rating"), 
			tuple<std::string>("Swiss", "Cheder"), 
			true
		)
	);
	depSheet1->addDependency(swissDep2);
	swissDep2->evaluate();
	TEST_ASSERT(!swissDep2->isDependentVisible());
	multiStringVisDepList.set("Favorite Cheese", "Cheder");
	swissDep2->evaluate();
	TEST_ASSERT(swissDep2->isDependentVisible());

	/*
	 * Another test of the NumberVisualDependency.
	 */
	ParameterList
    numberVisDepList = My_deplist->sublist(
		"Number Visual Dependency List", 
		false, 
		"Number Visual Dependency testing list."
	);
	numberVisDepList.set("Ice", 50, "Ice stuff");
	numberVisDepList.set("Room Temp", 10, "Room temperature");
  RCP<SubtractionFunction<int> > visFunc = rcp(new
    SubtractionFunction<int>(32));
	RCP<NumberVisualDependency<int> > 
	iceDep = RCP<NumberVisualDependency<int> >(
		new NumberVisualDependency<int>(
			numberVisDepList.getEntryRCP("Room Temp"),
			numberVisDepList.getEntryRCP("Ice"), 
      true,
			visFunc
		)
	);
	depSheet1->addDependency(iceDep);
	iceDep->evaluate();
	TEST_ASSERT(!iceDep->isDependentVisible());
	numberVisDepList.set("Room Temp", 33);
	iceDep->evaluate();
	TEST_ASSERT(iceDep->isDependentVisible());

	/*
	 * Test condition visual dependency
	 */
	RCP<ParameterList> conVisDepList = sublist(
    My_deplist,"Condition Visual Dependency List", false);
	conVisDepList->set("double param", 4.0, "double parameter");
	conVisDepList->set("bool param", true, "bool parameter");
	conVisDepList->set("string param", "blah", "a string parameter");
	RCP<NumberCondition<double> > numberCon = 
    rcp( new NumberCondition<double>(
      conVisDepList->getEntryRCP("double param")));
	RCP<BoolCondition> boolCon = 
    rcp(new BoolCondition(conVisDepList->getEntryRCP("bool param")));
	Condition::ConstConditionList conList = 
    tuple<RCP<const Condition> >(numberCon, boolCon);
	RCP<AndCondition> andCon = rcp(new AndCondition(conList));
	RCP<ConditionVisualDependency> conVisDep = 
    rcp(new ConditionVisualDependency(
      andCon, conVisDepList->getEntryRCP("string param"), true));
	depSheet1->addDependency(conVisDep);
	conVisDep->evaluate();
	TEST_ASSERT(conVisDep->isDependentVisible());
	conVisDepList->set("bool param", false);
	conVisDep->evaluate();
	TEST_ASSERT(!conVisDep->isDependentVisible());
}


/**
 * Test the ArrayLengthDependency.
 */
TEUCHOS_UNIT_TEST(Teuchos_Dependencies, testArrayLengthDep){
	RCP<ParameterList> My_deplist = RCP<ParameterList>(new ParameterList);
	RCP<DependencySheet> depSheet1 = 
    RCP<DependencySheet>(new DependencySheet);

	ParameterList
	numberArrayLengthDepList = My_deplist->sublist(
    "Number Array Length Dependency List", false,
    "Number Array Length Dependecy testing list.");
	numberArrayLengthDepList.set("Array Length", 10, "array length setter");
	Array<double> variableLengthArray(10,23.0);
	RCP<EnhancedNumberValidator<double> > 
	varLengthArrayVali = RCP<EnhancedNumberValidator<double> >(
  		new EnhancedNumberValidator<double>(10,50,4) 
	);
	numberArrayLengthDepList.set(
    "Variable Length Array", variableLengthArray, "variable length array",
	  RCP<ArrayNumberValidator<double> >(
      new ArrayNumberValidator<double>(varLengthArrayVali)));

	RCP<NumberArrayLengthDependency<int, double> >
	  arrayLengthDep(
  		new NumberArrayLengthDependency<int, double>(
			numberArrayLengthDepList.getEntryRCP("Array Length"),
			numberArrayLengthDepList.getEntryRCP("Variable Length Array") 
		)
	);
	depSheet1->addDependency(arrayLengthDep);
  Array<double> curArray = 
    numberArrayLengthDepList.get<Array<double> >("Variable Length Array");
	TEST_ASSERT(curArray.length() ==10);
	numberArrayLengthDepList.set("Array Length", 12);
	arrayLengthDep()->evaluate();
  curArray = 
    numberArrayLengthDepList.get<Array<double> >("Variable Length Array");
  out << curArray.length() << std::endl;
	TEST_ASSERT(curArray.length() ==12);
	numberArrayLengthDepList.set("Array Length", -1);
	TEST_THROW(arrayLengthDep()->evaluate(), 
    Exceptions::InvalidParameterValue);
}

/**
 * Tests the excpetions associated with Dependencies
 */
TEUCHOS_UNIT_TEST(Teuchos_Dependencies, testDepExceptions){
	RCP<ParameterList> list1 = RCP<ParameterList>(new ParameterList());
	RCP<ParameterList> list2 = RCP<ParameterList>(new ParameterList());

	list1->set("int parameter", 4, "int parameter");
	list1->set("double parameter", 6.0, "double parameter");
	list1->set("string parameter", "hahahaha", "string parameter");
	Array<double> doubleArray(10,23.0);
	list1->set("array parameter", doubleArray, "array parameter");
	list1->set("bool parameter", true, "bool parameter");

  RCP<AdditionFunction<int> > intFuncTester = rcp(new
    AdditionFunction<int>(10));
	TEST_THROW(RCP<NumberVisualDependency<int> > numValiDep = 
    rcp(
      new NumberVisualDependency<int>(
        list1->getEntryRCP("bool parameter"),
        list1->getEntryRCP("double parameter"), 
        true,
        intFuncTester)), 
    InvalidDependencyException);

	/*
	 * Testing StringVisualDepenendcy exceptions.
	 */
	RCP<StringVisualDependency> stringVisDep;
	TEST_THROW(stringVisDep = RCP<StringVisualDependency>(
    new StringVisualDependency(
      list1->getEntryRCP("double parameter"), 
      list1->getEntryRCP("int parameter"),
      "cheese", true)), 
    InvalidDependencyException);

	/*
	 * Testing BoolVisualDependency exceptions.
	 */
	TEST_THROW(RCP<BoolVisualDependency> boolVisDep = 
    RCP<BoolVisualDependency>(new BoolVisualDependency(
      list1->getEntryRCP("int parameter"), 
      list1->getEntryRCP("double parameter"), false)),
      InvalidDependencyException);

  /**
   * Tesint NumberArrayLengthDependency excpetions */
  RCP<NumberArrayLengthDependency<int, double> > numArrayLengthDep;
	TEST_THROW(numArrayLengthDep = 
      rcp(new NumberArrayLengthDependency<int, double>(
        list1->getEntryRCP("double parameter"), 
        list1->getEntryRCP("array parameter"))), 
      InvalidDependencyException);

	TEST_THROW(numArrayLengthDep = 
      rcp(new NumberArrayLengthDependency<int, double>(
        list1->getEntryRCP("int parameter"), 
        list1->getEntryRCP("double parameter"))), 
      InvalidDependencyException);

	/*
	 * Testing StringValidatorDependency exceptions.
	 */
	RCP<StringToIntegralParameterEntryValidator<int> >
    cheeseValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
   	  tuple<std::string>( "Swiss", "American", "Super Awesome Cheese"),
			"Food Selector"
		)
	);

	RCP<StringToIntegralParameterEntryValidator<int> >
	sodaValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( "Pepsi", "Coke", "Kurtis Cola", "Bad Cola" ),
			"Food Selector"
		)
	);

	RCP<StringToIntegralParameterEntryValidator<int> >
	chipsValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( "Lays", "Doritos", "Kurtis Super Awesome Brand"),
			"Food Selector"
		)
	);


	list1->set(
    "string 2 parameter", "Swiss", 
    "second string parameter", cheeseValidator);
	StringValidatorDependency::ValueToValidatorMap testValidatorMap1;
	testValidatorMap1["Cheese"] = cheeseValidator;
	testValidatorMap1["Soda"] = sodaValidator;
	testValidatorMap1["Chips"] = chipsValidator;
	TEST_THROW(RCP<StringValidatorDependency> stringValiDep = 
    RCP<StringValidatorDependency>(
      new StringValidatorDependency(
        list1->getEntryRCP("int parameter"),  
        list1->getEntryRCP("string 2 parameter"),
        testValidatorMap1)),
    InvalidDependencyException);
	RCP<EnhancedNumberValidator<int> > intVali = 
    rcp(new EnhancedNumberValidator<int>(0,20));
	testValidatorMap1["Candy"] = intVali;
	TEST_THROW(RCP<StringValidatorDependency> stringValiDep = 
    RCP<StringValidatorDependency>(
      new StringValidatorDependency(
        list1->getEntryRCP("string parameter"),
        list1->getEntryRCP("string 2 parameter"),
        testValidatorMap1)),
    InvalidDependencyException);

  StringValidatorDependency::ValueToValidatorMap emptyMap;
	TEST_THROW(RCP<StringValidatorDependency> stringValiDep = 
    RCP<StringValidatorDependency>(
      new StringValidatorDependency(
        list1->getEntryRCP("string parameter"),
        list1->getEntryRCP("string 2 parameter"),
        emptyMap)),
    InvalidDependencyException);
	
	/*
	 * Testing BoolValidatorDependency exceptions.
	 */
	RCP<EnhancedNumberValidator<double> > doubleVali1 = 
    rcp(new EnhancedNumberValidator<double>(0.0,20.0));
	RCP<EnhancedNumberValidator<double> > doubleVali2 =
    rcp(new EnhancedNumberValidator<double>(5.0,20.0));
	list1->set("double parameter", 6.0, "double parameter", doubleVali1);

	TEST_THROW(RCP<BoolValidatorDependency> boolValiDep = 
    RCP<BoolValidatorDependency>(
      new BoolValidatorDependency(
        list1->getEntryRCP("int parameter"),
        list1->getEntryRCP("double parameter"),
        doubleVali1, 
        doubleVali2)), 
    InvalidDependencyException);

	TEST_THROW(RCP<BoolValidatorDependency> boolValiDep = 
    RCP<BoolValidatorDependency>(
      new BoolValidatorDependency(
      list1->getEntryRCP("bool parameter"),
      list1->getEntryRCP("double parameter"), 
      intVali, 
      doubleVali2)), 
    InvalidDependencyException);

	TEST_THROW(RCP<BoolValidatorDependency> boolValiDep = 
    RCP<BoolValidatorDependency>(
      new BoolValidatorDependency(
        list1->getEntryRCP("bool parameter"),
        list1->getEntryRCP("double parameter"), 
        doubleVali1, 
        intVali)), 
    InvalidDependencyException);

	/*
	 * Testing RangeValidatorDependency exceptions.
	 */
	list1->set("Cheese to Fondue", "Swiss", "the cheese to fondue");
	RCP<StringToIntegralParameterEntryValidator<int> >
	lowTempCheeseValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>( "PepperJack", "Swiss", "American" ),
			"Cheese to Fondue"
		)
	);
	RCP<StringToIntegralParameterEntryValidator<int> >
	highTempCheeseValidator = rcp(
		new StringToIntegralParameterEntryValidator<int>(
			tuple<std::string>("Munster", "Provalone", 
        "Kurtis Super Awesome Cheese"),
			"Cheese to Fondue"
		)
	);

	list1->set(
    "Cheese to Fondue", "Swiss", "the cheese to fondue", 
    lowTempCheeseValidator);

	RangeValidatorDependency<double>::RangeToValidatorMap tempranges;
	tempranges[std::pair<double,double>(100,200)] = lowTempCheeseValidator;
	tempranges[std::pair<double,double>(200,300)] = highTempCheeseValidator;
	TEST_THROW(
		RCP<RangeValidatorDependency<double> > 
		cheeseTempDep = RCP<RangeValidatorDependency<double> >(
			new RangeValidatorDependency<double>(
			  list1->getEntryRCP("string parameter"),
				list1->getEntryRCP("Cheese to Fondue"), 
				tempranges
			)
		),
		InvalidDependencyException
	);

	tempranges[std::pair<double,double>(400,800)] = intVali;
	TEST_THROW(
		RCP<RangeValidatorDependency<double> > 
		cheeseTempDep = RCP<RangeValidatorDependency<double> >(
			new RangeValidatorDependency<double>(
			  list1->getEntryRCP("int parameter"),
				list1->getEntryRCP("Cheese to Fondue"), 
				tempranges
			)
		),
		InvalidDependencyException
	);

  RangeValidatorDependency<double>::RangeToValidatorMap emptyMap2;
	TEST_THROW(
		RCP<RangeValidatorDependency<double> > 
		emptyMapDep = RCP<RangeValidatorDependency<double> >(
			new RangeValidatorDependency<double>(
			  list1->getEntryRCP("double parameter"),
				list1->getEntryRCP("Cheese to Fondue"), 
				emptyMap2
			)
		),
		InvalidDependencyException
	);

	RangeValidatorDependency<int>::RangeToValidatorMap badRanges;
	tempranges[std::pair<int,int>(200,100)] = lowTempCheeseValidator;
	TEST_THROW(
		RCP<RangeValidatorDependency<int> > 
		cheeseTempDep = RCP<RangeValidatorDependency<int> >(
			new RangeValidatorDependency<int>(
			  list1->getEntryRCP("string parameter"),
				list1->getEntryRCP("Cheese to Fondue"), 
				badRanges
			)
		),
		InvalidDependencyException
	);
}

/**
 * Tests various DependencySheet functions.
 */
TEUCHOS_UNIT_TEST(Teuchos_Dependencies, DepSheetTest){
	RCP<DependencySheet> depSheet1 = rcp(new DependencySheet);
  TEST_ASSERT(depSheet1->empty());
  depSheet1->addDependency(
    DummyObjectGetter<BoolVisualDependency>::getDummyObject());
  TEST_ASSERT(!depSheet1->empty());


}


} //namespace Teuchos

