// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"


// 2010/07/30: rabartl: Here I just added all the unit tests to the Teuchos
// namespace to remove some clutter.


namespace Teuchos {


/**
 * Tests Number Validators.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, numberValidators)
{
	/*
	 * Testing Int Validator.
	 */
	RCP<ParameterList> intList =
    rcp(new ParameterList("Int List"));
	RCP<EnhancedNumberValidator<int> > intVali =
    rcp(new EnhancedNumberValidator<int>(0,10,4));
	TEST_ASSERT(intVali->getMin() == 0);
	TEST_ASSERT(intVali->getMax() == 10);
	TEST_ASSERT(intVali->getStep() == 4);
	TEST_ASSERT(intVali->hasMin());
	TEST_ASSERT(intVali->hasMax());
	RCP<EnhancedNumberValidator<int> > intVali2 =
    rcp(new EnhancedNumberValidator<int>());
	TEST_ASSERT(!intVali2->hasMin());
	TEST_ASSERT(!intVali2->hasMax());
	TEST_ASSERT(intVali2->getMin() == std::numeric_limits<int>::min());
	TEST_ASSERT(intVali2->getMax() == std::numeric_limits<int>::max());
	TEST_ASSERT(intVali2->getStep() == EnhancedNumberTraits<int>::defaultStep());
	intList->set("Int Parameter", 5, "int parameter", intVali);
	TEST_NOTHROW(intList->validateParameters(*intList));
	TEST_THROW(intList->set("Int Parameter", 11),
    Exceptions::InvalidParameterValue);
	TEST_THROW(intList->set("Double Parameter", 5.0, "double parameter", intVali),
    Exceptions::InvalidParameterType);

  // Test String Conversions with int
  RCP<ParameterList> validList = rcp(new ParameterList("Valid List"));
  RCP<ParameterList> userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(validList->set("Int Parameter", 4, "int parameter",
    intVali));
#ifdef HAVE_TEUCHOSCORE_CXX11
  TEST_NOTHROW(userList->set("Int Parameter", "x4"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
  TEST_NOTHROW(userList->set("Int Parameter", "4x"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
  TEST_NOTHROW(userList->set("Int Parameter", "12")); // ok string bad range
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
#endif
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(userList->set("Int Parameter", 4));
  TEST_NOTHROW(userList->validateParametersAndSetDefaults(*validList));
  TEST_NOTHROW(userList->set("Int Parameter", "8"));
  TEST_NOTHROW(userList->validateParametersAndSetDefaults(*validList));
  int readInt = userList->getEntry("Int Parameter").getValue<int>(&readInt);
  TEST_ASSERT(readInt == 8);

  // check string can generate out of range
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(userList->set("Int Parameter", "20"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Short Validator.
	 */
	RCP<ParameterList> shortList =
    rcp(new ParameterList("Short List"));
	RCP<EnhancedNumberValidator<short> > shortVali =
    rcp(new EnhancedNumberValidator<short>(0,10,4));
	TEST_ASSERT(shortVali->getMin() == 0);
	TEST_ASSERT(shortVali->getMax() == 10);
	TEST_ASSERT(shortVali->getStep() == 4);
	TEST_ASSERT(shortVali->hasMin());
	TEST_ASSERT(shortVali->hasMax());
	RCP<EnhancedNumberValidator<short> > shortVali2 =
    rcp(new EnhancedNumberValidator<short>());
	TEST_ASSERT(!shortVali2->hasMin());
	TEST_ASSERT(!shortVali2->hasMax());
	TEST_ASSERT(shortVali2->getMin() == std::numeric_limits<short>::min());
	TEST_ASSERT(shortVali2->getMax() == std::numeric_limits<short>::max());
	TEST_ASSERT(shortVali2->getStep() == EnhancedNumberTraits<short>::defaultStep());
	shortList->set("Short Parameter", (short)5, "short parameter", shortVali);
	TEST_NOTHROW(shortList->validateParameters(*shortList));
	TEST_THROW(shortList->set("Short Parameter", (short)11),
    Exceptions::InvalidParameterValue);
	TEST_THROW(shortList->set("Double Parameter", 5.0, "double parameter", shortVali),
    Exceptions::InvalidParameterType);

  // Test String Conversions with short
  validList = rcp(new ParameterList("Valid List"));
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(validList->set("Short Parameter", (short)4, "short parameter",
    shortVali));
#ifdef HAVE_TEUCHOSCORE_CXX11
  TEST_NOTHROW(userList->set("Short Parameter", "x4"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
  TEST_NOTHROW(userList->set("Short Parameter", "4x"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
  TEST_NOTHROW(userList->set("Short Parameter", "12")); // ok string bad range
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
#endif
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(userList->set("Short Parameter", (short)4));
  TEST_NOTHROW(userList->validateParametersAndSetDefaults(*validList));
  TEST_NOTHROW(userList->set("Short Parameter", "8"));
  TEST_NOTHROW(userList->validateParametersAndSetDefaults(*validList));
  short readShort = userList->getEntry("Short Parameter").getValue<short>(&readShort);
  TEST_ASSERT(readShort == 8);

  // check string can generate out of range
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(userList->set("Short Parameter", "20"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Float Validator.
	 */
	RCP<ParameterList> floatList = rcp(new ParameterList("Float List"));
	RCP<EnhancedNumberValidator<float> > floatVali =
    rcp(new EnhancedNumberValidator<float>(0,10.0,4.0,6));
	TEST_ASSERT(floatVali->getMin() == 0.0);
	TEST_ASSERT(floatVali->getMax() == 10.0);
	TEST_ASSERT(floatVali->getStep() == 4.0);
	TEST_ASSERT(floatVali->getPrecision() == 6);
	TEST_ASSERT(floatVali->hasMin());
	TEST_ASSERT(floatVali->hasMax());
	RCP<EnhancedNumberValidator<float> > floatVali2 =
    rcp(new EnhancedNumberValidator<float>());
	TEST_ASSERT(!floatVali2->hasMin());
	TEST_ASSERT(!floatVali2->hasMax());
	TEST_ASSERT(floatVali2->getMin() == EnhancedNumberTraits<float>::min());
	TEST_ASSERT(floatVali2->getMax() == EnhancedNumberTraits<float>::max());
	TEST_ASSERT(floatVali2->getStep() == EnhancedNumberTraits<float>::defaultStep());
	TEST_ASSERT(floatVali2->getPrecision() == EnhancedNumberTraits<float>::defaultPrecision());
	floatList->set("Float Parameter", (float)5.0, "float parameter", floatVali);
	TEST_NOTHROW(floatList->validateParameters(*floatList));
	TEST_THROW(floatList->set("Float Parameter", (float)11.0),
    Exceptions::InvalidParameterValue);
	TEST_THROW(floatList->set("Int Parameter", 5, "int parameter", floatVali),
    Exceptions::InvalidParameterType);

  // Test String Conversions with float
  validList = rcp(new ParameterList("Valid List"));
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(validList->set("Float Parameter", (float)4.0, "float parameter",
    floatVali));
#ifdef HAVE_TEUCHOSCORE_CXX11
  TEST_NOTHROW(userList->set("Float Parameter", "x4.0"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
  TEST_NOTHROW(userList->set("Float Parameter", "4.0x"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
  TEST_NOTHROW(userList->set("Float Parameter", "12.0")); // ok string bad range
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
#endif
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(userList->set("Float Parameter", (float)8.0));
  TEST_NOTHROW(userList->validateParametersAndSetDefaults(*validList));
  TEST_NOTHROW(userList->set("Float Parameter", "8.0"));
  TEST_NOTHROW(userList->validateParametersAndSetDefaults(*validList));
  float readFloat = userList->getEntry("Float Parameter").getValue<float>(&readFloat);
  TEST_ASSERT(readFloat == 8.0);

  // check string can generate out of range
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(userList->set("Float Parameter", "20.0"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Double Validator.
	 */
	RCP<ParameterList> doubleList = rcp(new ParameterList("Double List"));
	RCP<EnhancedNumberValidator<double> > doubleVali =
    rcp(new EnhancedNumberValidator<double>(0,10.0,4.0,6));
	TEST_ASSERT(doubleVali->getMin() == 0.0);
	TEST_ASSERT(doubleVali->getMax() == 10.0);
	TEST_ASSERT(doubleVali->getStep() == 4.0);
	TEST_ASSERT(doubleVali->getPrecision() == 6);
	TEST_ASSERT(doubleVali->hasMin());
	TEST_ASSERT(doubleVali->hasMax());
	RCP<EnhancedNumberValidator<double> > doubleVali2 =
    rcp(new EnhancedNumberValidator<double>());
	TEST_ASSERT(!doubleVali2->hasMin());
	TEST_ASSERT(!doubleVali2->hasMax());
	TEST_ASSERT(doubleVali2->getMin() == EnhancedNumberTraits<double>::min());
	TEST_ASSERT(doubleVali2->getMax() == EnhancedNumberTraits<double>::max());
	TEST_ASSERT(doubleVali2->getStep() == EnhancedNumberTraits<double>::defaultStep());
	TEST_ASSERT(doubleVali2->getPrecision() == EnhancedNumberTraits<double>::defaultPrecision());
	doubleList->set("Double Parameter", (double)5.0, "double parameter", doubleVali);
	TEST_NOTHROW(doubleList->validateParameters(*doubleList));
	TEST_THROW(doubleList->set("Double Parameter", (double)11.0),
    Exceptions::InvalidParameterValue);
	TEST_THROW(doubleList->set("Int Parameter", 5, "int parameter", doubleVali),
    Exceptions::InvalidParameterType);

  // Test String Conversions with double
  validList = rcp(new ParameterList("Valid List"));
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(validList->set("Double Parameter", 4.0, "double parameter",
    doubleVali));
#ifdef HAVE_TEUCHOSCORE_CXX11
  TEST_NOTHROW(userList->set("Double Parameter", "x4.0"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
  TEST_NOTHROW(userList->set("Double Parameter", "4.0x"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidArgument);
  TEST_NOTHROW(userList->set("Double Parameter", "12.0"));
  TEST_THROW(userList->validateParameters(*validList), // bad range
    Exceptions::InvalidArgument);
#endif
  userList = rcp(new ParameterList("Valid List"));
  TEST_NOTHROW(userList->set("Double Parameter", 8.0));
  TEST_NOTHROW(userList->validateParametersAndSetDefaults(*validList));
  TEST_NOTHROW(userList->set("Double Parameter", "8.0"));
  TEST_NOTHROW(userList->validateParametersAndSetDefaults(*validList));
  double readDouble = userList->getEntry("Double Parameter").getValue<double>(&readDouble);
  TEST_ASSERT(readDouble == 8.0);

  // check string can generate out of range
  userList = rcp(new ParameterList("User List"));
  TEST_NOTHROW(userList->set("Double Parameter", "20.0"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidParameterType);
 }

TEUCHOS_UNIT_TEST(Teuchos_Validators, anyNumberValidator)
{
	RCP<ParameterList> userList = rcp(new ParameterList("User List"));
	RCP<ParameterList> validList = rcp(new ParameterList("Valid List"));

  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes intDoubleTypes;
  intDoubleTypes.allowString(false);
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes intStringTypes;
  intStringTypes.allowDouble(false);
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes intTypes;
  intTypes.allowDouble(false);
  intTypes.allowString(false);

  // set up validators to test
  // default prefers double and allows string and int
  RCP<Teuchos::AnyNumberParameterEntryValidator> allValidator =
    Teuchos::rcp( new Teuchos::AnyNumberParameterEntryValidator() );
  RCP<Teuchos::AnyNumberParameterEntryValidator> intDoubleValidator =
    Teuchos::rcp( new Teuchos::AnyNumberParameterEntryValidator(
    Teuchos::AnyNumberParameterEntryValidator::PREFER_DOUBLE, intDoubleTypes) );
  RCP<Teuchos::AnyNumberParameterEntryValidator> intStringValidator =
    Teuchos::rcp( new Teuchos::AnyNumberParameterEntryValidator(
    Teuchos::AnyNumberParameterEntryValidator::PREFER_INT, intStringTypes) );
  RCP<Teuchos::AnyNumberParameterEntryValidator> intValidator =
    Teuchos::rcp( new Teuchos::AnyNumberParameterEntryValidator(
    Teuchos::AnyNumberParameterEntryValidator::PREFER_INT, intTypes) );

  // first check the 'good' setups which do not throw
  TEST_NOTHROW(validList->set( "allParameter", "1.1", "documentation",
    allValidator));
  TEST_NOTHROW(validList->set( "allParameter", 1.1, "documentation",
    allValidator));
  TEST_NOTHROW(validList->set( "allParameter", "1", "documentation",
    allValidator));
  TEST_NOTHROW(validList->set( "allParameter", 1, "documentation",
    allValidator));
  TEST_NOTHROW(validList->set( "intDoubleParameter", 1.1, "documentation",
    intDoubleValidator));
  TEST_NOTHROW(validList->set( "intDoubleParameter", 1, "documentation",
    intDoubleValidator));
  TEST_NOTHROW(validList->set( "intStringParameter", "1", "documentation",
    intStringValidator));
  TEST_NOTHROW(validList->set( "intStringParameter", 1, "documentation",
    intStringValidator));
  TEST_NOTHROW(validList->set( "intParameter", 1, "documentation",
    intValidator));

  // This was a special case that might warrant discussion.
  // The issue is for validators which accept string/int but not double.
  // In the original setup the validator would always call getDouble
  // internally and accept a string of "1.1" without error.
  TEST_NOTHROW(validList->set( "intStringParameter", "1.1", "documentation",
    intStringValidator));

  //
  // these are some cases which  throw independent of HAVE_TEUCHOSCORE_CXX11
  //

  // if string it not allowed you can't use a string ever
  TEST_THROW(validList->set( "intDoubleParameter", "1.1", "documentation",
    intDoubleValidator), Exceptions::InvalidParameterType);

   // it also throws for a double number - double not allowed
  TEST_THROW(validList->set( "intStringParameter", 1.1, "documentation",
    intStringValidator), Exceptions::InvalidArgument);

  // for int only it can't be a string - any string will throw
  TEST_THROW(validList->set( "intParameter", "1", "documentation",
    intValidator), Exceptions::InvalidParameter);

  // this int only it can't be a double because double is not allowed
  TEST_THROW(validList->set( "intParameter", 1.1, "documentation",
    intValidator), Exceptions::InvalidParameter);

  //
  // these behaviors now depend on HAVE_TEUCHOSCORE_CXX11
  // std::stod and std::stoi will be used for HAVE_TEUCHOSCORE_CXX11
  // std::atof and std::atoi will be used for no CXX11
  //
#ifdef HAVE_TEUCHOSCORE_CXX11
  // for double types we throw for badly formatted string on std::stod
  // this will check the double type first because it is PREFER_DOUBLE
  TEST_THROW(validList->set( "allParameter", "1.1x", "documentation",
    allValidator), Exceptions::InvalidArgument);
  TEST_THROW(validList->set( "intDoubleParameter", "1.1x", "documentation",
    allValidator), Exceptions::InvalidArgument);
  TEST_THROW(validList->set( "allParameter", "x1.1", "documentation",
    allValidator), Exceptions::InvalidArgument);
  TEST_THROW(validList->set( "intDoubleParameter", "x1.1", "documentation",
    allValidator), Exceptions::InvalidArgument);
  // for int/string but no double - std::stoi throws for invalid formatting
  TEST_THROW(validList->set( "intStringParameter", "1x", "documentation",
    intStringValidator), Exceptions::InvalidArgument);
  TEST_THROW(validList->set( "intStringParameter", "x1", "documentation",
    intStringValidator), Exceptions::InvalidArgument);
  TEST_THROW(validList->set( "intStringParameter", "1 x", "documentation",
    intStringValidator), Exceptions::InvalidArgument);
#else
  // for int/double/string std::atod does NOT throw - this is the old behavior
  // this is different now when HAVE_TEUCHOSCORE_CXX11 is ON - see above
  TEST_NOTHROW(validList->set( "allParameter", "1.1x", "documentation",
    allValidator));
  // for int/string std::atoi does NOT throw - this is the old behavior
  // this is different now when HAVE_TEUCHOSCORE_CXX11 is ON - see above
  TEST_NOTHROW(validList->set( "intStringParameter", "1.1x", "documentation",
    intStringValidator));
#endif
}

TEUCHOS_UNIT_TEST(Teuchos_Validators, boolValidator)
{
  RCP<ParameterList> userList = rcp(new ParameterList("User List"));
  RCP<ParameterList> validList = rcp(new ParameterList("Valid List"));

  // first without validator - accepts only true/false
  validList->set( "justBool", false, "documentation" );
  TEST_NOTHROW(userList->set( "justBool", false));
  TEST_NOTHROW(userList->validateParameters(*validList));
  TEST_NOTHROW(userList->set( "justBool", true));
  TEST_NOTHROW(userList->validateParameters(*validList));
  // this will not validate because we did not add a bool validator
  TEST_NOTHROW(userList->set( "justBool", "true"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidParameterType);
  // this will not validate because we did not add a bool validator
  TEST_NOTHROW(userList->set( "justBool", "false"));
  TEST_THROW(userList->validateParameters(*validList),
    Exceptions::InvalidParameterType);

  // now with BoolParameterEntryValidator validator
  // accepts true/false/"true"/"false"
  RCP<Teuchos::BoolParameterEntryValidator> boolValidator =
    Teuchos::rcp( new Teuchos::BoolParameterEntryValidator() );
  userList = rcp(new ParameterList("User List")); // make a new list
  validList = rcp(new ParameterList("Valid List")); // make a new list
  validList->set( "boolOrString", false, "documentation", boolValidator );
  TEST_NOTHROW(userList->set( "boolOrString", false));
  TEST_NOTHROW(userList->validateParameters(*validList));
  TEST_NOTHROW(userList->set( "boolOrString", true));
  TEST_NOTHROW(userList->validateParameters(*validList));
  // this will validate because we added a bool validator
  TEST_NOTHROW(userList->set( "boolOrString", "true"));
  TEST_NOTHROW(userList->validateParameters(*validList));
  // this will validate because we added a bool validator
  TEST_NOTHROW(userList->set( "boolOrString", "false"));
  TEST_NOTHROW(userList->validateParameters(*validList));
  // but only "false" and "true" work - anything else will not validate
  TEST_NOTHROW(userList->set( "boolOrString", "falsex")); // sets ok
  TEST_THROW(userList->validateParameters(*validList), // but throws
    Exceptions::InvalidParameterType);

  // now with BoolParameterEntryValidator validator
  // but consider what happens if we created it using "false" instead of false
  // this should still work identically to the previous case
  userList = rcp(new ParameterList("User List")); // make a new list
  validList = rcp(new ParameterList("Valid List")); // make a new list
  validList->set( "boolOrString", "false", "documentation", boolValidator );
  TEST_NOTHROW(userList->set( "boolOrString", false));
  TEST_NOTHROW(userList->validateParameters(*validList));
  TEST_NOTHROW(userList->set( "boolOrString", true ));
  TEST_NOTHROW(userList->validateParameters(*validList));
  // this will validate because we added a bool validator
  TEST_NOTHROW(userList->set( "boolOrString", "true"));
  TEST_NOTHROW(userList->validateParameters(*validList));
  // this will validate because we added a bool validator
  TEST_NOTHROW(userList->set( "boolOrString", "false"));
  TEST_NOTHROW(userList->validateParameters(*validList));
  // but only "false" and "true" work - anything else will not validate
  TEST_NOTHROW(userList->set( "boolOrString", "falsex")); // sets ok
  TEST_THROW(userList->validateParameters(*validList), // but throws
    Exceptions::InvalidParameterType);

  // do another test using validateParametersAndSetDefaults
  userList = rcp(new ParameterList("User List")); // make a new list
  validList = rcp(new ParameterList("Valid List")); // make a new list
  // Default values for parameters are bool
  validList->set("boolOne", true, "doc", boolValidator);
  validList->set("boolTwo", false, "doc", boolValidator);
  bool defOne = validList->getEntry("boolOne").getValue(&defOne);
  bool defTwo = validList->getEntry("boolTwo").getValue(&defTwo);

  // Create user parameter list
  userList->set("boolOne", false);   // User can provide bool value...
  userList->set("boolTwo", "true");  // or string "true"/"false"
  TEST_NOTHROW(userList->validateParametersAndSetDefaults(*validList));
}


/*
 * Testing StringValidator.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, stringValidator)
{
	RCP<ParameterList> stringList = rcp(new ParameterList("String List"));
	Array<std::string> stringVals = tuple<std::string>("str1", "str2", "str3");
	RCP<StringValidator> stringVali = rcp(new StringValidator(stringVals));
	RCP<const Array<std::string> > valiVals = stringVali->validStringValues();
  /*bool local_success = true;
  for(int i =0; i<valiVals.size() ++i){
	  TEST_ARRAY_ELE_EQUALITY(*valiVals, i, stringVals[i]);
  }
  if (local_success) out << "passed\n";
  else success = false;*/
  TEST_COMPARE_ARRAYS(*valiVals, stringVals);
	TEST_NOTHROW(stringList->set("String param1", "str1", "a string parameter", stringVali));
	TEST_THROW(stringList->set("String param2", "not in list", "a string parameter", stringVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("int param", 5, "a int parameter", stringVali),
    Exceptions::InvalidParameterType);
}


/*
 * Testing StringToIntegralParameterEntryValidator.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, StringToIntegralParameterEntryValidator) {
  Array<std::string> strVals = tuple<std::string>("str1", "str2", "str3");
  Array<std::string> strDocs = tuple<std::string>("a str1", "a str2", "a str3");
  Array<int> intVals = tuple<int>(1, 2, 3);
  bool caseSensitive = true;
  typedef StringToIntegralParameterEntryValidator<int> ret_type;
  // Note that validator1 maps the strings to {0, 1, 2} not {1, 2, 3} as in `intVals`
  RCP<ret_type> validator1 = rcp(new ret_type(strVals, "str1", caseSensitive));
  RCP<ret_type> validator2 = rcp(new ret_type(strVals, intVals, "str1", caseSensitive));
  RCP<ret_type> validator3 = rcp(new ret_type(strVals, strDocs, intVals, "str1", caseSensitive));
  TEST_EQUALITY(strDocs, *validator3->getStringDocs());
  ParameterList valid_pl = ParameterList();
  valid_pl.set("Param1", "str1", "Parameter 1", validator1);
  valid_pl.set("Param2", "str1", "Parameter 2", validator2);
  valid_pl.set("Param3", "str1", "Parameter 3", validator3);
  ParameterList user_pl = ParameterList();
  user_pl.set("Param1", "str1");
  user_pl.set("Param2", "str2");
  user_pl.set("Param3", "str3");
  // Test `getStringValue` and `getIntegralValue` before validation on `valid_pl`
  TEST_EQUALITY(0,          getIntegralValue<int>(valid_pl, "Param1"));
  TEST_EQUALITY(intVals[0], getIntegralValue<int>(valid_pl, "Param2"));
  TEST_EQUALITY(strVals[0], getStringValue<int>(valid_pl, "Param2"));
  // Test `getStringValue` and `getIntegralValue` after validation on `user_pl`
  user_pl.validateParametersAndSetDefaults(valid_pl);
  TEST_EQUALITY(0,          getIntegralValue<int>(user_pl, "Param1"));
  TEST_EQUALITY(intVals[1], getIntegralValue<int>(user_pl, "Param2"));
  TEST_EQUALITY(intVals[2], getIntegralValue<int>(user_pl, "Param3"));
  TEST_EQUALITY(strVals[0], getStringValue<int>(user_pl, "Param1"));
  TEST_EQUALITY(strVals[1], getStringValue<int>(user_pl, "Param2"));
  TEST_EQUALITY(strVals[2], getStringValue<int>(user_pl, "Param3"));
}


/*
 * Testing FileNameValidator.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, fileNameValidator)
{
	RCP<ParameterList> fileNameList = rcp(new ParameterList("Filename List"));
	RCP<FileNameValidator> fileNameVali = rcp(new FileNameValidator(true));
	TEST_ASSERT(fileNameVali->fileMustExist());
	fileNameVali->setFileMustExist(false);
	TEST_ASSERT(!fileNameVali->fileMustExist());
	TEST_NOTHROW(fileNameList->set("File name param", "../path", "file name parameter",
      fileNameVali));
	TEST_THROW(fileNameList->set("int param", 5, "int parameter", fileNameVali),
    Exceptions::InvalidParameterType);
	fileNameVali->setFileMustExist(true);
	TEST_NOTHROW(fileNameList->set("file name param", "testFile.txt", "a file name", fileNameVali));
	TEST_THROW(fileNameList->set("file name param", "doesntexist.txt", "a file name", fileNameVali),
    Exceptions::InvalidParameterValue);
}


/*
 * Testing Array Validators.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, arrayValidators)
{

	/*
	 * Testing StringArrayValidator.
	 */
	RCP<ParameterList> stringList = rcp(new ParameterList("String List"));
	Array<std::string> stringVals = tuple<std::string>("str1", "str2", "str3");
	RCP<StringValidator> stringVali = rcp(new StringValidator(stringVals));
	RCP<ArrayStringValidator> stringArrayVali = rcp(new ArrayStringValidator(stringVali));
	TEST_ASSERT(stringVali.get() == stringArrayVali->getPrototype().get());
	Array<std::string> stringArray = tuple<std::string>("str2","str3","str1","str3","str2");
	TEST_NOTHROW(stringList->set("String Array Param", stringArray, "string array parameter", stringArrayVali));
	Array<std::string> badStringArray = tuple<std::string>("not valid","str3","str1","str3","str2");
	TEST_THROW(stringList->set("String Array Param", badStringArray, "string array parameter", stringArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Int param", 5, "int parameter", stringArrayVali),
    Exceptions::InvalidParameterType);
	Array<long> longArray = tuple<long>((long)5,(long)5,(long)3);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", stringArrayVali),
    Exceptions::InvalidParameterType);
	
	/*
	 * Testing Int ArrayValidator.
	 */
	RCP<ParameterList> intList = rcp(new ParameterList("Int List"));
	RCP<EnhancedNumberValidator<int> > intVali = rcp(new EnhancedNumberValidator<int>(0, 10));
	RCP<ArrayNumberValidator<int> > intArrayVali = rcp(new ArrayNumberValidator<int>(intVali));
	TEST_ASSERT(intVali.get() == intArrayVali->getPrototype().get());
	Array<int> intArray = tuple<int>(1,4,2,5);
	TEST_NOTHROW(intList->set("int array param", intArray, "int array parameter", intArrayVali));
	Array<int> intBadArray = tuple<int>(11,4,2,5);
	TEST_THROW(intList->set("int bad array param", intBadArray, "int bad array parameter", intArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", intArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Short ArrayValidator.
	 */
	RCP<ParameterList> shortList = rcp(new ParameterList("Short List"));
	RCP<EnhancedNumberValidator<short> > shortVali =
    rcp(new EnhancedNumberValidator<short>(0, 10));
	RCP<ArrayNumberValidator<short> > shortArrayVali =
    rcp(new ArrayNumberValidator<short>(shortVali));
	TEST_ASSERT(shortVali.get() == shortArrayVali->getPrototype().get());
	Array<short> shortArray = tuple<short>(1,4,2,5);
	TEST_NOTHROW(shortList->set("short array param", shortArray, "short array parameter", shortArrayVali));
	Array<short> shortBadArray = tuple<short>(11,4,2,5);
	TEST_THROW(shortList->set("short bad array param", shortBadArray, "short bad array parameter", shortArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", shortArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Float ArrayValidator.
	 */
	RCP<ParameterList> floatList = rcp(new ParameterList("Float List"));
	RCP<EnhancedNumberValidator<float> > floatVali =
    rcp(new EnhancedNumberValidator<float>(0.0, 10.0));
	RCP<ArrayNumberValidator<float> > floatArrayVali =
    rcp(new ArrayNumberValidator<float>(floatVali));
	TEST_ASSERT(floatVali.get() == floatArrayVali->getPrototype().get());
	Array<float> floatArray = tuple<float>(1.0,4.0,2.0,5.0);
	TEST_NOTHROW(floatList->set("float array param", floatArray, "float array parameter", floatArrayVali));
	Array<float> floatBadArray = tuple<float>(11.0,4.0,2.0,5.0);
	TEST_THROW(floatList->set("float bad array param", floatBadArray, "float bad array parameter", floatArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", floatArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Double ArrayValidator.
	 */
	RCP<ParameterList> doubleList = rcp(new ParameterList("Double List"));
	RCP<EnhancedNumberValidator<double> > doubleVali =
    rcp(new EnhancedNumberValidator<double>(0.0, 10.0));
	RCP<ArrayNumberValidator<double> > doubleArrayVali =
    rcp(new ArrayNumberValidator<double>(doubleVali));
	TEST_ASSERT(doubleVali.get() == doubleArrayVali->getPrototype().get());
	Array<double> doubleArray = tuple<double>(1.0,4.0,2.0,5.0);
	TEST_NOTHROW(doubleList->set("double array param", doubleArray, "double array parameter", doubleArrayVali));
	Array<double> doubleBadArray = tuple<double>(11.0,4.0,2.0,5.0);
	TEST_THROW(doubleList->set("double bad array param", doubleBadArray, "double bad array parameter", doubleArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", doubleArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing FileName ArrayValidator.
	 */
	RCP<ParameterList> fileNameList = rcp(new ParameterList("Filename List"));
	RCP<FileNameValidator> fileNameVali = rcp(new FileNameValidator(true));
	RCP<ArrayFileNameValidator> arrayFileNameVali = rcp(new ArrayFileNameValidator(fileNameVali));
	TEST_ASSERT(arrayFileNameVali->getPrototype().get() == fileNameVali.get());
	Array<std::string> fileNameArray = tuple<std::string>("testFile.txt", "testFile2.txt", "testFile3.txt");
	Array<std::string> fileNameBadArray = tuple<std::string>("doesnexist.txt", "testFile2.txt", "testFile3.txt");
	TEST_NOTHROW(fileNameList->set("File name array", fileNameArray, "file name array parameter", arrayFileNameVali));
	TEST_THROW(fileNameList->set("Bad File name array", fileNameBadArray, "bad file name array parameter", arrayFileNameVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", arrayFileNameVali),
    Exceptions::InvalidParameterType);
}

/*
 * Testing TwoDArray Validators.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, twoDArrayValidators)
{

	/*
	 * Testing StringArrayValidator.
	 */
	RCP<ParameterList> stringList = rcp(new ParameterList("String List"));
	Array<std::string> stringVals = tuple<std::string>("str1", "str2", "str3");
	RCP<StringValidator> stringVali = rcp(new StringValidator(stringVals));
	RCP<TwoDArrayStringValidator> stringArrayVali =
    rcp(new TwoDArrayStringValidator(stringVali));
	TEST_ASSERT(stringVali.get() == stringArrayVali->getPrototype().get());
	TwoDArray<std::string> stringArray(2,2);
  stringArray(0,0) = "str2";
  stringArray(0,1) = "str1";
  stringArray(1,0) = "str3";
  stringArray(1,1) = "str2";
	TEST_NOTHROW(stringList->set("String Array Param", stringArray, "string array parameter", stringArrayVali));
	TwoDArray<std::string> badStringArray(2,2);
  badStringArray(0,0) = "str2";
  badStringArray(0,1) = "str1";
  badStringArray(1,0) = "str3";
  badStringArray(1,1) = "not valid";
	TEST_THROW(stringList->set("String Array Param", badStringArray, "string array parameter", stringArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Int param", 5, "int parameter", stringArrayVali),
    Exceptions::InvalidParameterType);
	TwoDArray<long> longArray(2,2);
  longArray(0,0) = (long)5;
  longArray(0,1) = (long)4;
  longArray(1,0) = (long)9;
  longArray(1,1) = (long)1;
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", stringArrayVali),
    Exceptions::InvalidParameterType);
	
	/*
	 * Testing Int ArrayValidator.
	 */
	RCP<ParameterList> intList = rcp(new ParameterList("Int List"));
	RCP<EnhancedNumberValidator<int> > intVali = rcp(new EnhancedNumberValidator<int>(0, 10));
	RCP<TwoDArrayNumberValidator<int> > intArrayVali =
    rcp(new TwoDArrayNumberValidator<int>(intVali));
	TEST_ASSERT(intVali.get() == intArrayVali->getPrototype().get());
	TwoDArray<int> intArray(2,2);
  intArray(0,0) = 1;
  intArray(0,1) = 4;
  intArray(1,0) = 2;
  intArray(1,1) = 5;
	TEST_NOTHROW(intList->set("int array param", intArray, "int array parameter", intArrayVali));
	TwoDArray<int> intBadArray(2,2);
  intBadArray(0,0) = 11;
  intBadArray(0,1) = 4;
  intBadArray(1,0) = 2;
  intBadArray(1,1) = 5;
	TEST_THROW(intList->set("int bad array param", intBadArray, "int bad array parameter", intArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", intArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Short ArrayValidator.
	 */
	RCP<ParameterList> shortList = rcp(new ParameterList("Short List"));
	RCP<EnhancedNumberValidator<short> > shortVali =
    rcp(new EnhancedNumberValidator<short>(0, 10));
	RCP<TwoDArrayNumberValidator<short> > shortArrayVali =
    rcp(new TwoDArrayNumberValidator<short>(shortVali));
	TEST_ASSERT(shortVali.get() == shortArrayVali->getPrototype().get());
	TwoDArray<short> shortArray(2,2);
  shortArray(0,0) = 1;
  shortArray(0,1) = 4;
  shortArray(1,0) = 2;
  shortArray(1,1) = 5;
	TEST_NOTHROW(shortList->set("short array param", shortArray, "short array parameter", shortArrayVali));
	TwoDArray<short> shortBadArray(2,2);
  shortBadArray(0,0) = 11;
  shortBadArray(0,1) = 4;
  shortBadArray(1,0) = 2;
  shortBadArray(1,1) = 5;
	TEST_THROW(shortList->set("short bad array param", shortBadArray, "short bad array parameter", shortArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", shortArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Float ArrayValidator.
	 */
	RCP<ParameterList> floatList = rcp(new ParameterList("Float List"));
	RCP<EnhancedNumberValidator<float> > floatVali =
    rcp(new EnhancedNumberValidator<float>(0.0, 10.0));
	RCP<TwoDArrayNumberValidator<float> > floatArrayVali =
    rcp(new TwoDArrayNumberValidator<float>(floatVali));
	TEST_ASSERT(floatVali.get() == floatArrayVali->getPrototype().get());
	TwoDArray<float> floatArray(2,2);
  floatArray(0,0) = 1.0;
  floatArray(0,1) = 4.0;
  floatArray(1,0) = 5.0;
  floatArray(1,1) = 2.0;
	TEST_NOTHROW(floatList->set("float array param", floatArray, "float array parameter", floatArrayVali));
	TwoDArray<float> floatBadArray(2,2);
  floatBadArray(0,0) = 11.0;
  floatBadArray(0,1) = 4.0;
  floatBadArray(1,0) = 5.0;
  floatBadArray(1,1) = 2.0;
	TEST_THROW(floatList->set("float bad array param", floatBadArray, "float bad array parameter", floatArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", floatArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Double ArrayValidator.
	 */
	RCP<ParameterList> doubleList = rcp(new ParameterList("Double List"));
	RCP<EnhancedNumberValidator<double> > doubleVali =
    rcp(new EnhancedNumberValidator<double>(0.0, 10.0));
	RCP<TwoDArrayNumberValidator<double> > doubleArrayVali =
    rcp(new TwoDArrayNumberValidator<double>(doubleVali));
	TEST_ASSERT(doubleVali.get() == doubleArrayVali->getPrototype().get());
	TwoDArray<double> doubleArray(2,2);
  doubleArray(0,0) = 1.0;
  doubleArray(0,1) = 4.0;
  doubleArray(1,0) = 5.0;
  doubleArray(1,1) = 2.0;
	TEST_NOTHROW(doubleList->set("double array param", doubleArray, "double array parameter", doubleArrayVali));
	TwoDArray<double> doubleBadArray(2,2);
  doubleBadArray(0,0) = 11.0;
  doubleBadArray(0,1) = 4.0;
  doubleBadArray(1,0) = 5.0;
  doubleBadArray(1,1) = 2.0;
	TEST_THROW(doubleList->set("double bad array param", doubleBadArray, "double bad array parameter", doubleArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", doubleArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing FileName ArrayValidator.
	 */
	RCP<ParameterList> fileNameList = rcp(new ParameterList("Filename List"));
	RCP<FileNameValidator> fileNameVali = rcp(new FileNameValidator(true));
	RCP<TwoDArrayFileNameValidator> arrayFileNameVali =
    rcp(new TwoDArrayFileNameValidator(fileNameVali));
	TEST_ASSERT(arrayFileNameVali->getPrototype().get() == fileNameVali.get());
	TwoDArray<std::string> fileNameArray(2,2);
  fileNameArray(0,0) = "testFile.txt";
  fileNameArray(0,1) = "testFile2.txt";
  fileNameArray(1,0) = "testFile3.txt";
  fileNameArray(1,1) = "testFile.txt";
	TwoDArray<std::string> fileNameBadArray(2,2);
  fileNameBadArray(0,0) = "doesntexist.txt";
  fileNameBadArray(0,1) = "testFile2.txt";
  fileNameBadArray(1,0) = "testFile3.txt";
  fileNameBadArray(1,1) = "testFile.txt";
	TEST_NOTHROW(fileNameList->set("File name array", fileNameArray, "file name array parameter", arrayFileNameVali));
	TEST_THROW(fileNameList->set("Bad File name array", fileNameBadArray, "bad file name array parameter", arrayFileNameVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", arrayFileNameVali),
    Exceptions::InvalidParameterType);
}


} // namespace Teuchos

