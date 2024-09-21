// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ParameterListModifier.hpp"
#include "Teuchos_UnitTestHarness.hpp"


//
// Utilities
//


namespace {


class DummyValidator : public Teuchos::ParameterEntryValidator
{
public:

  const std::string getXMLTypeName() const { return ""; }
  virtual void printDoc(std::string const& docString, std::ostream &out) const {}
  virtual ValidStringsList validStringValues() const { return Teuchos::null; }
  virtual void validate(
    Teuchos::ParameterEntry const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const
    {}
};


class SimpleModifier : public Teuchos::ParameterListModifier
{
public:

  SimpleModifier() : Teuchos::ParameterListModifier("Simple Modifier"){}

  void modify(Teuchos::ParameterList &pl, Teuchos::ParameterList &valid_pl) const
  {
    expandSublistsUsingBaseName("SubA", pl, valid_pl);
  }

  void reconcile(Teuchos::ParameterList &pl) const
  {
    // If A and B are less than 0.0 then throw an error
    TEUCHOS_TEST_FOR_EXCEPTION(pl.get<double>("A") < 0.0 && pl.get<double>("B") < 0.0,
        std::logic_error, "Parameters A and B can't both be less than 0.0");
  }
};


class SimpleSubModifier : public Teuchos::ParameterListModifier {

public:

  SimpleSubModifier() : Teuchos::ParameterListModifier("Simple Sub Modifier"){}

  void modify(Teuchos::ParameterList &pl, Teuchos::ParameterList &valid_pl) const
  {
    expandSublistsUsingBaseName("SubB", pl, valid_pl);
  }
  void reconcile(Teuchos::ParameterList &pl) const
  {
    // If E and F are less than 10 then throw an error
    const int max_CD = 10;
    TEUCHOS_TEST_FOR_EXCEPTION(pl.get<int>("C") > max_CD && pl.get<int>("D") > max_CD,
        std::logic_error, "Parameters C and D can't both be greater than 10")
  }
};


class ReconciliationModifier1 : public Teuchos::ParameterListModifier
{
public:
  ReconciliationModifier1() : Teuchos::ParameterListModifier("Reconciliation Modifier 1"){}
  void reconcile(Teuchos::ParameterList &pl) const
  {
    // This reconciliation routine needs the ReconciliationModifier2's reconcile method
    // to be run first to create the "e" parameter.
    Teuchos::ParameterList &subA = pl.sublist("A");
    pl.set("b", subA.get<int>("e"));
  }
};


class ReconciliationModifier2 : public Teuchos::ParameterListModifier
{
public:
  ReconciliationModifier2() : Teuchos::ParameterListModifier("Reconciliation Modifier 2"){}
  void reconcile(Teuchos::ParameterList &pl) const
  {
    // Add a convenience parameter
    pl.set("e", pl.get<int>("c") + pl.get<int>("d"));
  }
};


} // namespace


namespace Teuchos {


//
// Test help utilities
//


ParameterList createMainPL()
{
  ParameterList PL_Main("PL_Main");
  const std::string Direction_Doc = "This sublist controls how direction is computed.";
  ParameterList &PL_Direction = PL_Main.sublist("Direction", false, Direction_Doc);
  ParameterList &PL_Newton = PL_Direction.sublist("Newton");
  PL_Newton.sublist("Linear Solver");
  PL_Main.sublist("Line Search");
  return PL_Main;
}



ParameterList createValidMainPL()
{

  ParameterList PL_Main_valid("PL_Main_valid");
  PL_Main_valid.setParameters(createMainPL());

  // Create a validator for the "Nonlinear Solver" parameter
  setStringToIntegralParameter<int>(
    "Nonlinear Solver",
    "Line Search Based",
    "Selects the type of nonlinear solver to use",
    tuple<std::string>("Line Search Based","Trust Region Based"),
    &PL_Main_valid
    );

  // Create a validator for the parameter "Line Search"->"Polynomial"->"Max Iters"
  // that accepts an 'int', a 'double' or a 'std::string' value!
  typedef Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes AcceptedTypes;
  Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
    linesearchMaxItersValidator = rcp(
      new Teuchos::AnyNumberParameterEntryValidator(
        Teuchos::AnyNumberParameterEntryValidator::PREFER_INT, // Not used here!
        AcceptedTypes(false).allowInt(true).allowDouble(true).allowString(true)
        )
      );
  PL_Main_valid.sublist("Line Search").sublist("Polynomial").set(
    "Max Iters",3
    ,"The maximum number of inner linear search iterations allowed."
    ,linesearchMaxItersValidator
    );

  // Create a validator for the parameter "Direction"->"Newton"->"Linear Solver"->"Tol"
  // that accepts a 'double' or a 'std::string' value!
  Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
    linSolveTolValidator = rcp(
      new Teuchos::AnyNumberParameterEntryValidator(
        Teuchos::AnyNumberParameterEntryValidator::PREFER_INT, // Not used here!
        AcceptedTypes(false).allowDouble(true).allowString(true)
        )
      );
  PL_Main_valid.sublist("Direction",true).sublist("Newton",true)
    .sublist("Linear Solver",true).set(
      "Tol", double(1e-5)
      ,"Select the linear solve tolerance"
    ,linSolveTolValidator
    );

  // Create a validator for the parameter "Elements"
  // that accepts an 'int', a 'long long' or a 'std::string' value!
  Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
    elementsValidator = rcp(
      new Teuchos::AnyNumberParameterEntryValidator(
        Teuchos::AnyNumberParameterEntryValidator::PREFER_LONG_LONG, // Not used here!
        AcceptedTypes(false).allowInt(true).allowLongLong(true).allowString(true)
        )
      );
  typedef long long LL;
  PL_Main_valid.set(
      "Elements", LL(72057594037927936ll) // 2^56
      ,"Number of finite elements to generate"
    ,elementsValidator
    );

  return PL_Main_valid;

}


//
// Unit tests
//


TEUCHOS_UNIT_TEST( ParameterList, construct_default )
{
  ParameterList pl;
  TEST_EQUALITY_CONST(pl.name(), "ANONYMOUS");
  TEST_EQUALITY_CONST(pl.numParams(), 0);
}


TEUCHOS_UNIT_TEST( ParameterList, construct_withName )
{
  ParameterList pl("someName");
  TEST_EQUALITY_CONST(pl.name(), "someName");
  TEST_EQUALITY_CONST(pl.numParams(), 0);
}


TEUCHOS_UNIT_TEST( ParameterList, createParameterList_empty )
{
  RCP<ParameterList> pl = createParameterList();
  TEST_ASSERT(nonnull(pl));
  TEST_EQUALITY_CONST(pl->name(), "ANONYMOUS");
}


TEUCHOS_UNIT_TEST( ParameterList, createParameterList_withName )
{
  RCP<ParameterList> pl = createParameterList("dummyName");
  TEST_ASSERT(nonnull(pl));
  TEST_EQUALITY_CONST(pl->name(), "dummyName");
}


TEUCHOS_UNIT_TEST( ParameterList, set_get_int )
{
  ParameterList pl;

  out << "\n";
  ECHO(pl.set("my int", 3));

  out << "\n";
  ECHO(const ParameterEntry& my_int_c_param = getConst(pl).getEntry("my int"));
  TEST_EQUALITY_CONST(my_int_c_param.isUsed(), false);
  TEST_EQUALITY_CONST(my_int_c_param.isList(), false);
  TEST_EQUALITY_CONST(my_int_c_param.isDefault(), false);
  TEST_EQUALITY_CONST(my_int_c_param.docString(), "");
  TEST_ASSERT(is_null(my_int_c_param.validator()));
  TEST_EQUALITY_CONST(getValue<int>(my_int_c_param), 3);
  ECHO(const bool param_isType_int1 = my_int_c_param.isType<int>());
  TEST_EQUALITY_CONST(param_isType_int1, true);
  ECHO(const bool param_isType_double1 = my_int_c_param.isType<double>());
  TEST_EQUALITY_CONST(param_isType_double1, false);

  out << "\n";
  ECHO(const ParameterEntry& my_int_param = pl.getEntry("my int"));
  TEST_EQUALITY_CONST(my_int_param.isUsed(), true);

  out << "\n";
  ECHO(const int my_int = pl.get<int>("my int"));
  TEST_EQUALITY_CONST(my_int, 3);

}


TEUCHOS_UNIT_TEST( ParameterList, setParametersWithModifier )
{
  using Teuchos::ParameterListModifier;
  RCP<ParameterListModifier> modifier1 = rcp(new ParameterListModifier("Modifier 1"));
  RCP<ParameterListModifier> modifier2 = rcp(new ParameterListModifier("Modifier 2"));
  //pl1:
  //  A: 1.0
  //  SubA: # with `modifier1`
  //    B: 2
  ParameterList pl1("pl");
  pl1.set("A", 1.0);
  pl1.sublist("SubA", modifier1).set("B", 2);
  ParameterList pl2("pl");
  pl2.setParameters(pl1);
  // Show that all values and the modifier were copied from `pl1` to `pl2`
  TEST_EQUALITY(*pl2.sublist("SubA").getModifier(), *modifier1);
  TEST_EQUALITY(pl1, pl2);
  pl1.sublist("SubA").setModifier(Teuchos::null);
  TEST_ASSERT(is_null(pl1.sublist("SubA").getModifier()));
  TEST_EQUALITY(*pl2.sublist("SubA").getModifier(), *modifier1);
  // Now test that `setParametersNotAlreadySet` has the correct behavior
  pl1.sublist("SubA").setModifier(modifier1);
  pl1.sublist("SubB", modifier2).set("C", 3);
  pl2.setParametersNotAlreadySet(pl1);
  TEST_EQUALITY(*pl2.sublist("SubB").getModifier(), *modifier2);
  TEST_EQUALITY(pl1, pl2);
  // Test that sublists with their modifiers and parameters are correctly
  // overwritten in `setParameters`
  pl1 = ParameterList();
  pl1.sublist("SubA", modifier1).set("A", 1);
  pl2 = ParameterList();
  pl2.sublist("SubA", modifier2).set("B", 2);
  pl2.setParameters(pl1);
  // pl2 should look just like pl1 except with the extra "B" parameter
  ParameterList pl_expected = ParameterList();
  pl_expected.sublist("SubA", modifier1).set("B", 2).set("A", 1);
  TEST_EQUALITY(pl_expected, pl2);
}


TEUCHOS_UNIT_TEST( ParameterList, param_isParameter_isSublist_isType )
{
  ParameterList pl;
  ECHO(pl.set("my int", 3));
  ECHO(const int my_int = pl.get<int>("my int"));
  TEST_EQUALITY_CONST(my_int, 3);
  TEST_EQUALITY_CONST(pl.isParameter("my int"), true);
  TEST_EQUALITY_CONST(pl.isParameter("Does not Exist"), false);
  TEST_EQUALITY_CONST(pl.isSublist("my int"), false);
  TEST_EQUALITY_CONST(pl.isSublist("Does not exist"), false);
  TEST_EQUALITY_CONST(pl.isType<int>("my int"), true);
  TEST_EQUALITY_CONST(pl.isType<double>("my int"), false);
  TEST_EQUALITY_CONST(pl.isType("my int", static_cast<int*>(0)), true);
  TEST_EQUALITY_CONST(pl.isType("my int", static_cast<double*>(0)), false);
}


TEUCHOS_UNIT_TEST( ParameterList, sublist_isParameter_isSublist_isType )
{
  ParameterList pl;
  ECHO(pl.sublist("my sublist").set("my int", 3));
  ECHO(const int my_int = getConst(pl).sublist("my sublist").get<int>("my int"));
  TEST_EQUALITY_CONST(my_int, 3);
  TEST_EQUALITY_CONST(pl.isParameter("my sublist"), true); // Should be false, but backward compatiable!
  TEST_EQUALITY_CONST(pl.isParameter("Does not Exist"), false);
  TEST_EQUALITY_CONST(pl.isSublist("my sublist"), true);
  TEST_EQUALITY_CONST(pl.isType<ParameterList>("my sublist"), true);
  TEST_EQUALITY_CONST(pl.isType<double>("my sublist"), false);
  TEST_EQUALITY_CONST(pl.isType("my sublist", static_cast<ParameterList*>(0)), true);
  TEST_EQUALITY_CONST(pl.isType("my sublist", static_cast<double*>(0)), false);
}


TEUCHOS_UNIT_TEST( ParameterList, set_doc )
{
  ParameterList pl;
  ECHO(pl.set("my int", 3, "Some documentation"));
  ECHO(const ParameterEntry& my_int_param = getConst(pl).getEntry("my int"));
  TEST_EQUALITY_CONST(my_int_param.docString(), "Some documentation");
  TEST_ASSERT(is_null(my_int_param.validator()));
}


TEUCHOS_UNIT_TEST( ParameterList, set_doc_validator )
{
  ParameterList pl;
  ECHO(pl.set("my int", 3, "Some documentation", rcp(new DummyValidator)));
  ECHO(const ParameterEntry& my_int_param = getConst(pl).getEntry("my int"));
  TEST_EQUALITY_CONST(my_int_param.docString(), "Some documentation");
  TEST_NOTHROW(rcp_dynamic_cast<const DummyValidator>(my_int_param.validator(), true));
}


TEUCHOS_UNIT_TEST( ParameterList, set_invalid_int_first )
{
  ParameterList pl;
  ECHO(const RCP<ParameterEntryValidator>
    validator(new Teuchos::EnhancedNumberValidator<int>(0, 1)));
  TEST_THROW(pl.set("my int", -1, "", validator),
    Exceptions::InvalidParameterValue);
  TEST_EQUALITY_CONST(pl.numParams(), 0);
}


TEUCHOS_UNIT_TEST( ParameterList, set_invalid_int_second )
{
  ParameterList pl;
  ECHO(const RCP<ParameterEntryValidator>
    validator(new Teuchos::EnhancedNumberValidator<int>(0, 1)));
  TEST_NOTHROW(pl.set("my int", 1, "", validator));
  TEST_EQUALITY_CONST(pl.numParams(), 1);
  TEST_EQUALITY_CONST(pl.get<int>("my int"), 1);
  TEST_THROW(pl.set("my int", -1), Exceptions::InvalidParameterValue);
  TEST_EQUALITY_CONST(pl.get<int>("my int"), 1);
}


TEUCHOS_UNIT_TEST( ParameterList, set_get_entry )
{
  ParameterList pl;
  ECHO(pl.setEntry("my int", ParameterEntry(as<int>(3), true, true, "Some doc", rcp(new DummyValidator))));
  ECHO(const ParameterEntry& my_int_param = getConst(pl).getEntry("my int"));
  TEST_EQUALITY_CONST(my_int_param.docString(), "Some doc");
  ECHO(const int my_int_1 = my_int_param.getValue<int>(0));
  TEST_EQUALITY_CONST(my_int_1, 3);
  TEST_EQUALITY_CONST(my_int_param.isUsed(), true);
  TEST_EQUALITY_CONST(my_int_param.isList(), false); // The isList entry is ignored!
  TEST_INEQUALITY_CONST(rcp_dynamic_cast<const DummyValidator>(my_int_param.validator(), true), null);
}


TEUCHOS_UNIT_TEST( ParameterList, set_int_twice_keep_validator )
{
  ParameterList pl;
  ECHO(pl.setEntry("my int", ParameterEntry(as<int>(3), true, true, "Some doc", rcp(new DummyValidator))));
  {
    ECHO(const ParameterEntry& my_int_param = getConst(pl).getEntry("my int"));
    TEST_INEQUALITY_CONST(rcp_dynamic_cast<const DummyValidator>(my_int_param.validator(), true), null);
  }
  TEST_EQUALITY_CONST(pl.get<int>("my int"), 3);
  ECHO(pl.set("my int", 4));
  TEST_EQUALITY_CONST(pl.get<int>("my int"), 4);
  {
    ECHO(const ParameterEntry& my_int_param = getConst(pl).getEntry("my int"));
    TEST_INEQUALITY_CONST(rcp_dynamic_cast<const DummyValidator>(my_int_param.validator(), true), null);
  }
}


TEUCHOS_UNIT_TEST( ParameterList, set_get_char_str )
{
  ParameterList pl;

  ECHO(char dummy_str_1[] = "dummy str 1");
  ECHO(pl.set("dummy 1", dummy_str_1));
  ECHO(const std::string dummy_1 = pl.get<std::string>("dummy 1"));
  TEST_EQUALITY_CONST(dummy_1, "dummy str 1");

  ECHO(const char dummy_str_const_2[] = "dummy str 2");
  ECHO(pl.set("dummy 2", dummy_str_const_2));
  ECHO(const std::string dummy_2 = pl.get<std::string>("dummy 2"));
  TEST_EQUALITY_CONST(dummy_2, "dummy str 2");

}


TEUCHOS_UNIT_TEST( ParameterList, set_get_string )
{
  ParameterList pl;

  ECHO(const std::string dummy_str = "dummy str");
  ECHO(pl.set("my str", dummy_str));
  ECHO(const std::string my_str = pl.get<std::string>("my str"));
  TEST_EQUALITY_CONST(my_str, "dummy str");

}

TEUCHOS_UNIT_TEST( ParameterList, set_string_move_semantics)
{
  ParameterList pl;

  ECHO(std::string my_str_1{"my text 1"});
  ECHO(pl.set("my string 1", std::move(my_str_1)));

  // Check that the parameter value was moved by checking that my_str_1 is now empty.
  TEST_ASSERT(my_str_1.empty());

  TEST_EQUALITY_CONST(pl.get<std::string>("my string 1"), "my text 1");
}

TEUCHOS_UNIT_TEST( ParameterList, set_string_specified_template_argument)
{
  // Check the templated set method and its overload when the template argument is specified.

   ParameterList pl;

  // The parameter value can be passed by rvalue reference.
  // The main templated set method is called, and the parameter value is moved.
  ECHO(std::string my_str_2{"my text 2"});
  ECHO(pl.set<std::string>("my string 2", std::move(my_str_2)));
  TEST_ASSERT(my_str_2.empty());
  TEST_EQUALITY_CONST(pl.get<std::string>("my string 2"), "my text 2");
  
  // The parameter value cannot be passed by rvalue reference.
  // The overload of the templated set method is called, and the parameter value is not moved.
  ECHO(std::string my_str_3{"my text 3"});
  ECHO(pl.set<std::string>("my string 3", my_str_3));
  TEST_ASSERT( ! my_str_3.empty());
  TEST_EQUALITY_CONST(pl.get<std::string>("my string 3"), "my text 3");

  ECHO(const std::string my_str_4{"my text 4"});
  ECHO(pl.set<std::string>("my string 4", my_str_4));
  TEST_ASSERT( ! my_str_4.empty());
  TEST_EQUALITY_CONST(pl.get<std::string>("my string 4"), "my text 4");

  ECHO(pl.set<std::string>("my string 5", "my text 5"));
  TEST_EQUALITY_CONST(pl.get<std::string>("my string 5"), "my text 5");
}

struct MyWrappedInt
{
    operator const int&() const { return value; }

    int value;
};

TEUCHOS_UNIT_TEST( ParameterList, set_int_user_defined_conversion_function)
{
  // Check the templated set method and its overload when the template argument is specified
  // and the parameter is set via a user defined conversion function. 
  ParameterList pl;

  ECHO(MyWrappedInt my_wrapped_int{42});
  ECHO(pl.set<int>("my int", my_wrapped_int));
  TEST_EQUALITY_CONST(pl.get<int>("my int"), my_wrapped_int.value);
}

TEUCHOS_UNIT_TEST( ParameterList, get_nonexisting_param )
{
  ParameterList pl;
  TEST_THROW(pl.getEntry("Does not exist 1"), Exceptions::InvalidParameterName);
  TEST_THROW(pl.get<int>("Does not exist 2"), Exceptions::InvalidParameterName);
  TEST_THROW(getConst(pl).get<int>("Does not exist 3"), Exceptions::InvalidParameterName);
  TEST_EQUALITY(pl.getPtr<int>("Does not exist 4"), static_cast<int*>(0));
  TEST_EQUALITY(getConst(pl).getPtr<int>("Does not exist 5"), static_cast<const int*>(0));
  ECHO(char raw_str[] = "dummy");
  TEST_EQUALITY_CONST(pl.get("Does not exist 6", raw_str), "dummy");
  ECHO(const char raw_c_str[] = "dummy");
  TEST_EQUALITY_CONST(pl.get("Does not exist 7", raw_c_str), "dummy");
  ECHO(const std::string str = "dummy");
  TEST_EQUALITY_CONST(pl.get("Does not exist 8", str), "dummy");
  TEST_THROW(pl.getEntry("Does not exist 9"), Exceptions::InvalidParameterName);
  TEST_THROW(getConst(pl).getEntry("Does not exist 10"), Exceptions::InvalidParameterName);
  TEST_EQUALITY(pl.getEntryPtr("Does not exist 11"), static_cast<ParameterEntry*>(0));
  TEST_EQUALITY(getConst(pl).getEntryPtr("Does not exist 12"), static_cast<const ParameterEntry*>(0));
  TEST_EQUALITY(pl.getEntryRCP("Does not exist 13"), RCP<ParameterEntry>());
  TEST_EQUALITY(getConst(pl).getEntryRCP("Does not exist 14"), RCP<const ParameterEntry>());
}


TEUCHOS_UNIT_TEST( ParameterList, get_existing_incorrect_type )
{
  ParameterList pl;
  pl.set("my int", 4);
  TEST_THROW(pl.get<double>("my int"), Exceptions::InvalidParameterType);
  // ToDo: Assert the contents of the error message
}


TEUCHOS_UNIT_TEST( ParameterList, getPtr )
{
  ParameterList pl;
  pl.set("my int", 4);
  TEST_EQUALITY_CONST(pl.getPtr<int>("Does not Exist"), static_cast<int*>(0));
  TEST_INEQUALITY_CONST(pl.getPtr<int>("my int"), static_cast<int*>(0));
  TEST_EQUALITY_CONST(*pl.getPtr<int>("my int"), 4);
  TEST_EQUALITY_CONST(pl.getPtr<double>("my int"), static_cast<double*>(0));
  TEST_EQUALITY_CONST(getConst(pl).getPtr<int>("Does not Exist"), static_cast<const int*>(0));
  TEST_INEQUALITY_CONST(getConst(pl).getPtr<int>("my int"), static_cast<int*>(0));
  TEST_EQUALITY_CONST(*getConst(pl).getPtr<int>("my int"), 4);
  TEST_EQUALITY_CONST(getConst(pl).getPtr<double>("my int"), static_cast<const double*>(0));
}


TEUCHOS_UNIT_TEST( ParameterList, getEntryRCP )
{
  ParameterList pl;
  pl.set("my int", 4);
  TEST_EQUALITY_CONST(pl.getEntryRCP("Does not Exist"), null);
  TEST_INEQUALITY_CONST(pl.getEntryRCP("my int"), null);
  TEST_EQUALITY_CONST(pl.getEntryRCP("my int")->getValue<int>(0), 4);
  TEST_EQUALITY_CONST(getConst(pl).getEntryRCP("Does not Exist"), null);
  TEST_INEQUALITY_CONST(getConst(pl).getEntryRCP("my int"), null);
  TEST_EQUALITY_CONST(getConst(pl).getEntryRCP("my int")->getValue<int>(0), 4);
}


// Test nonconstFind()

// Test find()


TEUCHOS_UNIT_TEST( ParameterList, get_default_then_change )
{
  ParameterList pl;
  ECHO(int &my_int = pl.get("my int", 3));
  TEST_EQUALITY_CONST(my_int, 3);
  TEST_EQUALITY_CONST(pl.get<int>("my int"), 3);
  ECHO(my_int = 5);
  TEST_EQUALITY_CONST(pl.get<int>("my int"), 5);
}


TEUCHOS_UNIT_TEST( ParameterList, remove_1 )
{
  ParameterList pl;
  TEST_EQUALITY_CONST(pl.numParams(), 0);
  ECHO(pl.set("my int", 2));
  TEST_EQUALITY_CONST(pl.numParams(), 1);
  TEST_EQUALITY_CONST(pl.get<int>("my int"), 2);
  ECHO(const bool param_was_removed_1 = pl.remove("my int"));
  TEST_EQUALITY_CONST(param_was_removed_1, true);
  TEST_EQUALITY_CONST(pl.numParams(), 0);
  TEST_THROW(pl.get<int>("my int"), Exceptions::InvalidParameterName);
  TEST_THROW(pl.remove("my int"), Exceptions::InvalidParameterName);
  ECHO(const bool param_was_removed_2 = pl.remove("my int", false));
  TEST_EQUALITY_CONST(param_was_removed_2, false);
}


TEUCHOS_UNIT_TEST( ParameterList, get_nonexisting_sublist_default )
{
  ParameterList pl("Base");
  ECHO(pl.sublist("my sublist"));
  ECHO(const ParameterEntry &sublistParam = pl.getEntry("my sublist"));
  TEST_EQUALITY_CONST(sublistParam.isUsed(), false);
  TEST_EQUALITY_CONST(sublistParam.isList(), true);
  TEST_EQUALITY_CONST(sublistParam.isDefault(), false);
  TEST_EQUALITY_CONST(sublistParam.docString(), "");
  TEST_EQUALITY_CONST(sublistParam.getValue<ParameterList>(0).name(), "Base->my sublist");
}


TEUCHOS_UNIT_TEST( ParameterList, get_nonexisting_sublist_docString )
{
  ParameterList pl("Base");
  ECHO(pl.sublist("my sublist", false, "My great sublist"));
  ECHO(const ParameterEntry &sublistParam = pl.getEntry("my sublist"));
  TEST_EQUALITY_CONST(sublistParam.isUsed(), false);
  TEST_EQUALITY_CONST(sublistParam.isList(), true);
  TEST_EQUALITY_CONST(sublistParam.isDefault(), false);
  TEST_EQUALITY_CONST(sublistParam.docString(), "My great sublist");
  TEST_EQUALITY_CONST(sublistParam.getValue<ParameterList>(0).name(), "Base->my sublist");
}


TEUCHOS_UNIT_TEST( ParameterList, get_nonexisting_sublist_mustAlreadyExist )
{
  ParameterList pl("Base");
  TEST_THROW(pl.sublist("my sublist", true), Exceptions::InvalidParameterName);
  // ToDo: Examine the actual structure of the error message
}


TEUCHOS_UNIT_TEST( ParameterList, get_existing_sublist_nonsublist )
{
  ParameterList pl("Base");
  ECHO(pl.set("my sublist", 1)); // Not a sublist!
  TEST_THROW(pl.sublist("my sublist"), Exceptions::InvalidParameterType);
  // ToDo: Examine the actual structure of the error message
}


TEUCHOS_UNIT_TEST( ParameterList, get_existing_sublist_nonconst )
{
  ParameterList pl("Base");
  ECHO(pl.sublist("my sublist").set("my int", 2));
  ECHO(const int my_int = pl.sublist("my sublist").get<int>("my int"));
  TEST_EQUALITY_CONST(my_int, 2);
}


TEUCHOS_UNIT_TEST( ParameterList, get_existing_sublist_const )
{
  ParameterList pl("Base");
  ECHO(pl.sublist("my sublist").set("my int", 2));
  ECHO(const int my_int = getConst(pl).sublist("my sublist").get<int>("my int"));
  TEST_EQUALITY_CONST(my_int, 2);
}


TEUCHOS_UNIT_TEST( ParameterList, get_nonexisting_sublist_const )
{
  ParameterList pl("Base");
  TEST_THROW(getConst(pl).sublist("my sublist"), Exceptions::InvalidParameterName);
  // ToDo: Examine the actual structure of the error message
}


TEUCHOS_UNIT_TEST( ParameterList, get_existing_sublist_const_nonsublist )
{
  ParameterList pl("Base");
  ECHO(pl.set("my sublist", 1)); // Not a sublist!
  TEST_THROW(getConst(pl).sublist("my sublist"), Exceptions::InvalidParameterType);
  // ToDo: Examine the actual structure of the error message
}


TEUCHOS_UNIT_TEST( ParameterList, sublist_add_2 )
{
  ParameterList PL_Main("PL_Main");
  const std::string Direction_Doc = "This sublist controls how direction is computed.";
  ParameterList& PL_Direction = PL_Main.sublist("Direction", false, Direction_Doc);
  ParameterList& PL_LineSearch = PL_Main.sublist("Line Search");
  out << "PL_Main=\n" << PL_Main << "\n";
  TEST_EQUALITY_CONST(PL_Main.name(), "PL_Main");
  TEST_EQUALITY_CONST(PL_Main.isSublist("Direction"), true);
  TEST_EQUALITY_CONST(PL_Main.isSublist("Line Search"), true);
  ECHO(const ParameterList& PL_Direction_2 = getConst(PL_Main).sublist("Direction"));
  TEST_EQUALITY(&PL_Direction, &PL_Direction_2);
  ECHO(const ParameterList& PL_LineSearch_2 = getConst(PL_Main).sublist("Line Search"));
  TEST_EQUALITY(&PL_LineSearch, &PL_LineSearch_2);
  TEST_EQUALITY_CONST(getConst(PL_Main).sublist("Direction").name(), "PL_Main->Direction");
  TEST_EQUALITY_CONST(PL_Direction.name(), "PL_Main->Direction");
  TEST_EQUALITY_CONST(PL_LineSearch.name(), "PL_Main->Line Search");
}


TEUCHOS_UNIT_TEST( ParameterList, sublist_scenario_1 )
{
  // This is the scenario in the orginal testing program
  ParameterList PL_Main("PL_Main");
  const std::string Direction_Doc = "This sublist controls how direction is computed.";
  ParameterList &PL_Direction = PL_Main.sublist("Direction", false, Direction_Doc);
  ParameterList &PL_Newton = PL_Direction.sublist("Newton");
  ParameterList &PL_LinSol = PL_Newton.sublist("Linear Solver");
  ParameterList &PL_LineSearch = PL_Main.sublist("Line Search");
  out << "PL_Main=\n" << PL_Main << "\n";
  TEST_EQUALITY_CONST(PL_Main.name(), "PL_Main");
  TEST_EQUALITY_CONST(PL_Main.isSublist("Direction"), true);
  ECHO(const ParameterList& PL_Direction_2 = getConst(PL_Main).sublist("Direction"));
  TEST_EQUALITY(&PL_Direction, &PL_Direction_2);
  TEST_EQUALITY_CONST(getConst(PL_Main).sublist("Direction").name(), "PL_Main->Direction");
  TEST_EQUALITY_CONST(PL_Direction.name(), "PL_Main->Direction");
  TEST_EQUALITY_CONST(PL_Direction.isSublist("Newton"), true);
  TEST_EQUALITY_CONST(PL_Newton.isSublist("Linear Solver"), true);
  TEST_EQUALITY_CONST(PL_Newton.name(), "PL_Main->Direction->Newton");
  TEST_EQUALITY_CONST(PL_Newton.isSublist("Linear Solver"), true);
  TEST_EQUALITY_CONST(PL_LinSol.name(), "PL_Main->Direction->Newton->Linear Solver");
  TEST_EQUALITY_CONST(PL_Main.isSublist("Line Search"), true);
  TEST_EQUALITY_CONST(PL_LineSearch.name(), "PL_Main->Line Search");
}


TEUCHOS_UNIT_TEST( ParameterList, copy_constructor )
{
  ECHO(ParameterList pl1("A"));
  ECHO(pl1.set("my int", 2));
  ECHO(ParameterList pl2(pl1));
  TEST_EQUALITY_CONST(pl2.name(), "A");
  TEST_EQUALITY_CONST(pl2.get<int>("my int"), 2);
}


TEUCHOS_UNIT_TEST( ParameterList, assignment_operator )
{
  ECHO(ParameterList pl1("A"));
  ECHO(pl1.set("my int", 2));
  ECHO(ParameterList pl2);
  ECHO(const ParameterList &pl2_ref = pl2 = pl1);
  TEST_EQUALITY_CONST(&pl2_ref, &pl2);
  TEST_EQUALITY_CONST(pl2.name(), "A");
  TEST_EQUALITY_CONST(pl2.get<int>("my int"), 2);
}


TEUCHOS_UNIT_TEST( ParameterList, iterator_params )
{
  typedef ParameterList::ConstIterator ConstIter;
  ParameterList pl;
  pl.set("c", 1);
  pl.set("a", 2);
  pl.set("b", 3);
  ConstIter pl_itr = pl.begin();
  TEST_EQUALITY_CONST(pl_itr->first, "c");
  TEST_EQUALITY_CONST(pl_itr->second.getValue<int>(0), 1);
  ECHO(++pl_itr);
  TEST_EQUALITY_CONST(pl_itr->first, "a");
  TEST_EQUALITY_CONST(pl_itr->second.getValue<int>(0), 2);
  ECHO(++pl_itr);
  TEST_EQUALITY_CONST(pl_itr->first, "b");
  TEST_EQUALITY_CONST(pl_itr->second.getValue<int>(0), 3);
  ECHO(++pl_itr);
  TEST_ITER_EQUALITY(pl_itr, pl.end());
}


TEUCHOS_UNIT_TEST( ParameterList, iterator_params_sublists )
{
  typedef ParameterList::ConstIterator ConstIter;
  ParameterList pl("base");
  pl.set("c", 1);
  pl.sublist("a");
  pl.set("b", 3);
  ConstIter pl_itr = pl.begin();
  TEST_EQUALITY_CONST(pl_itr->first, "c");
  TEST_EQUALITY_CONST(pl_itr->second.getValue<int>(0), 1);
  ECHO(++pl_itr);
  TEST_EQUALITY_CONST(pl_itr->first, "a");
  TEST_EQUALITY_CONST(pl_itr->second.getValue<ParameterList>(0).name(), "base->a");
  ECHO(++pl_itr);
  TEST_EQUALITY_CONST(pl_itr->first, "b");
  TEST_EQUALITY_CONST(pl_itr->second.getValue<int>(0), 3);
  ECHO(++pl_itr);
  TEST_ITER_EQUALITY(pl_itr, pl.end());
}

// Test iterator access after removing params

// Test iterator access after removing sublists


TEUCHOS_UNIT_TEST( ParameterList, operatorEqualityWithEmpty )
{
  // An empty list should not be equal to a full list
  ParameterList A;
  ParameterList B;
  TEST_ASSERT( A == B );
  A.set("Hello","World");
  TEST_ASSERT( A != B );
  B.set("Hello","World");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( ParameterList, operatorEqualityDifferentSublistNames )
{
  // Sublists with different names should not be equal
  ParameterList A;
  ParameterList B;
  A.sublist("Bob");
  B.sublist("Tom");
  TEST_ASSERT( A != B );
}


TEUCHOS_UNIT_TEST( ParameterList, operatorEqualityDifferentLengths )
{
  ParameterList A;
  ParameterList B;
  A.set("A","a");
  A.set("B","b");
  A.set("C","c");
  A.print(out);

  B.set("A","a");
  B.set("B","b");
  B.print(out);

  TEST_ASSERT( A != B );

  B.set("C","c");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( ParameterList, haveSameValuesWithEmpty )
{
  ParameterList A;
  ParameterList B;
  TEST_ASSERT( haveSameValues(A,B) );
  A.set("a",1);
  TEST_ASSERT( !haveSameValues(A,B) );
  A.set("b",2);
  B.set("a",1);
  B.set("b",2);
  TEST_ASSERT( haveSameValues(A,B) );
}


TEUCHOS_UNIT_TEST( ParameterList, haveSameValuesDifferentSublistNames )
{
  ParameterList A;
  ParameterList B;
  A.sublist("Smith").set("People",4);
  B.sublist("Jones").set("People",4);
  TEST_ASSERT( !haveSameValues(A,B) ); // sublist names matter
}


TEUCHOS_UNIT_TEST( ParameterList, haveSameValuesSortedReversedOrder )
{
  ParameterList A, B;
  A.set("a",1);
  A.set("b",2);
  // Create second list with the same entries but different order
  B.set("b",2);
  B.set("a",1);
  TEST_ASSERT( haveSameValuesSorted(A,B) );
  B.set("c",3);
  TEST_ASSERT( !haveSameValuesSorted(A,B) ); // check for length
}


TEUCHOS_UNIT_TEST( ParameterList, haveSameValuesSortedNested)
{
  ParameterList A, B;
  ParameterList &asublist = A.sublist("A");
  asublist.set("a",1);
  asublist.set("b",2);
  ParameterList &bsublist = B.sublist("A");
  bsublist.set("a",1);
  bsublist.set("b",2);
  TEST_ASSERT( haveSameValuesSorted(A,B) );
  asublist.set("c",3);
  bsublist.set("c",4);
  TEST_ASSERT( !haveSameValuesSorted(A,B) );
}


TEUCHOS_UNIT_TEST( ParameterList, validateAgainstSelf )
{
  ParameterList PL_Main = createMainPL();
  ParameterList PL_Main_valid = createValidMainPL();
  TEST_NOTHROW(PL_Main.validateParameters(PL_Main_valid));
}


TEUCHOS_UNIT_TEST( ParameterList, validateParametersAndSetDefaults_default )
{
  // Test for proper behavior when the user doesn't set `Nonlinear Solver`
  ParameterList PL_Main = createMainPL();
  ParameterList PL_Main_valid = createValidMainPL();
  ECHO(PL_Main.validateParametersAndSetDefaults(PL_Main_valid));
  TEST_NOTHROW(
    rcp_dynamic_cast<const StringToIntegralParameterEntryValidator<int> >(
      PL_Main.getEntry("Nonlinear Solver").validator(), true ) );
  // Make sure the parameter entry is set to default and unused after validation
  const ParameterEntry &default_entry = PL_Main.getEntry("Nonlinear Solver");
  TEST_EQUALITY(default_entry.isDefault(), true);
  TEST_EQUALITY(default_entry.isUsed(), false);
  // Make sure the value is stored as an integer after validation
#if defined(HAVE_TEUCHOS_MODIFY_DEFAULTS_DURING_VALIDATION)
  TEST_NOTHROW(Teuchos::any_cast<int>(default_entry.getAny()));
#endif
}


TEUCHOS_UNIT_TEST( ParameterList, validateParametersAndSetDefaults_noDefault )
{
  // Now make sure we have the correct behavior when not using a default value
  ParameterList PL_Main = createMainPL();
  PL_Main.set("Nonlinear Solver", "Trust Region Based");
  ParameterList PL_Main_valid = createValidMainPL();
  PL_Main.validateParametersAndSetDefaults(PL_Main_valid);
  const ParameterEntry &entry = PL_Main.getEntry("Nonlinear Solver");
  TEST_EQUALITY(entry.isDefault(), false);
  TEST_EQUALITY(entry.isUsed(), false);
#if defined(HAVE_TEUCHOS_MODIFY_DEFAULTS_DURING_VALIDATION)
  TEST_NOTHROW(Teuchos::any_cast<int>(entry.getAny()));
#endif
}


TEUCHOS_UNIT_TEST( ParameterList, getIntegralValue_int )
{
  ParameterList PL_Main = createMainPL();
  ParameterList PL_Main_valid = createValidMainPL();
  ECHO(PL_Main.set("Nonlinear Solver", "Line Search Based"));
  ECHO(PL_Main.validateParametersAndSetDefaults(PL_Main_valid));
  ECHO(const int lineSearchValue = getIntegralValue<int>(PL_Main, "Nonlinear Solver"));
  TEST_EQUALITY_CONST(lineSearchValue, 0);
  ECHO(PL_Main.set("Nonlinear Solver", "Trust Region Based"));
  ECHO(const int trustRegionValue  = getIntegralValue<int>(PL_Main, "Nonlinear Solver"));
  TEST_EQUALITY_CONST(trustRegionValue, 1);
}


TEUCHOS_UNIT_TEST( ParameterList, replaceScalarParameterWithArray ) {
  ParameterList pl = ParameterList("Parameter List with Scalar Parameter");
  const int a_val = 2, b_val = 3;
  pl.set("A", a_val);
  pl.set("B", b_val);
  replaceParameterWithArray<int>("A", "A array", pl);
  replaceParameterWithArray<int>("B", "B", pl);
  ParameterList expected_pl = ParameterList("Parameter List with Array Parameter");
  Array<int> a_array = tuple<int>(a_val), b_array = tuple<int>(b_val);
  expected_pl.set("A array", a_array);
  expected_pl.set("B", b_array);
  TEST_ASSERT(haveSameValuesSorted(expected_pl, pl, true));
  // Throw an error when trying to overwrite a parameter that already exists but
  // doesn't have the same name.
  pl.set("C", 1);
  TEST_THROW(replaceParameterWithArray<int>("C", "B", pl), std::logic_error);
  pl.print();
}


TEUCHOS_UNIT_TEST( ParameterList, simpleModifierModifyReconcile )
{
  RCP<SimpleModifier> modifier = rcp(new SimpleModifier());
  ParameterList valid_pl("My Valid Parameter List with a Modifier", modifier);
  //valid_pl before modification
  //  A: 1.0
  //  B: 0.1
  //  SubA:
  //    C: 1
  valid_pl.set("A", 1.0);
  valid_pl.set("B", 0.1);
  valid_pl.sublist("SubA").set("C", 1);
  ParameterList pl("My Parameter List");
  pl.set("A", 5.0);
  pl.set("B", -0.1);
  pl.sublist("SubA 1").set("C", 3);
  pl.sublist("SubA 2").set("C", 4);
  ParameterList expected_valid_pl(valid_pl);
  expected_valid_pl.remove("SubA");
  expected_valid_pl.sublist("SubA 1").set("C", 1);
  expected_valid_pl.sublist("SubA 2").set("C", 1);
  pl.modifyParameterList(valid_pl);
  //valid_pl after modification
  //  A: 1.0
  //  B: 0.1
  //  SubA 1:
  //    C: 1
  //  SubA 2:
  //    C: 1
  TEST_EQUALITY(valid_pl, expected_valid_pl);
//  std::cout << haveSameValuesSorted(expected_valid_pl, valid_pl, true) << std::endl;
  pl.validateParametersAndSetDefaults(valid_pl);
  TEST_NOTHROW(pl.reconcileParameterList(valid_pl));
  pl.set("A", -1.0);
  TEST_THROW(pl.reconcileParameterList(valid_pl), std::logic_error);
  // Test the copy constructor
  ParameterList copy_valid_pl(valid_pl);
  TEST_EQUALITY(valid_pl, copy_valid_pl);
}


TEUCHOS_UNIT_TEST( ParameterList, modify_CopiesModifiers ) {
  RCP<ParameterListModifier> mod_top = rcp<ParameterListModifier>(new ParameterListModifier("Modifier Top"));
  RCP<ParameterListModifier> mod_A = rcp<ParameterListModifier>(new ParameterListModifier("Modifier A"));
  RCP<ParameterListModifier> mod_B = rcp<ParameterListModifier>(new ParameterListModifier("Modifier B"));
  ParameterList input{"Plist"}, valid{"Plist", mod_top};
  input.sublist("A").sublist("B");
  valid.sublist("A", mod_A);
  valid.sublist("A").sublist("B", mod_B);
  input.modifyParameterList(valid);
  // This calls `haveSameModifiers` to make sure they are all copied to `input`
  TEST_EQUALITY(input, valid);
}


TEUCHOS_UNIT_TEST( ParameterList, nestedSublistExpansion ) {
  Teuchos::RCP<SimpleModifier> modifier = Teuchos::rcp(new SimpleModifier());
  Teuchos::RCP<SimpleSubModifier> sub_modifier = Teuchos::rcp(new SimpleSubModifier());
  // The unmodified (template-like) validation parameter list
  ParameterList valid_pl("valid_pl", modifier);
  valid_pl.set("A", 1.0);
  valid_pl.set("B", 1.0);
  valid_pl.sublist("SubA").setModifier(sub_modifier);
  valid_pl.sublist("SubA").set("C", 3);
  valid_pl.sublist("SubA").set("D", 4);
  valid_pl.sublist("SubA").sublist("SubB").set("E", 10);
  valid_pl.sublist("SubA").sublist("SubB").set("F", 11);
  // The user's input parameter list
  ParameterList pl("pl");
  pl.set("A", 1.0);
  pl.set("B", 2.0);
  pl.sublist("SubA 1").set("C", 3);
  pl.sublist("SubA 1").set("D", 4);
  pl.sublist("SubA 1").sublist("SubB 1").set("E", 51);
  pl.sublist("SubA 1").sublist("SubB 1").set("F", 61);
  pl.sublist("SubA 1").sublist("SubB 2").set("E", 52);
  pl.sublist("SubA 1").sublist("SubB 2").set("F", 62);
  pl.sublist("SubA 2").set("C", 3);
  pl.sublist("SubA 2").set("D", 4);
  pl.sublist("SubA 2").sublist("SubB 3").set("E", 53);
  pl.sublist("SubA 2").sublist("SubB 3").set("F", 63);
  // The expanded valid parameter list after modification
  ParameterList expected_valid_pl("valid_pl_expanded");
  expected_valid_pl.set("A", 1.0);
  expected_valid_pl.set("B", 1.0);
  expected_valid_pl.sublist("SubA 1").set("C", 3);
  expected_valid_pl.sublist("SubA 1").set("D", 4);
  expected_valid_pl.sublist("SubA 1").sublist("SubB 1").set("E", 10);
  expected_valid_pl.sublist("SubA 1").sublist("SubB 1").set("F", 11);
  expected_valid_pl.sublist("SubA 1").sublist("SubB 2").set("E", 10);
  expected_valid_pl.sublist("SubA 1").sublist("SubB 2").set("F", 11);
  expected_valid_pl.sublist("SubA 2").set("C", 3);
  expected_valid_pl.sublist("SubA 2").set("D", 4);
  expected_valid_pl.sublist("SubA 2").sublist("SubB 3").set("E", 10);
  expected_valid_pl.sublist("SubA 2").sublist("SubB 3").set("F", 11);
  // Expand the validation parameter list based on the user's input parameter list
  pl.modifyParameterList(valid_pl);
  // Modified parameter lists aren't equal because they don't have the same modifiers
  TEST_ASSERT(valid_pl != expected_valid_pl);
  // Test that they are the same except for the modifiers
  TEST_ASSERT(haveSameValuesSorted(expected_valid_pl, valid_pl, true));
  // Check the equality of the modifiers
  expected_valid_pl.setModifier(modifier);
  expected_valid_pl.sublist("SubA 1", true).setModifier(sub_modifier);
  expected_valid_pl.sublist("SubA 2", true).setModifier(sub_modifier);
  TEST_ASSERT(haveSameModifiers(valid_pl, expected_valid_pl));
  // Now test the recursive reconciliation
  TEST_NOTHROW(pl.reconcileParameterList(valid_pl));
  pl.sublist("SubA 1").set("C", 11);
  pl.sublist("SubA 1").set("D", 11);
  TEST_THROW(pl.reconcileParameterList(valid_pl), std::logic_error);
}


TEUCHOS_UNIT_TEST( ParameterList, disableRecursion ) {
  Teuchos::RCP<SimpleModifier> modifier = Teuchos::rcp(new SimpleModifier());
  Teuchos::RCP<SimpleSubModifier> sub_modifier = Teuchos::rcp(new SimpleSubModifier());
  // The unmodified (template-like) validation parameter list
  ParameterList valid_pl("valid_pl", modifier);
  valid_pl.set("A", 1.0);
  valid_pl.set("B", 1.0);
  valid_pl.sublist("SubA").setModifier(sub_modifier);
  valid_pl.sublist("SubA").set("C", 3.0);
  valid_pl.sublist("SubA").set("D", 4);
  valid_pl.sublist("SubA").sublist("SubB").set("E", 10);
  valid_pl.sublist("SubA").sublist("SubB").set("F", 11);
  // The user's input parameter list
  ParameterList pl("pl");
  pl.set("A", 1.0);
  pl.set("B", 2.0);
  pl.sublist("SubA 1").set("C", 3);
  pl.sublist("SubA 1").set("D", 4);
  pl.sublist("SubA 1").sublist("SubB").set("E", 53);
  pl.sublist("SubA 1").sublist("SubB").set("E", 63);
  // The expanded valid parameter list after modification
  ParameterList expected_valid_pl("valid_pl");
  expected_valid_pl.set("A", 1.0);
  expected_valid_pl.set("B", 1.0);
  expected_valid_pl.sublist("SubA 1").set("C", 3.0);
  expected_valid_pl.sublist("SubA 1").set("D", 4);
  expected_valid_pl.sublist("SubA 1").sublist("SubB").set("E", 10);
  expected_valid_pl.sublist("SubA 1").sublist("SubB").set("F", 11);
  // Make a copy of the user's input parameter list before it is validated
  ParameterList copy_pl(pl);
  // The float validator will cast integers in `pl` to floats
  RCP<AnyNumberParameterEntryValidator> float_validator = rcp(
    new AnyNumberParameterEntryValidator(AnyNumberParameterEntryValidator::PREFER_DOUBLE,
        AnyNumberParameterEntryValidator::AcceptedTypes(false).allowInt(true).allowDouble(true)));
  valid_pl.sublist("SubA").getEntry("C").setValidator(float_validator);
  // Don't modify `SubA`
  valid_pl.sublist("SubA").disableRecursiveModification();
  pl.modifyParameterList(valid_pl);
  TEST_ASSERT(haveSameValuesSorted(expected_valid_pl, valid_pl, true));
  // Don't validate `SubA 1`
  valid_pl.sublist("SubA 1").disableRecursiveValidation();
  pl.validateParametersAndSetDefaults(valid_pl);
  // If we were to validate `SubA 1` then parameter C would turn into a float and the following test would fail
  TEST_ASSERT(haveSameValuesSorted(pl, copy_pl, true));
}


TEUCHOS_UNIT_TEST( ParameterList, recursiveValidation ) {
  ParameterList valid_pl("valid_pl");
  valid_pl.set("A", 1);
  valid_pl.sublist("SubA").set("B", 1);
  ParameterList pl("pl");
  pl.set("A", 1.0);
  pl.sublist("SubA").set("B", 2);
  ParameterList validated_pl("valid_pl");
  validated_pl.set("A", 1.0);
  validated_pl.sublist("SubA").set("B", 2.0);
  // The float validator will cast integers in `pl` to floats
  RCP<AnyNumberParameterEntryValidator> float_validator = rcp(
    new AnyNumberParameterEntryValidator(AnyNumberParameterEntryValidator::PREFER_DOUBLE,
        AnyNumberParameterEntryValidator::AcceptedTypes(false).allowInt(true).allowDouble(true)));
  valid_pl.getEntry("A").setValidator(float_validator);
  valid_pl.sublist("SubA").getEntry("B").setValidator(float_validator);
  pl.validateParametersAndSetDefaults(valid_pl);
  // All of the integers in `pl` should be casted to floats in `validated_pl`
  TEST_ASSERT(haveSameValuesSorted(validated_pl, pl, true));
}


TEUCHOS_UNIT_TEST( ParameterList, recursiveReconciliation ) {
  Teuchos::RCP<ReconciliationModifier1> modifier1 = Teuchos::rcp(new ReconciliationModifier1());
  Teuchos::RCP<ReconciliationModifier2> modifier2 = Teuchos::rcp(new ReconciliationModifier2());
  ParameterList valid_pl("valid_pl");
  valid_pl.set("a", 1);
  valid_pl.setModifier(modifier1);
  valid_pl.sublist("A").setModifier(modifier2);
  valid_pl.sublist("A").set("c", 1);
  valid_pl.sublist("A").set("d", 1);
  ParameterList pl("pl");
  pl.set("a", 1);
  pl.sublist("A").set("c", 2);
  pl.sublist("A").set("d", 3);
  ParameterList reconciled_pl("reconciled_pl");
  reconciled_pl.set("a", 1);
  reconciled_pl.set("b", 5);
  reconciled_pl.sublist("A").set("c", 2);
  reconciled_pl.sublist("A").set("d", 3);
  reconciled_pl.sublist("A").set("e", 5);
  pl.reconcileParameterList(valid_pl);
  TEST_ASSERT(haveSameValuesSorted(reconciled_pl, pl, true));
}


TEUCHOS_UNIT_TEST( ParameterList, attachValidatorRecursively ) {
  ParameterList valid_pl("valid_pl");
  valid_pl.set("a", 0.);
  valid_pl.sublist("A").set("b", 0.);
  valid_pl.sublist("A").set("c", 0.);
  valid_pl.sublist("A").sublist("AA").set("d", 0.);
  ParameterList pl("pl");
  pl.set("a", 1);
  pl.sublist("A").set("b", 2);
  pl.sublist("A").set("c", 3);
  pl.sublist("A").sublist("AA").set("d", 4);
  ParameterList validated_pl("validated_pl");
  validated_pl.set("a", 1.);
  validated_pl.sublist("A").set("b", 2.);
  validated_pl.sublist("A").set("c", 3.);
  validated_pl.sublist("A").sublist("AA").set("d", 4.);
  // The float validator will cast integers in `pl` to floats
  RCP<AnyNumberParameterEntryValidator> float_validator = rcp(
    new AnyNumberParameterEntryValidator(AnyNumberParameterEntryValidator::PREFER_DOUBLE,
        AnyNumberParameterEntryValidator::AcceptedTypes(false).allowInt(true).allowDouble(true)));
  valid_pl.recursivelySetValidator<double>(float_validator, 1);
  // This should fail since we only set the float validator on the top level of `valid_pl`
  TEST_THROW(pl.validateParametersAndSetDefaults(valid_pl), std::logic_error);
  // Now attach the validator to every double
  valid_pl.recursivelySetValidator<double>(float_validator);
  pl.validateParametersAndSetDefaults(valid_pl);
  TEST_ASSERT(haveSameValuesSorted(validated_pl, pl, true));
}

// TODO: test printing of modifiers
// valid_pl.print(std::cout, ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true));

TEUCHOS_UNIT_TEST( ParameterList, print ) {
  ParameterList valid_pl("valid_pl");
  valid_pl.set("a", 0.);

  {
    ParameterList pl("pl");
    pl.set("a", 1.);
    pl.validateParametersAndSetDefaults(valid_pl);
    ParameterList::PrintOptions printOptions;
    printOptions.showDefault(false);
    std::ostringstream ss;
    pl.print(ss, printOptions);
    std::cout << ss.str();
    TEST_ASSERT(ss.str().size() > 0);
  }

  {
    ParameterList pl("pl");
    pl.validateParametersAndSetDefaults(valid_pl);
    ParameterList::PrintOptions printOptions;
    std::ostringstream ss;
    pl.print(ss, printOptions);
    std::cout << ss.str();
    TEST_ASSERT(ss.str().size() > 0);

    ss.str("");
    printOptions.showDefault(false);
    pl.print(ss, printOptions);
    std::cout << ss.str();
    TEST_ASSERT(ss.str().size() == 0);
  }
}

// define enum class Shape in anonymous namespace outside of unittest
// "NonPrintableParameterEntries" to avoid polution of class type in string comparison
namespace {
    enum class Shape : int { CIRCLE, SQUARE, TRIANGLE };
}

TEUCHOS_UNIT_TEST( ParameterList, NonPrintableParameterEntries){
  // test printing std::vector<int> from a parameter list
  {
    std::vector<int> testVec = {1};
    ParameterList paramList = ParameterList("std::vector test");
    paramList.set("My std::vector<int>", testVec);

    try {
      paramList.print();  // Should throw!
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "If you get here then the test failed!");
    }
    catch (const NonprintableTypeException &except) {
      std::string actualMessage = except.what();
      std::string expectedMessage = "Trying to print type std::vector<int, std::allocator<int> > "
                                    "which is not printable (i.e. does not have operator<<() defined)!";
      TEST_ASSERT(actualMessage.find(expectedMessage) != std::string::npos);
    }
  }

  // test printing enum class from a parameter list
  {
    ParameterList paramList = ParameterList("enum class test");
    paramList.set("My enum class", Shape::SQUARE);

    try {
      paramList.print();  // Should throw!
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "If you get here then the test failed!" );
    }
    catch (const NonprintableTypeException &except) {
      std::string actualMessage = except.what();
      std::string expectedMessage =
              "Trying to print type Teuchos::(anonymous namespace)::Shape which is not printable "
              "(i.e. does not have operator<<() defined)!";
      TEST_ASSERT(actualMessage.find(expectedMessage) != std::string::npos);
    }
  }
}




} // namespace Teuchos



