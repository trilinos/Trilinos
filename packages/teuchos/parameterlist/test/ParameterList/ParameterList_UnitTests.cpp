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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
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
    linesearchMaxItersValiator = rcp(
      new Teuchos::AnyNumberParameterEntryValidator(
        Teuchos::AnyNumberParameterEntryValidator::PREFER_INT, // Not used here!
        AcceptedTypes(false).allowInt(true).allowDouble(true).allowString(true)
        )
      );
  PL_Main_valid.sublist("Line Search").sublist("Polynomial").set(
    "Max Iters",3
    ,"The maximum number of inner linear search iterations allowed."
    ,linesearchMaxItersValiator
    );

  // Create a validator for the parameter "Direction"->"Newton"->"Linear Solver"->"Tol"
  // that accepts a 'double' or a 'std::string' value!
  typedef Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes AcceptedTypes;
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
  A.set("Hello","World");
  TEST_ASSERT( !haveSameValues(A,B) );
  B.set("Hello","World");
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


TEUCHOS_UNIT_TEST( ParameterList, validateAgainstSelf )
{
  ParameterList PL_Main = createMainPL();
  ParameterList PL_Main_valid = createValidMainPL();
  TEST_NOTHROW(PL_Main.validateParameters(PL_Main_valid));
}


TEUCHOS_UNIT_TEST( ParameterList, validateParametersAndSetDefaults )
{
  ParameterList PL_Main = createMainPL();
  ParameterList PL_Main_valid = createValidMainPL();
  ECHO(PL_Main.validateParametersAndSetDefaults(PL_Main_valid));
  TEST_NOTHROW(
    rcp_dynamic_cast<const StringToIntegralParameterEntryValidator<int> >(
      PL_Main.getEntry("Nonlinear Solver").validator(), true ) );
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


} // namespace Teuchos



