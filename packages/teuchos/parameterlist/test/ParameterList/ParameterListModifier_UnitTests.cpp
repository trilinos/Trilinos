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
#include "Teuchos_ParameterListAcceptorHelpers.hpp"
#include "Teuchos_ParameterListModifier.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_UnitTestHarness.hpp"

//
// Utilities
//


namespace {


class EmptyModifier : public Teuchos::ParameterListModifier {

public:

  EmptyModifier() : Teuchos::ParameterListModifier("Empty Modifier"){}

};


} // namespace


namespace Teuchos {


TEUCHOS_UNIT_TEST( ParameterListModifier, findMatchingBaseNames ){
  ParameterList pl = ParameterList("Parameter List with Matching Base Names");
  const std::string base_name = "foo", non_base_name = "bar";
  const std::string param1 = base_name, param2 = base_name + " 1";
  const std::string sub1 = base_name + "_sublist";
  pl.set(non_base_name, 1);
  pl.set(non_base_name + base_name, 2);
  pl.set(param1, 1);
  pl.set(param2, 2);
  pl.sublist(sub1).set("A", 1);
  RCP<EmptyModifier> empty_modifier = rcp(new EmptyModifier());
  Array<std::string> matches = empty_modifier->findMatchingBaseNames(pl, base_name);
  Array<std::string> expected = tuple(param1, param2, sub1);
  TEST_EQUALITY(matches, expected);
  matches = empty_modifier->findMatchingBaseNames(pl, base_name, true, false);
  expected = Teuchos::tuple(param1, param2);
  TEST_EQUALITY(matches, expected);
  matches = empty_modifier->findMatchingBaseNames(pl, base_name, false, true);
  expected = Teuchos::tuple(sub1);
  TEST_EQUALITY(matches, expected);
}


TEUCHOS_UNIT_TEST( ParameterListModifier, expandParameters ) {
  RCP<EmptyModifier> empty_modifier = rcp(new EmptyModifier());
  ParameterList pl = ParameterList("Parameter List with Expanded Parameters");
  const std::string param_template_name = "Template Parameter";
  ParameterList valid_pl = ParameterList("Parameter List with Template for Parameter Expansion");
  valid_pl.set(param_template_name, 1);
  auto valid_pl_copy = valid_pl;
  pl.set("A", 2);
  pl.set("B", 3);
  empty_modifier->expandParameters(param_template_name, pl, valid_pl);
  ParameterList expected_valid_pl = ParameterList("Parameter List with Template for Parameter Expansion");
  expected_valid_pl.set("A", 1);
  expected_valid_pl.set("B", 1);
  TEST_ASSERT(haveSameValuesSorted(expected_valid_pl, valid_pl, true));
  // Now test for excluding a parameter
  valid_pl = valid_pl_copy;
  pl.set("C", 4);
  pl.set("D", 5);
  empty_modifier->expandParameters(param_template_name, pl, valid_pl, tuple<std::string>("C", "D"));
  TEST_ASSERT(haveSameValuesSorted(expected_valid_pl, valid_pl, true));
}


TEUCHOS_UNIT_TEST( ParameterListModifier, expandSublists ) {
  RCP<EmptyModifier> empty_modifier = rcp(new EmptyModifier());
  ParameterList pl = ParameterList("Parameter List with Expanded Sublists");
  const std::string sublist_template_name = "Template Sublist";
  ParameterList valid_pl = ParameterList("Parameter List with Template for Sublist Expansion");
  // Make sure `expandSublists` ignores normal parameters
  valid_pl.set("var", 1);
  valid_pl.sublist(sublist_template_name);
  auto valid_pl_copy = valid_pl;
  pl.sublist("A");
  pl.sublist("B");
  pl.set("var", 1);
  empty_modifier->expandSublists(sublist_template_name, pl, valid_pl);
  ParameterList expected_valid_pl = ParameterList("Parameter List with Template for Parameter Expansion");
  expected_valid_pl.sublist("A");
  expected_valid_pl.sublist("B");
  expected_valid_pl.set("var", 1);
  TEST_ASSERT(haveSameValuesSorted(expected_valid_pl, valid_pl, true));
  // Now test for excluding a parameter
  valid_pl = valid_pl_copy;
  pl.sublist("C");
  pl.sublist("D");
  empty_modifier->expandSublists(sublist_template_name, pl, valid_pl, tuple<std::string>("C", "D"));
  TEST_ASSERT(haveSameValuesSorted(expected_valid_pl, valid_pl, true));
}


TEUCHOS_UNIT_TEST( ParameterListModifier, expandSublistsUsingBaseName ){
  RCP<EmptyModifier> empty_modifier = rcp(new EmptyModifier());
  ParameterList pl = ParameterList("Parameter List with Expanded Sublists");
  ParameterList valid_pl = ParameterList("Parameter List with Template for Expanding Sublists");
  const std::string base_name = "A";
  auto &vsub = valid_pl.sublist(base_name, empty_modifier, "Sublist A2");
  vsub.set("Val1", 1);
  ParameterList copy_valid_pl(valid_pl);
  pl.sublist("A1");
  pl.sublist("A2");
  empty_modifier->expandSublistsUsingBaseName(base_name, pl, valid_pl);
  ParameterList expected_valid_pl = ParameterList("Parameter List with Template for Expanding Sublists");
  expected_valid_pl.sublist("A1", empty_modifier, "Sublist A2").set("Val1", 1);
  expected_valid_pl.sublist("A2", empty_modifier, "Sublist A2").set("Val1", 1);
  TEST_ASSERT(haveSameValuesSorted(expected_valid_pl, valid_pl, true));
  // Test for enabling the `allow_base_name` option
  pl.sublist(base_name);
  valid_pl = copy_valid_pl;
  expected_valid_pl.sublist(base_name, empty_modifier, "Sublist A2").set("Val1", 1);
  empty_modifier->expandSublistsUsingBaseName(base_name, pl, valid_pl, true);
  TEST_ASSERT(haveSameValuesSorted(expected_valid_pl, valid_pl, true));
  // Test for disabling the `allow_base_name` option
  valid_pl = copy_valid_pl;
  TEST_THROW(empty_modifier->expandSublistsUsingBaseName(base_name, pl, valid_pl, false), std::logic_error);
}


TEUCHOS_UNIT_TEST( ParameterListModifier, setDefaultsInSublists ) {
  RCP<EmptyModifier> empty_modifier = rcp(new EmptyModifier());
  ParameterList pl = ParameterList("Parameter List with Default Parameter");
  const int a_val = 2;
  pl.sublist("AA");
  pl.sublist("AB").set("A", 3);
  ParameterList expected_pl(pl);
  pl.set("A", a_val);
  expected_pl.sublist("AA").set("A", a_val);
  empty_modifier->setDefaultsInSublists("A", pl, tuple<std::string>("AA", "AB"));
  // The AA sublist should change and the "A" parameter should be deleted
  // but the AB sublist should remain the same
  TEST_ASSERT(haveSameValuesSorted(expected_pl, pl, true));
}


} // namespace Teuchos



