// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TwoDArray.hpp"

namespace Teuchos {

//! @test Check that the move constructor of @ref Teuchos::ParameterEntry works as expected.
TEUCHOS_UNIT_TEST(Teuchos_ParameterEntry, move_constructor)
{
  ParameterEntry source(std::string("nice-entry"));
  TEST_EQUALITY(Teuchos::any_cast<std::string>(source.getAny()),"nice-entry");

  //! After the move-construct, @c source has no value (its holding a @c nullptr when looking at how @c Teuchos::any holds values).
  ParameterEntry move_constructed(std::move(source));
  TEST_ASSERT  (source.getAny().access_content() == nullptr);
  TEST_EQUALITY(source.getAny().has_value(),false);
  TEST_EQUALITY(Teuchos::any_cast<std::string>(move_constructed.getAny()),"nice-entry");

  //! Same behavior is expected for the move-assign.
  ParameterEntry move_assigned;
  move_assigned = std::move(move_constructed);
  TEST_ASSERT  (move_constructed.getAny().access_content() == nullptr);
  TEST_EQUALITY(move_constructed.getAny().has_value(),false);
  TEST_EQUALITY(Teuchos::any_cast<std::string>(move_assigned.getAny()),"nice-entry");
}

//! @test This test checks that it is possible to move-extract the value held by a parameter entry.
TEUCHOS_UNIT_TEST(Teuchos_ParameterList, move_extract_value_from_any)
{
  ParameterEntry source(std::string("nice-entry"));

  Teuchos::any::holder<std::string>* dyn_cast_content = dynamic_cast<Teuchos::any::holder<std::string>*>(
    source.getAny().access_content());

  TEST_EQUALITY(dyn_cast_content->held, std::string("nice-entry"));

  //! Simple copy assignment from a reference to the value is still possible and leaves the parameter entry untouched.
  std::string copy_extracted = Teuchos::any_cast<std::string>(source.getAny());
  TEST_EQUALITY(dyn_cast_content->held, std::string("nice-entry"));
  TEST_EQUALITY(copy_extracted        , std::string("nice-entry"));

  //! Move-extract the value will make the string held by the parameter entry empty (default for @c std::string).
  std::string move_extracted = Teuchos::any_cast<std::string>(std::move(source.getAny()));
  TEST_EQUALITY(dyn_cast_content->held, std::string());
  TEST_EQUALITY(move_extracted        , std::string("nice-entry"));
}

TEUCHOS_UNIT_TEST( Teuchos_ParameterEntry, testTypeFunctions )
{
  ParameterEntry intEntry;
  intEntry.setValue(1);
  TEST_ASSERT(intEntry.isType<int>());
  TEST_ASSERT(!intEntry.isArray());
  TEST_ASSERT(!intEntry.isTwoDArray());

  ParameterEntry doubleEntry;
  doubleEntry.setValue(1.0);
  TEST_ASSERT(doubleEntry.isType<double>());
  TEST_ASSERT(!doubleEntry.isArray());
  TEST_ASSERT(!doubleEntry.isTwoDArray());

  Array<int> intArray = tuple<int>(1,2,3);
  ParameterEntry arrayEntry;
  arrayEntry.setValue(intArray);
  TEST_ASSERT(arrayEntry.isType<Array<int> >());
  TEST_ASSERT(arrayEntry.isArray());
  TEST_ASSERT(!arrayEntry.isTwoDArray());

  TwoDArray<double> twoDArray(3,3,3.0);
  ParameterEntry twoDEntry(twoDArray);
  TEST_ASSERT(twoDEntry.isType<TwoDArray<double> >());
  TEST_ASSERT(twoDEntry.isTwoDArray());
  TEST_ASSERT(!twoDEntry.isArray());
}


} // namespace Teuchos



