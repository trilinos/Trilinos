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



