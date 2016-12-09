/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/
#include <iostream>
#include <string>
#include <cstdlib>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Details_Environment.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_CommHelpers.hpp>


namespace {

  using Tpetra::Details::Environment;
  using std::endl;

  // Constants
  std::string TEST_VAR_SET = "TPETRA_ENVIRONMENT_TEST_HELPER_VAR_SET";
  std::string TEST_VAR_NOT_SET = "TPETRA_ENVIRONMENT_TEST_HELPER_VAR_NOT_SET";

  TEUCHOS_STATIC_SETUP()
  {
    setenv(TEST_VAR_SET.c_str(), TEST_VAR_SET.c_str(), 1);
  }

  //
  // UNIT TESTS
  //
  TEUCHOS_UNIT_TEST(Environment, Check_variableExists_1) {
    out << "Check existence of variable that *should* exist" << endl;
    bool env_exists = Environment::getInstance().variableExists(TEST_VAR_SET);
    TEUCHOS_TEST_ASSERT(env_exists, out, success);
  }

  TEUCHOS_UNIT_TEST(Environment, Check_variableExists_2) {
    out << "Check existence of variable that *should not* exist" << endl;
    bool env_does_not_exist =
      ! Environment::getInstance().variableExists(TEST_VAR_NOT_SET);
    TEUCHOS_TEST_ASSERT(env_does_not_exist, out, success);
  }

  TEUCHOS_UNIT_TEST(Environment, Check_getValue_1) {
    out << "Check variable value that *should* exist" << endl;
    std::string envar = Environment::getInstance().getValue(TEST_VAR_SET);
    bool envar_is_correct = envar == TEST_VAR_SET;
    TEUCHOS_TEST_ASSERT(envar_is_correct, out, success);
  }

  TEUCHOS_UNIT_TEST(Environment, Check_getValue_withCacheCheck) {
    out << "Check that nonexistent variable returns an empty string" << endl;
    std::string name = "__AN_OBVIOSULY_FAKE_NAME__";
    std::string envar = Environment::getInstance().getValue(name);
    bool envar_is_empty = envar.empty();
    TEUCHOS_TEST_ASSERT(envar_is_empty, out, success);
    out << "Check that variable is cached (even though it does not exist)";
    bool envar_is_cached = Environment::getInstance().variableIsCached(name);
    TEUCHOS_TEST_ASSERT(envar_is_cached, out, success);
  }

  TEUCHOS_UNIT_TEST(Environment, Check_getValue_WithDefault) {
    out << "Check that nonexistent variable returns requested default" << endl;
    std::string envar = Environment::getInstance().getValue(TEST_VAR_NOT_SET,
                                                            TEST_VAR_NOT_SET);
    bool env_used_requested_default = envar == TEST_VAR_NOT_SET;
    TEUCHOS_TEST_ASSERT(env_used_requested_default, out, success);
  }

  TEUCHOS_UNIT_TEST(Environment, Check_variableIsCached) {
    out << "Check that environment variable cached";
    std::string name = "TPETRA_DEBUG";
    bool envar_is_cached = Environment::getInstance().variableIsCached(name);
    TEUCHOS_TEST_ASSERT(envar_is_cached, out, success);
    name = "TPETRA_USE_BLAS";
    envar_is_cached = Environment::getInstance().variableIsCached(name);
    TEUCHOS_TEST_ASSERT(envar_is_cached, out, success);
  }

  TEUCHOS_UNIT_TEST(Environment, Check_getBooleanValue_1) {
    out << "Check existing variable returns true" << endl;
    bool envar_is_true = Environment::getInstance().getBooleanValue(TEST_VAR_SET);
    TEUCHOS_TEST_ASSERT(envar_is_true, out, success);
  }

  TEUCHOS_UNIT_TEST(Environment, Check_getBooleanValue_2) {
    out << "Check nonexistent variable returns false" << endl;
    bool envar_is_false =
      Environment::getInstance().getBooleanValue(TEST_VAR_NOT_SET);
    TEUCHOS_TEST_ASSERT(!envar_is_false, out, success);
  }

  // TEUCHOS_UNIT_TEST(Environment, Check_setValue) {
  //   out << "Check variable setting" << endl;
  //   std::string envar;
  //   std::string value1 = "X1";
  //   std::string value2 = "X2";

  //   out << "Make sure nonexistent envar does not exist" << endl;
  //   bool env_does_not_exist =
  //     ! Environment::getInstance().variableExists(TEST_VAR_NOT_SET);
  //   TEUCHOS_TEST_ASSERT(env_does_not_exist, out, success);

  //   out << "Now set the value" << endl;
  //   // setValue overwrites by default
  //   Environment::getInstance().setValue(TEST_VAR_NOT_SET, value1);

  //   out << "Check that the value was set correctly" << endl;
  //   envar = Environment::getInstance().getValue(TEST_VAR_NOT_SET);
  //   bool envar_set_correctly = envar == value1;
  //   TEUCHOS_TEST_ASSERT(envar_set_correctly, out, success);

  //   out << "Try to set the value again, but with overwrite=0" << endl;
  //   int overwrite = 0;
  //   Environment::getInstance().setValue(TEST_VAR_NOT_SET, value2, overwrite);

  //   out << "Check that the value was NOT set" << endl;
  //   envar = Environment::getInstance().getValue(TEST_VAR_NOT_SET);
  //   bool envar_not_set = envar != value2;
  //   TEUCHOS_TEST_ASSERT(envar_not_set, out, success);

  //   out << "Check that the value is still correct" << endl;
  //   bool envar_still_correct = envar == value1;
  //   TEUCHOS_TEST_ASSERT(envar_still_correct, out, success);

  //   // Unset the environment variable
  //   unsetenv(TEST_VAR_NOT_SET.c_str());
  // }

} // namespace (anonymous)
