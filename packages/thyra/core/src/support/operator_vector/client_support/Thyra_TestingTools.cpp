// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_TestingTools.hpp"


bool Thyra::testBoolExpr(
  const std::string &boolExprName,
  const bool &boolExpr,
  const bool &boolExpected,
  const Ptr<std::ostream> &out,
  const std::string &li
  )
{
  const bool success = ( boolExpr == boolExpected );
  if (nonnull(out)) {
    *out
      << std::endl
      << li << "Check: " << boolExprName << " = " << boolExpr << " == " << boolExpected
      << " : " << passfail(success) << std::endl;
  }
  return success;
}


void Thyra::printTestResults(
  const bool result,
  const std::string &test_summary,
  const bool show_all_tests,
  const Ptr<bool> &success,
  const Ptr<std::ostream> &out
  )
{
  if (!result) *success = false;
  if (nonnull(out)) {
    if (!result || show_all_tests) {
      *out << std::endl << test_summary;
    }
    else {
      *out << "passed!\n";
    }
  }
}


// TestResultsPrinter


namespace Thyra {


TestResultsPrinter::TestResultsPrinter(
  const RCP<FancyOStream> &out, const bool show_all_tests)
  : out_(out.assert_not_null()), show_all_tests_(show_all_tests),
    printedTestResults_(false)
{
  if (show_all_tests_) {
    oss_ = out_;
  }
  else {
    oss_ = Teuchos::fancyOStream(Teuchos::rcpFromRef(ossStore_));
    ossStore_.copyfmt(*out_);
  }
}


TestResultsPrinter::~TestResultsPrinter()
{
  using Teuchos::inoutArg;
  if (!printedTestResults_) {
    // If we get here, either someone made a mistake in not calling
    // printTestResults() or an exception was thrown.  Either way, we are
    // going to assume failure and dump everything.
    try {
      bool dummy_success = true;
      this->printTestResults(false, inoutArg(dummy_success));
    }
    catch(...) {
      // Need to eat any exceptions in case an exception is already active
      // which is calling this destructor.
    }
  }
  printedTestResults_ = true;
}


RCP<FancyOStream>
TestResultsPrinter::replaceOStream(const RCP<FancyOStream> &out)
{
  const RCP<FancyOStream> oldOut = out_;
  out_ = out;
  return oldOut;
}


RCP<FancyOStream> TestResultsPrinter::getTestOStream()
{
  return oss_;
}


void TestResultsPrinter::printTestResults(const bool this_result,
  const Ptr<bool> &success)
{
  if (!show_all_tests_) {
    Thyra::printTestResults(this_result, ossStore_.str(), false,
      success, out_.ptr());
  }
  else {
    if (!this_result) {
      *success = false;
    }
  }
  printedTestResults_ = true;
}


} // namespace Thyra
