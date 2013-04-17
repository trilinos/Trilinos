// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
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
  bool *success,
  std::ostream *out
  )
{
  if (!result) *success = false;
  if (out) {
    if( !result || show_all_tests ) {
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
    Thyra::printTestResults(this_result, ossStore_.str(), false, &*success, &*out_);
  }
  else {
    if (!this_result) {
      *success = false;
    }
  }
  printedTestResults_ = true;
}


} // namespace Thyra
