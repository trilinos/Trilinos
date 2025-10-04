// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_TestingTools.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace Thyra {

//
// Helper code and declarations
//

using Teuchos::as;
using Teuchos::fancyOStream;
using Teuchos::get_extra_data;
using Teuchos::inoutArg;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::Utils;

//
// Unit Tests
//

TEUCHOS_UNIT_TEST(TestResultsPrinter, show_all_tests_pass) {
  out << "*** Testing that output is send directly to main output stream!\n";
  ECHO(const bool show_all_tests = true);
  ECHO(std::ostringstream myOut);
  ECHO(std::ostringstream myOut2);
  {
    TestResultsPrinter testResultsPrinter(fancyOStream(rcpFromRef(myOut)), show_all_tests);
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();
    ECHO(bool this_result = true);
    ECHO(int num = 5);
    TEUCHOS_TEST_EQUALITY_CONST(num, 5, *testOut, this_result);
    out << "*** Make sure test output is being directly pritned to myOut!\n";
    const std::string myOut_str = Utils::trimWhiteSpace(myOut.str());
    out << "myOut.str() = '" << myOut_str << "'\n";
    TEST_EQUALITY_CONST(myOut_str, std::string("num = 5 == 5 : passed"));
    out << "*** Make sure that we don't see any addition output printed.\n";
    ECHO(std::ostringstream finalOut);
    ECHO(testResultsPrinter.replaceOStream(fancyOStream(rcpFromRef(finalOut))));
    ECHO(bool final_success = true);
    ECHO(testResultsPrinter.printTestResults(this_result, inoutArg(final_success)));
    TEST_EQUALITY_CONST(final_success, true);
    TEST_EQUALITY_CONST(finalOut.str(), "");
    out << "*** Make sure that nothing gets printed out from the destructor\n";
    ECHO(testResultsPrinter.replaceOStream(fancyOStream(rcpFromRef(myOut2))));
  }
  TEST_EQUALITY_CONST(myOut2.str(), std::string(""));
}

TEUCHOS_UNIT_TEST(TestResultsPrinter, show_all_tests_fail) {
  out << "*** Testing that output is send directly to main output stream!\n";
  ECHO(const bool show_all_tests = true);
  ECHO(std::ostringstream myOut);
  ECHO(std::ostringstream myOut2);
  {
    TestResultsPrinter testResultsPrinter(
        fancyOStream(rcpFromRef(myOut)), show_all_tests);
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();
    ECHO(bool this_result = true);
    ECHO(int num = 5);
    TEUCHOS_TEST_EQUALITY_CONST(num, 4, *testOut, this_result);
    TEST_EQUALITY_CONST(this_result, false);
    out << "*** Make sure test output is being directly pritned to myOut!\n";
    const std::string myOut_str = Utils::trimWhiteSpace(myOut.str());
    out << "myOut.str() = '" << myOut_str << "'\n";
    // NOTE: Tetsing myOut.str() is hard becuase it has file path in failed msg.
    out << "*** Make sure that we don't see any addition output printed.\n";
    ECHO(std::ostringstream finalOut);
    ECHO(testResultsPrinter.replaceOStream(fancyOStream(rcpFromRef(finalOut))));
    ECHO(bool final_success = true);
    testResultsPrinter.printTestResults(this_result, inoutArg(final_success));
    TEST_EQUALITY_CONST(final_success, false);
    TEST_EQUALITY_CONST(finalOut.str(), "");
    out << "*** Make sure that nothing gets printed out from the destructor\n";
    ECHO(testResultsPrinter.replaceOStream(fancyOStream(rcpFromRef(myOut2))));
  }
  TEST_EQUALITY_CONST(myOut2.str(), std::string(""));
}

TEUCHOS_UNIT_TEST(TestResultsPrinter, no_show_all_tests_pass) {
  out << "*** Testing that output is send directly to main output stream!\n";
  ECHO(const bool show_all_tests = false);
  ECHO(std::ostringstream myOut);
  ECHO(std::ostringstream myOut2);
  {
    TestResultsPrinter testResultsPrinter(fancyOStream(rcpFromRef(myOut)), show_all_tests);
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();
    ECHO(bool this_result = true);
    ECHO(int num = 5);
    TEUCHOS_TEST_EQUALITY_CONST(num, 5, *testOut, this_result);
    out << "*** Make sure test output is *not* directly printing myOut!\n";
    TEST_EQUALITY_CONST(myOut.str(), std::string(""));
    out << "*** Make sure that the only thing printed is 'passed'.\n";
    ECHO(bool final_success = true);
    ECHO(testResultsPrinter.printTestResults(this_result, inoutArg(final_success)));
    TEST_EQUALITY_CONST(final_success, true);
    TEST_EQUALITY_CONST(myOut.str(), "passed!\n");
    out << "*** Make sure that nothing gets printed out from the destructor\n";
    ECHO(testResultsPrinter.replaceOStream(fancyOStream(rcpFromRef(myOut2))));
  }
  TEST_EQUALITY_CONST(myOut2.str(), std::string(""));
}

TEUCHOS_UNIT_TEST(TestResultsPrinter, no_show_all_tests_fail) {
  out << "*** Testing that output is send directly to main output stream!\n";
  ECHO(const bool show_all_tests = false);
  ECHO(std::ostringstream myOut);
  ECHO(std::ostringstream myOut2);
  {
    TestResultsPrinter testResultsPrinter(fancyOStream(rcpFromRef(myOut)), show_all_tests);
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();
    ECHO(bool this_result = true);
    ECHO(int num = 5);
    TEUCHOS_TEST_EQUALITY_CONST(num, 4, *testOut, this_result);
    out << "*** Make sure test output is *not* directly printing myOut!\n";
    TEST_EQUALITY_CONST(myOut.str(), std::string(""));
    out << "*** Make sure that the test results details are printed because test failed\n";
    ECHO(bool final_success = true);
    ECHO(testResultsPrinter.printTestResults(this_result, inoutArg(final_success)));
    TEST_EQUALITY_CONST(final_success, false);
    const std::string myOut_str = Utils::trimWhiteSpace(myOut.str());
    TEST_EQUALITY_CONST(myOut_str.substr(0, 25), "num = 5 == 4 : FAILED ==>");
    out << "*** Make sure that nothing gets printed out from the destructor\n";
    ECHO(testResultsPrinter.replaceOStream(fancyOStream(rcpFromRef(myOut2))));
  }
  TEST_EQUALITY_CONST(myOut2.str(), std::string(""));
}

TEUCHOS_UNIT_TEST(TestResultsPrinter, no_show_all_tests_pass_throws) {
  out << "*** Testing that output is send directly to main output stream!\n";
  ECHO(const bool show_all_tests = false);
  ECHO(std::ostringstream myOut);
  try {
    TestResultsPrinter testResultsPrinter(fancyOStream(rcpFromRef(myOut)), show_all_tests);
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();
    ECHO(bool this_result = true);
    ECHO(int num = 5);
    TEUCHOS_TEST_EQUALITY_CONST(num, 5, *testOut, this_result);
    TEST_EQUALITY_CONST(this_result, true);
    out << "*** Make sure test output is *not* directly printing myOut!\n";
    TEST_EQUALITY_CONST(myOut.str(), std::string(""));
    out << "*** Throw before we can print the final test results!\n";
    TEUCHOS_TEST_FOR_EXCEPTION(!(num == 10), std::logic_error,
                               "Unexpected exception from my test!");
    TEST_ASSERT(false);  // If we get here, we failed!
  } catch (const std::logic_error &except) {
    out << "*** Check to make sure test result got printed in destructor!\n";
    const std::string myOut_str = Utils::trimWhiteSpace(myOut.str());
    TEST_EQUALITY_CONST(myOut_str, std::string("num = 5 == 5 : passed"));
    out << "Caught except.what() = '" << except.what() << "'\n";
  }
  out << "*** Make double sure exception got thrown and printed out myOut!\n";
  TEST_ASSERT(myOut.str().length() > 0);
}

}  // namespace Thyra
