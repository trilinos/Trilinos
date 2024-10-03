// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include <algorithm> // std::transform
#include <cstdlib> // std::getenv
#include <cctype> // std::toupper
#include <string>

namespace { // (anonymous)

  //
  // UNIT TESTS
  //

#ifdef HAVE_TPETRA_DEBUG
  constexpr bool defaultDebug = true;
#else
  constexpr bool defaultDebug = false;
#endif // HAVE_TPETRA_DEBUG

  constexpr bool defaultVerbose = false;

  constexpr size_t defaultThreshold = 200;

  // See example here:
  //
  // http://en.cppreference.com/w/cpp/string/byte/toupper
  std::string stringToUpper (std::string s)
  {
    std::transform (s.begin (), s.end (), s.begin (),
                    [] (unsigned char c) { return std::toupper (c); });
    return s;
  }

  bool
  getEnvironmentVariableAsBool (Teuchos::FancyOStream& out,
                                const char environmentVariableName[],
                                const bool defaultValue)
  {
    out << "Attempt to get value of environment variable \""
        << environmentVariableName << "\"" << std::endl;
    Teuchos::OSTab tab1 (out);

    const char* varVal = std::getenv (environmentVariableName);

    bool retVal = defaultValue;
    if (varVal != NULL) {
      out << "Got value \"" << varVal << "\"" << std::endl;
      const std::string varStr (stringToUpper (std::string (varVal)));
      out << "Upper-case value: \"" << varStr << "\"" << std::endl;
      if (varStr == "1" || varStr == "YES" || varStr == "TRUE") {
        out << "Value is true" << std::endl;
        retVal = true;
      }
      else if (varStr == "0" || varStr == "NO" || varStr == "FALSE") {
        out << "Value is false" << std::endl;
        retVal = false;
      }
      else {
        out << "Value is neither true nor false, as far as I can tell"
            << std::endl;
      }
      // Otherwise, use the default value.
    }
    else {
      out << "getenv returned NULL; use default value "
          << (defaultValue ? "true" : "false") << std::endl;
    }
    return retVal;
  }

  TEUCHOS_UNIT_TEST( TpetraDetailsBehavior, DebugAndVerbose )
  {
    using Tpetra::Details::Behavior;
    using std::endl;

    out << "Test Tpetra::Details::Behavior" << endl;
    Teuchos::OSTab tab1 (out);

    {
      out << "Test Behavior::debug()" << endl;
      Teuchos::OSTab tab2 (out);

      const char* varVal = std::getenv ("TPETRA_DEBUG");
      if (varVal == NULL) {
        out << "TPETRA_DEBUG is not already set in environment" << endl;
        Teuchos::OSTab tab3 (out);

        const bool debug = Behavior::debug ();
        TEST_EQUALITY( debug, defaultDebug );
      }
      else {
        out << "TPETRA_DEBUG is already set in environment" << endl;
        Teuchos::OSTab tab3 (out);

        const bool expectedDebug =
          getEnvironmentVariableAsBool (out, "TPETRA_DEBUG", defaultDebug);
        const bool debug = Behavior::debug ();
        TEST_EQUALITY( debug, expectedDebug );
      }
    }

    {
      out << "Test Behavior::verbose()" << endl;
      Teuchos::OSTab tab2 (out);

      const char* varVal = std::getenv ("TPETRA_VERBOSE");
      if (varVal == NULL) {
        out << "TPETRA_VERBOSE is not already set in environment" << endl;
        Teuchos::OSTab tab3 (out);

        const bool verbose = Behavior::verbose ();
        TEST_EQUALITY( verbose, defaultVerbose );
      }
      else {
        out << "TPETRA_VERBOSE is already set in environment" << endl;
        Teuchos::OSTab tab3 (out);

        const bool expectedVerbose =
          getEnvironmentVariableAsBool (out, "TPETRA_VERBOSE", defaultVerbose);
        const bool verbose = Behavior::verbose ();
        TEST_EQUALITY( verbose, expectedVerbose );
      }
    }
  }

  TEUCHOS_UNIT_TEST( TpetraDetailsBehavior, verbosePrintCountThreshold )
  {
    using Tpetra::Details::Behavior;
    using std::endl;

    out << "Test Tpetra::Details::Behavior" << endl;
    Teuchos::OSTab tab1 (out);

    {
      out << "Test Behavior::verbosePrintCountThreshold()" << endl;
      Teuchos::OSTab tab2 (out);

      const char* varVal = std::getenv ("TPETRA_VERBOSE_PRINT_COUNT_THRESHOLD");
      if (varVal == NULL) {
        out << "TPETRA_VERBOSE_PRINT_COUNT_THRESHOLD is not already set in environment" << endl;
        Teuchos::OSTab tab3 (out);

        const size_t threshold = Behavior::verbosePrintCountThreshold ();
        TEST_EQUALITY( threshold, defaultThreshold );
      }
      else {
        out << "TPETRA_VERBOSE_PRINT_COUNT_THRESHOLD is already set in environment" << endl;
        Teuchos::OSTab tab3 (out);

        size_t threshold = std::stoi(std::string(varVal));
        size_t threshold2 = Behavior::verbosePrintCountThreshold ();
        TEST_EQUALITY(threshold, threshold2);
      }
    }

  } //verbosePrintCountThreshold

} // namespace (anonymous)


