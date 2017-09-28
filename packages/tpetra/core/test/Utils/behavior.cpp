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

} // namespace (anonymous)


