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

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <sstream>

namespace {

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP ();
    clp.addOutputSetupOptions (true);
  }

  // "Local string" used in the test below.
  std::string localString (const int myRank, const int numProcs) {
    std::ostringstream os;
    os << "Hello from Process " << myRank << "!";
    if (myRank + 1 < numProcs) {
      os << "  ";
    }
    return os.str ();
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( TpetraUtils, gathervPrint )
  {
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    using std::endl;

    int lclSuccess = 1;
    int gblSuccess = 1;

    auto comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    // On each process, construct the local string, which we will
    // gather.
    std::string lclString = localString (myRank, numProcs);
    const size_t lclSize = lclString.size ();

    // Compute expected total output length.
    size_t gblSize = 0;
    reduceAll<int, size_t> (*comm, REDUCE_SUM, lclSize, outArg (gblSize));

    std::ostringstream gatheredString;
    try {
      Tpetra::Details::gathervPrint (gatheredString, lclString, *comm);
    }
    catch (std::exception& e) {
      lclSuccess = 0;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "FAILED: gathervPrint threw an exception on at least one process!" << endl;
      success = false; // just to be sure
      return;
    }

    std::string gblString = gatheredString.str ();
    if (myRank == 0) {
      TEST_EQUALITY( gblString.size (), gblSize );
      if (gblString.size () != gblSize) {
        lclSuccess = 0;
      }
    }
    else {
      // Any other process should not have written to the output stream.
      TEST_EQUALITY( gblString.size (), static_cast<size_t> (0) );
      if (gblString.size () != static_cast<size_t> (0)) {
        lclSuccess = 0;
      }
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "FAILED to get correct string length on at least one process!" << endl;
      success = false; // just to be sure
      return;
    }

    // On Process 0, check whether the received string is correct.
    if (myRank == 0) {
      out << "Result of gathervPrint: " << endl << gblString << endl;

      // Construct the comparison string, that we will use to test the
      // gathered string.
      std::ostringstream comparisonStream;
      for (int p = 0; p < numProcs; ++p) {
        comparisonStream << localString (p, numProcs);
      }
      std::string comparisonString = comparisonStream.str ();
      const size_t comparisonSize = comparisonString.size ();

      out << "Expected result: " << endl << comparisonString << endl;

      TEST_EQUALITY( comparisonSize, gblSize );
      if (comparisonSize != gblSize) {
        lclSuccess = 0;
      }

      const bool stringsSame = (comparisonString == gblString);
      TEST_EQUALITY_CONST( stringsSame, true );
      if (! stringsSame) {
        lclSuccess = 0;
      }
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "FAILED: Gathered string did not match expected string!" << endl;
      success = false; // just to be sure
      return;
    }
  }

} // namespace (anonymous)


