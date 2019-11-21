//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Tsqr_KokkosNodeTsqrTest.hpp"
#include "Kokkos_Core.hpp"

#ifdef HAVE_KOKKOSTSQR_COMPLEX
#  include <complex>
#endif // HAVE_KOKKOSTSQR_COMPLEX

namespace {
  //
  // The documentation string for this test executable to print out at
  // the command line on request.
  //
  const char docString[] = "This program tests TSQR::KokkosNodeTsqr, "
    "which implements an intranode parallel version of TSQR for "
    "Kokkos::DefaultHostExecutionSpace.  Accuracy and performance "
    "tests are included.";

  //
  // TestParameters encapsulates values of command-line parameters, as
  // well as state that may change from one benchmark / verify
  // invocation to the next.
  //
  class TestParameters {
  public:
    TestParameters () = default;
    TestParameters (const std::vector<int> /* theSeed */);

    bool verify = true;
    bool benchmark = false;
    int numRows = 100000;
    int numCols = 10;
    int numTrials = 1;
    bool testReal = true;
#ifdef HAVE_KOKKOSTSQR_COMPLEX
    bool testComplex = true;
#endif // HAVE_KOKKOSTSQR_COMPLEX
    int numPartitions = 16;
    int cacheSizeHint = 0;
    bool contiguousCacheBlocks = false;
    bool printFieldNames = true;
    bool humanReadable = true;
    bool debug = false;
  };

  // Run the test(s) for a particular Scalar type T.
  // Used by Cons, which in turn is used by runTests().
  template<class T>
  class Dispatcher {
  public:
    typedef T dispatch_type;

    static void
    benchmark (std::vector<int>&,
               const TestParameters& params,
               bool& printFieldNames)
    {
      using TSQR::Test::benchmarkKokkosNodeTsqr;
      benchmarkKokkosNodeTsqr<int, T> (params.numTrials,
                                       params.numRows,
                                       params.numCols,
                                       params.numPartitions,
                                       params.cacheSizeHint,
                                       params.contiguousCacheBlocks,
                                       printFieldNames,
                                       params.humanReadable);
      printFieldNames = false;
    }

    static void
    verify (std::vector<int>& seed,
            const TestParameters& params,
            bool& printFieldNames)
    {
      TSQR::Random::NormalGenerator<int, T> gen (seed);
      using TSQR::Test::verifyKokkosNodeTsqr;
      verifyKokkosNodeTsqr<int, T> (gen,
                                    params.numRows,
                                    params.numCols,
                                    params.numPartitions,
                                    params.cacheSizeHint,
                                    params.contiguousCacheBlocks,
                                    printFieldNames,
                                    params.humanReadable,
                                    params.debug);
      printFieldNames = false;
      // Save the seed for next time, since we can't use the same
      // NormalGenerator for a different Scalar type T.
      gen.getSeed (seed);
    }
  };

  //
  // Class for executing a template function over a compile-time
  // fixed-length list of types.  See runTests() for an example.
  //
  template<class CarType, class CdrType>
  class Cons {
  public:
    static void
    verify (std::vector<int>& seed,
            const TestParameters& params,
            bool& printFieldNames)
    {
      Dispatcher<CarType>::verify (seed, params, printFieldNames);
      CdrType::verify (seed, params, printFieldNames);
    }

    static void
    benchmark (std::vector<int>& seed,
               const TestParameters& params,
               bool& printFieldNames)
    {
      Dispatcher<CarType>::benchmark (seed, params, printFieldNames);
      CdrType::benchmark (seed, params, printFieldNames);
    }
  };

  // Base case for Cons template recursion.
  class NullCons {
  public:
    static void
    verify (std::vector<int>&,
            const TestParameters&,
            bool& printFieldNames) {}

    static void
    benchmark (std::vector<int>&,
               const TestParameters&,
               bool& printFieldNames) {}
  };

  // Run the tests for all types of interest.
  void
  runTests (const TestParameters& params)
  {
    using real_tests = Cons<float, Cons<double, NullCons>>;
#ifdef HAVE_KOKKOSTSQR_COMPLEX
    using complex_tests =
      Cons<std::complex<float>, Cons<std::complex<double>, NullCons>>;
#endif // HAVE_KOKKOSTSQR_COMPLEX

    // Length-4 seed for the pseudorandom number generator.  The last
    // entry must be an odd number.  There are other restrictions on
    // these values; see the LAPACK documentation for details.  (0, 0,
    // 0, 1) is a typical initial seed if you want reproducible
    // results, but don't actually care much about randomness.
    std::vector<int> seed {{0, 0, 0, 1}};

    bool printFieldNames = params.printFieldNames;
    if (params.verify) {
      if (params.testReal) {
        real_tests::verify (seed, params, printFieldNames);
      }
#ifdef HAVE_KOKKOSTSQR_COMPLEX
      if (params.testComplex) {
        complex_tests::verify (seed, params, printFieldNames);
      }
#endif // HAVE_KOKKOSTSQR_COMPLEX
    }
    // Reset this, since the first call of verify() sets it to false.
    printFieldNames = params.printFieldNames;
    if (params.benchmark) {
      if (params.testReal) {
        real_tests::benchmark (seed, params, printFieldNames);
      }
#ifdef HAVE_KOKKOSTSQR_COMPLEX
      if (params.testComplex) {
        complex_tests::benchmark (seed, params, printFieldNames);
      }
#endif // HAVE_KOKKOSTSQR_COMPLEX
    }
  }

  // Parse command-line options for this test.
  //
  // argc [in] As usual in C(++)
  //
  // argv [in] As usual in C(++)
  //
  // allowedToPrint [in] Whether this (MPI) process is allowed
  //   to print to stdout/stderr.  Different per (MPI) process.
  //
  // printedHelp [out] Whether this (MPI) process printed the
  //   "help" display (summary of command-line options).
  //
  // Return an encapsulation of the command-line options.
  TestParameters
  parseOptions (int argc,
                char* argv[],
                const bool allowedToPrint,
                bool& printedHelp)
  {
    using std::cerr;
    using std::endl;

    printedHelp = false;

    // Command-line parameters, set to their default values.
    TestParameters params;
    /// We really want the cache size hint as a size_t, but
    /// Teuchos::CommandLineProcessor doesn't offer that option.  So
    /// we read it in as an int, which means negative inputs are
    /// possible.  We check for those below in the input validation
    /// phase.
    //
    // Fetch default value of cacheSizeHint.
    int cacheSizeHint = params.cacheSizeHint;
    try {
      using Teuchos::CommandLineProcessor;

      CommandLineProcessor cmdLineProc (/* throwExceptions=*/ true,
                                        /* recognizeAllOptions=*/ true);
      cmdLineProc.setDocString (docString);
      cmdLineProc.setOption ("verify",
                             "noverify",
                             &params.verify,
                             "Test accuracy");
      cmdLineProc.setOption ("benchmark",
                             "nobenchmark",
                             &params.benchmark,
                             "Test performance");
      cmdLineProc.setOption ("numRows",
                             &params.numRows,
                             "Number of rows in the test matrix");
      cmdLineProc.setOption ("numCols",
                             &params.numCols,
                             "Number of columns in the test matrix");
      cmdLineProc.setOption ("numTrials",
                             &params.numTrials,
                             "Number of trials (only used when \"--benchmark\"");
      cmdLineProc.setOption ("testReal",
                             "noTestReal",
                             &params.testReal,
                             "Test real arithmetic");
#ifdef HAVE_KOKKOSTSQR_COMPLEX
      cmdLineProc.setOption ("testComplex",
                             "noTestComplex",
                             &params.testComplex,
                             "Test complex arithmetic");
#endif // HAVE_KOKKOSTSQR_COMPLEX
      params.numPartitions = Kokkos::DefaultHostExecutionSpace::concurrency();
      cmdLineProc.setOption ("numPartitions",
                             &params.numPartitions,
                             "Number of partitions to use (max available parallelism)");
      cmdLineProc.setOption ("cacheSizeHint",
                             &cacheSizeHint,
                             "Cache size hint in bytes (0 means pick a reasonable default)");
      cmdLineProc.setOption ("contiguousCacheBlocks",
                             "noncontiguousCacheBlocks",
                             &params.contiguousCacheBlocks,
                             "Whether cache blocks should be stored contiguously");
      cmdLineProc.setOption ("printFieldNames",
                             "noPrintFieldNames",
                             &params.printFieldNames,
                             "Print field names (for machine-readable output only)");
      cmdLineProc.setOption ("humanReadable",
                             "machineReadable",
                             &params.humanReadable,
                             "If set, make output easy to read by humans "
                             "(but hard to parse)");
      cmdLineProc.setOption ("debug",
                             "noDebug",
                             &params.debug,
                             "Print debugging information");
      cmdLineProc.parse (argc, argv);
    }
    catch (Teuchos::CommandLineProcessor::UnrecognizedOption& e) {
      if (allowedToPrint)
        cerr << "Unrecognized command-line option: " << e.what() << endl;
      throw e;
    }
    catch (Teuchos::CommandLineProcessor::HelpPrinted& e) {
      printedHelp = true;
      return params; // Don't verify parameters in this case
    }

    // Validate command-line options.  We provide default values
    // for unset options, so we don't have to validate those.
    if (params.numRows <= 0) {
      throw std::invalid_argument ("Number of rows must be positive");
    } else if (params.numCols <= 0) {
      throw std::invalid_argument ("Number of columns must be positive");
    } else if (params.numRows < params.numCols) {
      throw std::invalid_argument ("Number of rows must be >= number of columns");
    } else if (params.benchmark && params.numTrials < 1) {
      throw std::invalid_argument ("\"--benchmark\" option requires numTrials >= 1");
    } else if (params.numPartitions < 1) {
      throw std::invalid_argument ("\"--numPartitions\" option must be >= 1");
    } else if (params.cacheSizeHint < 0) {
      throw std::invalid_argument ("Cache size hint must be nonnegative");
    }
    return params;
  }
} // namespace (anonymous)

//
// The "main" test driver.
//
int
main (int argc, char *argv[])
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  bool performingTests = true;
  const bool allowedToPrint = true;
  std::ostream& out = std::cout;

  // Fetch command-line parameters.
  bool printedHelp = false;
  TestParameters params =
    parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp) {
    return EXIT_SUCCESS;
  }

  bool success = false;
  bool verbose = false;
  try {
    if (performingTests) {
      Kokkos::ScopeGuard kokkosScope (argc, argv);
      runTests (params);
      success = true;
      // The Trilinos test framework expects a message like this.
      out << "\nEnd Result: TEST PASSED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
