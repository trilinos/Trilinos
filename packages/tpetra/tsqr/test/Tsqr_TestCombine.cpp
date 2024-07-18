// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Time.hpp"
#include "Tsqr_CombineBenchmark.hpp"
#include "Tsqr_CombineTest.hpp"

#ifdef HAVE_TPETRATSQR_COMPLEX
#  include <complex>
#endif // HAVE_TPETRATSQR_COMPLEX

#include "Kokkos_Core.hpp"
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
  using Teuchos::RCP;

  //
  // Short description of this program, printed if the --help option
  // is given at the command line.
  //
  static const char docString[] =
    "This program tests accuracy and performance of TSQR::Combine.";

  //
  // The TestParameters struct encapsulates values of command-line
  // parameters.
  //
  struct TestParameters {
    // Whether to run the accuracy test.
    bool verify = true;
    // Whether to run the performance test.
    bool benchmark = false;
    // Number of rows in the test matrix.
    int numRows = 100;
    // Number of columns in the test matrix.
    int numCols = 5;
    // Number of trials (benchmark only).
    int numTrials = 3;
    // Whether to pick the number of trials automatically, using an
    // iterative calibration process (benchmark only).
    bool calibrate = false;
    // Whether to print averaged timings over all trials (true), or
    // the cumulative timing over all trials (false).
    bool averageTimings = true;
    // Whether to test real-arithmetic routines.
    bool testReal = true;
    // Whether to test complex-arithmetic routines.  If TSQR was not
    // built with complex arithmetic support, then this must always be
    // false.
#ifdef HAVE_TPETRATSQR_COMPLEX
    bool testComplex = true;
#else
    bool testComplex = false;
#endif // HAVE_TPETRATSQR_COMPLEX
    // Whether to print column (field) names.
    bool printFieldNames = true;
    // Whether to print output that the Trilinos test framework
    // expects, in order to judge a test as passed or failed.
    bool printTrilinosTestStuff = true;
    // Whether the benchmark should fail if performance of
    // TSQR::CombineNative relative to that of TSQR::CombineDefault is
    // not good enough.
    bool strictPerfTests = false;
    // If strictPerfTests is true: how much slower CombineNative (and
    // CombineFortran, if applicable) is allowed to be, relative to
    // CombineDefault.
    double allowance = 1.2;
    // Whether to print verbose status output.
    bool verbose = true;
    // Whether to print debugging output to stderr.
    bool debug = false;

    std::string additionalFieldNames;
    std::string additionalData;
  };

  // Benchmark TSQR::Combine.
  //
  // out [out] output stream for benchmark results.
  //   It will only be used on rank 0.
  //
  // params [in] test parameter struct.  This method reads
  //   the following fields: numRows, numCols, numTrials,
  //   testReal, testComplex.
  //
  // Warning: Call only on (MPI) Process 0.  Otherwise, you'll run the
  //   test routine on every MPI process simultaneously, but only
  //   report results on Process 0.
  void
  benchmark(std::ostream& out,
            const TestParameters& params)
  {
    std::vector<int> seed(4);
    const bool useSeedValues = false; // Fill in seed with defaults.

    TSQR::Test::CombineBenchmarkParameters testParams;
    testParams.numRows = params.numRows;
    testParams.numCols = params.numCols;
    testParams.testReal = params.testReal;
    testParams.testComplex = params.testComplex;
    testParams.numTrials = params.numTrials;
    testParams.calibrate = params.calibrate;
    testParams.averageTimings = params.averageTimings;
    testParams.strictPerfTests = params.strictPerfTests;
    testParams.allowance = params.allowance;
    testParams.seed = seed;
    testParams.useSeedValues = useSeedValues;
    testParams.additionalFieldNames = params.additionalFieldNames;
    testParams.additionalData = params.additionalData;
    testParams.printFieldNames = params.printFieldNames;
    testParams.debug = params.debug;

    using timer_type = Teuchos::Time;
    TSQR::Test::benchmarkCombine<timer_type>(out, testParams);
  }

  // Test accuracy of TSQR::Combine.
  //
  // out [out] output stream for benchmark results.  It will only be
  //   used on Process 0.
  //
  // params [in] test parameter struct.  This method reads the
  //   following fields: numRows, numCols, numTrials, testReal,
  //   testComplex.
  //
  // Warning: Call only on (MPI) Process 0.  Otherwise, you'll run the
  //   test routine on every MPI process simultaneously, but only
  //   report results on Process 0.
  void
  verify(std::ostream& out, const TestParameters& params)
  {
    constexpr bool simulateSequentialTsqr = false;
    constexpr bool debug = false;

    using TSQR::Test::verifyCombine;
    verifyCombine(params.numRows, params.numCols, params.testReal,
                  params.testComplex, params.printFieldNames,
                  simulateSequentialTsqr, debug);
  }

  // \brief Parse command-line options for this test
  //
  // argc [in] As usual in C(++).
  // argv [in] As usual in C(++).
  //
  // allowedToPrint [in] Whether this (MPI) process is allowed
  //   to print to stdout/stderr.  Different per (MPI) process.
  //
  // printedHelp [out] Whether this (MPI) process printed the
  //   "help" display (summary of command-line options)
  //
  // Return: Encapsulation of command-line options.
  TestParameters
  parseOptions(int argc,
               char* argv[],
               std::ostream& err,
               bool& printedHelp)
  {
    using std::endl;

    printedHelp = false;

    // Command-line parameters, set to their default values.
    TestParameters params {};
    try {
      constexpr bool throwExceptions = true;
      constexpr bool recognizeAllOptions = true;
      using CLP = Teuchos::CommandLineProcessor;
      CLP cmdLineProc(throwExceptions, recognizeAllOptions);
      cmdLineProc.setDocString(docString);
      cmdLineProc.setOption("verify",
                            "noverify",
                            &params.verify,
                            "Test accuracy of TSQR::Combine implementations.");
      cmdLineProc.setOption("benchmark",
                            "nobenchmark",
                            &params.benchmark,
                            "Test performance of TSQR::Combine implementations.");
      cmdLineProc.setOption("debug",
                            "nodebug",
                            &params.debug,
                            "Print copious debugging information to stderr.");
      cmdLineProc.setOption("numRows",
                            &params.numRows,
                            "Number of rows in the cache block test.");
      cmdLineProc.setOption("numCols",
                            &params.numCols,
                            "Number of columns in the cache block test, and "
                            "number of rows and columns in each upper triangular "
                            "matrix in the pair test.");
      cmdLineProc.setOption("numTrials",
                            &params.numTrials,
                            "For benchmarks: Number of trials.  "
                            "Ignored if --calibrate option is set.");
      cmdLineProc.setOption("calibrate",
                            "noCalibrate",
                            &params.calibrate,
                            "For benchmarks: ignore numTrials, and calibrate "
                            "the number of trials based on computed timer "
                            "resolution and problem size (numRows and "
                            "numCols).");
      cmdLineProc.setOption("meanTimings",
                            "sumTimings",
                            &params.averageTimings,
                            "For benchmarks: whether timings should be "
                            "computed as an arithmetic mean (true) or as a "
                            "sum (false) over all trials.");
      cmdLineProc.setOption("testReal",
                            "noTestReal",
                            &params.testReal,
                            "Test real-arithmetic routines.");
      cmdLineProc.setOption("testComplex",
                            "noTestComplex",
                            &params.testComplex,
                            "Test complex-arithmetic routines.  This option "
                            "may only be true if Trilinos was built with "
                            "complex arithmetic support.");
      cmdLineProc.setOption("strictPerfTests",
                            "noStrictPerfTests",
                            &params.strictPerfTests,
                            "For benchmarks: whether the test should fail if "
                            "run time of TSQR::CombineNative / run time of "
                            "TSQR::CombineDefault (both for the cache block "
                            "benchmark) is greater than the given slowdown "
                            "allowance.  Ditto for TSQR::CombineFortran, if "
                            "TSQR was built with Fortran support.");
      cmdLineProc.setOption("allowance",
                            &params.allowance,
                            "For benchmarks: if strictPerfTests is true: "
                            "allowed slowdown factor.  If exceeded, the test "
                            "fails.");
      cmdLineProc.setOption("additionalFieldNames",
                            &params.additionalFieldNames,
                            "Any additional field name(s) (comma-delimited "
                            "string) to add to the benchmark output.  Empty "
                            "by default.  Good for things known when invoking "
                            "the benchmark executable, but not (easily) known "
                            "inside the benchmark -- e.g., environment "
                            "variables.");
      cmdLineProc.setOption("additionalData",
                            &params.additionalData,
                            "Any additional data to add to the output, "
                            "corresponding to the above field name(s). "
                            "Empty by default.");
      cmdLineProc.setOption("printFieldNames",
                            "noPrintFieldNames",
                            &params.printFieldNames,
                            "Print field names for benchmark output (including "
                            "any arguments to --fieldNames).");
      cmdLineProc.setOption("printTrilinosTestStuff",
                            "noPrintTrilinosTestStuff",
                            &params.printTrilinosTestStuff,
                            "Print output that makes the Trilinos test "
                            "framework happy (but makes benchmark results "
                            "parsing scripts unhappy)");
      cmdLineProc.parse(argc, argv);
    }
    catch(Teuchos::CommandLineProcessor::UnrecognizedOption& e) {
      err << "Unrecognized command-line option: " << e.what() << endl;
      throw e;
    }
    catch(Teuchos::CommandLineProcessor::HelpPrinted& e) {
      printedHelp = true;
      return params; // Don't verify parameters in this case
    }

    TEUCHOS_TEST_FOR_EXCEPTION
      (params.numRows <= 0, std::invalid_argument, "Number of "
       "rows must be positive, but you set --numRows=" <<
       params.numRows << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (params.numCols <= 0, std::invalid_argument, "Number of "
       "columns must be positive, but you set --numCols=" <<
       params.numCols << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (params.numRows < params.numCols, std::invalid_argument,
       "Number of rows must be >= number of columns, but "
       "--numRows=" << params.numRows << " and --numCols=" <<
       params.numCols << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (params.benchmark && params.numTrials < 1,
       std::invalid_argument, "If you set --benchmark, then the "
       "number of trials must be positive, but you set --numTrials="
       << params.numTrials << ".");
#ifndef HAVE_TPETRATSQR_COMPLEX
    TEUCHOS_TEST_FOR_EXCEPTION
      (params.testComplex, std::invalid_argument, "Complex "
       "arithmetic support was not enabled at configure time, "
       "but you set --testComplex.");
#endif // HAVE_TPETRATSQR_COMPLEX
    return params;
  }
} // namespace (anonymous)

int
main(int argc, char *argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;

  // Fetch command-line parameters.
  bool printedHelp = false;
  auto params = parseOptions(argc, argv, cerr, printedHelp);
  if(printedHelp) {
    return EXIT_SUCCESS;
  }
  bool success = false;
  constexpr bool actually_print_caught_exceptions = true;
  try {
    Kokkos::ScopeGuard kokkosScope(argc, argv);
    if(params.benchmark) {
      benchmark(cout, params);
    }
    // We allow the same run to do both benchmark and verify.
    if(params.verify) {
      verify(cout, params);
    }
    success = true;
    if(params.printTrilinosTestStuff) {
      // The Trilinos test framework expects a message like this.
      cout << "\nEnd Result: TEST PASSED" << endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS
    (actually_print_caught_exceptions, cerr, success);
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
