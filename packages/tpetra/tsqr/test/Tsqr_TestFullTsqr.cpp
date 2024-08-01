// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_FullTsqrTest.hpp"
#include "Tsqr_Test_MpiAndKokkosScope.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_TPETRATSQR_COMPLEX
#  include <complex>
#endif // HAVE_TPETRATSQR_COMPLEX

namespace {
  using Teuchos::parameterList;

  // Documentation string to print out if --help is a command-line
  // argument.
  const char docString[] = "This program tests correctness, "
    "accuracy, and/or performance of TSQR::Tsqr.  TSQR::Tsqr "
    "is the full implementation of TSQR, including both NodeTsqr "
    "(within a single MPI process) and DistTsqr (across MPI "
    "processes).";

  // Encapsulation of all command-line parameters.
  struct CmdLineOptions {
    // Given a default valid parameter list from
    // FullTsqrVerifierCaller, fill in the command-line options with
    // their default values.
    CmdLineOptions(const Teuchos::RCP<const Teuchos::ParameterList>& testParams) :
      cacheSizeHint(testParams->get<size_t>("Cache Size Hint")),
      numRowsLocal(testParams->get<int>("numRowsLocal")),
      numCols(testParams->get<int>("numCols")),
      numTrials(testParams->get<int>("numTrials")),
      contiguousCacheBlocks(testParams->get<bool>("contiguousCacheBlocks")),
      testFactorExplicit(testParams->get<bool>("testFactorExplicit")),
      testRankRevealing(testParams->get<bool>("testRankRevealing")),
      printFieldNames(testParams->get<bool>("printFieldNames")),
      printResults(testParams->get<bool>("printResults")),
      failIfInaccurate(testParams->get<bool>("failIfInaccurate")),
      nodeTsqr(testParams->get<std::string>("NodeTsqr")),
      verbose(testParams->get<bool>("verbose"))
    {}

    size_t cacheSizeHint = 0;
    int numRowsLocal = 10000;
    int numCols = 5;
    int numTrials = 100;
    bool contiguousCacheBlocks = false;
    bool testFactorExplicit = true;
    bool testRankRevealing = true;
    bool printFieldNames = true;
    bool printResults = true;
    bool failIfInaccurate = true;
    std::string nodeTsqr {"Default"};
#ifdef HAVE_TPETRATSQR_COMPLEX
    bool testComplex = true;
#else
    bool testComplex = false;
#endif // HAVE_TPETRATSQR_COMPLEX
    bool testReal = true;
    bool verbose = false;
    bool verify = true;
    bool benchmark = false;

    // \brief Read command-line options.
    //
    // We use Doxygen notation to document this function, but don't tell
    // Doxygen to generate documentation, since this method is local to
    // this test.
    //
    // \param argc [in] As usual in C(++).
    //
    // \param argv [in] As usual in C(++).
    //
    // \param testParams [in] List of test parameters for the
    //   FullTsqrVerifierCaller.
    //
    // \param err [out] Output stream to which to print error
    //   messages.  Different per (MPI) process.
    //
    // \return Whether help was printed.
    bool
    read(int argc,
         char* argv[],
         const Teuchos::RCP<const Teuchos::ParameterList>& defaultParams,
         std::ostream& err)
    {
      using Teuchos::CommandLineProcessor;
      using std::endl;

      try {
        const bool throwExceptions = true;
        const bool recognizeAllOptions = true;
        CommandLineProcessor cmdLineProc(throwExceptions,
                                         recognizeAllOptions);
        cmdLineProc.setDocString(docString);
        cmdLineProc.setOption("testReal",
                              "noTestReal",
                              &testReal,
                              "Test real Scalar types");
        cmdLineProc.setOption("testComplex",
                              "noTestComplex",
                              &testComplex,
                              "Test complex Scalar types; must be "
                              "false if complex Scalar types were "
                              "disabled at configure (pre-build) "
                              "time");
        // CommandLineProcessor takes int arguments, but not size_t
        // arguments, so we have to read in the argument as an int and
        // convert back to size_t later.
        int cacheSizeHintAsInt = cacheSizeHint;
        cmdLineProc.setOption("cacheSizeHint",
                              &cacheSizeHintAsInt,
                              defaultParams->getEntry
                              ("Cache Size Hint").docString().c_str());
        cmdLineProc.setOption("numRowsLocal",
                              &numRowsLocal,
                              defaultParams->getEntry
                              ("numRowsLocal").docString().c_str());
        cmdLineProc.setOption("numCols",
                              &numCols,
                              defaultParams->getEntry
                              ("numCols").docString().c_str());
        cmdLineProc.setOption("numTrials",
                              &numTrials,
                              defaultParams->getEntry
                              ("numTrials").docString().c_str());
        cmdLineProc.setOption("contiguousCacheBlocks",
                              "noContiguousCacheBlocks",
                              &contiguousCacheBlocks,
                              defaultParams->getEntry
                              ("contiguousCacheBlocks").docString().c_str());
        cmdLineProc.setOption("testFactorExplicit",
                              "noTestFactorExplicit",
                              &testFactorExplicit,
                              defaultParams->getEntry
                              ("testFactorExplicit").docString().c_str());
        cmdLineProc.setOption("testRankRevealing",
                              "noTestRankRevealing",
                              &testRankRevealing,
                              defaultParams->getEntry
                              ("testRankRevealing").docString().c_str());
        cmdLineProc.setOption("printFieldNames",
                              "noPrintFieldNames",
                              &printFieldNames,
                              defaultParams->getEntry
                              ("printFieldNames").docString().c_str());
        cmdLineProc.setOption("printResults",
                              "noPrintResults",
                              &printResults,
                              defaultParams->getEntry
                              ("printResults").docString().c_str());
        cmdLineProc.setOption("failIfInaccurate",
                              "noFailIfInaccurate",
                              &failIfInaccurate,
                              defaultParams->getEntry
                              ("failIfInaccurate").docString().c_str());
        cmdLineProc.setOption("NodeTsqr",
                              &nodeTsqr,
                              defaultParams->getEntry
                              ("NodeTsqr").docString().c_str());
        cmdLineProc.setOption("verbose",
                              "quiet",
                              &verbose,
                              defaultParams->getEntry
                              ("verbose").docString().c_str());
        cmdLineProc.setOption("verify",
                              "noverify",
                              &verify,
                              "Test accuracy");
        cmdLineProc.setOption("benchmark",
                              "nobenchmark",
                              &benchmark,
                              "Test performance");

        cmdLineProc.parse(argc, argv);
        cacheSizeHint = size_t(cacheSizeHintAsInt);
      }
      catch(Teuchos::CommandLineProcessor::UnrecognizedOption& e) {
        err << "Unrecognized command-line option: " << e.what()
            << endl;
        throw e;
      }
      catch(Teuchos::CommandLineProcessor::HelpPrinted& e) {
        return true;
      }

      // Validate command-line options.  We provide default values
      // for unset options, so we don't have to validate those.
      TEUCHOS_TEST_FOR_EXCEPTION
        (numRowsLocal <= 0, std::invalid_argument,
         "Number of rows per process must be positive.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (numCols <= 0, std::invalid_argument,
         "Number of columns must be positive, but you specified "
         "--numCols=" << numCols << ".");
      TEUCHOS_TEST_FOR_EXCEPTION
        (numTrials <= 0, std::invalid_argument,
         "Number of trials must be positive, but you specified "
         "--numTrials=" << numTrials << ".");
      return false; // Did not print help
    }
  };

  //
  // Given a default valid parameter list from FullTsqrVerifierCaller,
  // and the values of command-line options (that were read in from
  // the command line), return a parameter list describing the test.
  //
  Teuchos::RCP<Teuchos::ParameterList>
  testParameters(const Teuchos::RCP<const Teuchos::ParameterList>& validParams,
                 const CmdLineOptions& options)
  {
    auto testParams = parameterList("FullTsqrVerifier");
    testParams->set("Cache Size Hint", options.cacheSizeHint);
    testParams->set("numRowsLocal", options.numRowsLocal);
    testParams->set("numCols", options.numCols);
    testParams->set("numTrials", options.numTrials);
    testParams->set("testFactorExplicit",
                    options.testFactorExplicit);
    testParams->set("testRankRevealing", options.testRankRevealing);
    testParams->set("contiguousCacheBlocks",
                    options.contiguousCacheBlocks);
    testParams->set("printFieldNames", options.printFieldNames);
    testParams->set("printResults", options.printResults);
    testParams->set("failIfInaccurate", options.failIfInaccurate);
    testParams->set("NodeTsqr", options.nodeTsqr);
    testParams->set("verbose", options.verbose);

    testParams->validateParametersAndSetDefaults(*validParams);
    return testParams;
  }

  // Return true if all tests were successful, else false.
  bool
  test(int argc,
       char* argv[],
       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
       std::ostream& err)
  {
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // The Caller iterates the test over all Scalar types.
    using caller_type = TSQR::Test::FullTsqrVerifierCaller;
    caller_type caller(comm, caller_type::defaultRandomSeed());

    // Read command-line options
    auto defaultParams = caller.getValidParameterList();
    CmdLineOptions cmdLineOpts(defaultParams);
    const bool printedHelp =
      cmdLineOpts.read(argc, argv, defaultParams, err);
    // Don't run the tests (and do succeed) if help was printed.
    if(printedHelp) {
      return true;
    }

    // Use read-in command-line options to set up test parameters.
    auto testParams = testParameters(defaultParams, cmdLineOpts);
    defaultParams = null; // save a little space

    bool success = true;
    if (cmdLineOpts.verify) {
      // Verify accuracy of "full" TSQR.  If the tests are set up to
      // fail on insufficiently inaccurate results, run() will throw
      // an exception in that case.  Otherwise, the tests return
      // nothing, and "succeed" if they don't crash or throw an
      // exception.
      if (cmdLineOpts.testReal) {
        const bool ok = caller.verify<float, double>(testParams);
        success = success && ok;
      }
#ifdef HAVE_TPETRATSQR_COMPLEX
      if (cmdLineOpts.testComplex) {
        const bool ok =
          caller.verify<std::complex<float>,
                        std::complex<double>>(testParams);
        success = success && ok;
      }
#endif // HAVE_TPETRATSQR_COMPLEX
    }

    if (cmdLineOpts.benchmark) {
      if (cmdLineOpts.testReal) {
        caller.benchmark<float, double>(testParams);
      }
#ifdef HAVE_TPETRATSQR_COMPLEX
      if (cmdLineOpts.testComplex) {
        caller.benchmark<std::complex<float>,
                         std::complex<double>>(testParams);
      }
#endif // HAVE_TPETRATSQR_COMPLEX
    }

    return success;
  }
} // namespace (anonymous)


int
main(int argc, char* argv[])
{
  using std::endl;
  TSQR::Test::MpiAndKokkosScope testScope(&argc, &argv);
  auto comm = testScope.getComm();
  std::ostream& out = testScope.outStream();
  std::ostream& err = testScope.errStream();

  constexpr bool actually_print_caught_exceptions = true;
  bool success = false; // hopefully this will be true later
  try {
    success = test(argc, argv, comm, err);
    if (success) {
      // The Trilinos test framework expects a message like this.
      out << "\nEnd Result: TEST PASSED" << endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS
    (actually_print_caught_exceptions, err, success);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
