// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include <Tsqr_FullTsqrTest.hpp>

#ifdef HAVE_MPI
#  include <Teuchos_GlobalMPISession.hpp>
#  include <Teuchos_oblackholestream.hpp>
#endif // HAVE_MPI
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_TSQR_COMPLEX
#  include <complex>
#endif // HAVE_TSQR_COMPLEX

typedef int ordinal_type;

namespace {
  //
  // Documentation string to print out if --help is a command-line argument.
  //
  const char docString[] = "This program tests correctness and accuracy of "
    "TSQR::Tsqr, which is the full implementation of TSQR.";

  //
  // Encapsulation of all command-line parameters.
  //
  struct CmdLineOptions {
    //
    // Given a default valid parameter list from FullTsqrVerifierCaller,
    // fill in the command-line options with their default values.
    //
    CmdLineOptions (const Teuchos::RCP<const Teuchos::ParameterList>& testParams) :
      cacheSizeHint (testParams->get<size_t> ("cacheSizeHint")),
      numRowsLocal (testParams->get<int> ("numRowsLocal")),
      numCols (testParams->get<int> ("numCols")),
      contiguousCacheBlocks (testParams->get<bool> ("contiguousCacheBlocks")),
      testFactorExplicit (testParams->get<bool> ("testFactorExplicit")),
      testRankRevealing (testParams->get<bool> ("testRankRevealing")),
      printFieldNames (testParams->get<bool> ("printFieldNames")),
      printResults (testParams->get<bool> ("printResults")),
      failIfInaccurate (testParams->get<bool> ("failIfInaccurate")),
      debug (testParams->get<bool> ("debug")),
#ifdef HAVE_TSQR_COMPLEX
      testComplex (false),
#endif // HAVE_TSQR_COMPLEX
      testReal (false) // default is not to test _anything_
    {}

    size_t cacheSizeHint;
    ordinal_type numRowsLocal;
    ordinal_type numCols;
    bool contiguousCacheBlocks;
    bool testFactorExplicit;
    bool testRankRevealing;
    bool printFieldNames;
    bool printResults;
    bool failIfInaccurate;
    bool debug;
#ifdef HAVE_TSQR_COMPLEX
    bool testComplex;
#endif // HAVE_TSQR_COMPLEX
    bool testReal;

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
    // \param allowedToPrint [in] Whether this (MPI) process is allowed
    //   to print to stdout/stderr.  Different per (MPI) process.
    //
    // \param printedHelp [out] Whether this (MPI) process printed the
    //   "help" display (summary of command-line options)
    //
    // \param testParams [in] List of test parameters for the
    //   FullTsqrVerifierCaller.
    //
    // \return Whether help was printed.
    bool
    read (int argc, 
	  char* argv[], 
	  const Teuchos::RCP<const Teuchos::ParameterList>& defaultParams,
	  const bool allowedToPrint)
    {
      using std::cerr;
      using std::endl;

      try {
	Teuchos::CommandLineProcessor cmdLineProc (/* throwExceptions=*/ true, 
						   /* recognizeAllOptions=*/ true);
	cmdLineProc.setDocString (docString);
	cmdLineProc.setOption ("testReal",
			       "noTestReal",
			       &testReal,
			       "Test real Scalar types");
#ifdef HAVE_TSQR_COMPLEX
	cmdLineProc.setOption ("testComplex",
			       "noTestComplex",
			       &testComplex,
			       "Test complex Scalar types");
#endif // HAVE_TSQR_COMPLEX
	// CommandLineProcessor takes int arguments, but not size_t
	// arguments, so we have to read in the argument as an int and
	// convert back to size_t later.
	int cacheSizeHintAsInt = cacheSizeHint;
	cmdLineProc.setOption ("cacheSizeHint",
			       &cacheSizeHintAsInt,
			       defaultParams->getEntry("cacheSizeHint").docString().c_str());
	cmdLineProc.setOption ("numRowsLocal",
			       &numRowsLocal,
			       defaultParams->getEntry("numRowsLocal").docString().c_str());
	cmdLineProc.setOption ("numCols",
			       &numCols,
			       defaultParams->getEntry("numCols").docString().c_str());
	cmdLineProc.setOption ("contiguousCacheBlocks",
			       "noContiguousCacheBlocks",
			       &contiguousCacheBlocks,
			       defaultParams->getEntry("contiguousCacheBlocks").docString().c_str());
	cmdLineProc.setOption ("testFactorExplicit",
			       "noTestFactorExplicit",
			       &testFactorExplicit,
			       defaultParams->getEntry("testFactorExplicit").docString().c_str());
	cmdLineProc.setOption ("testRankRevealing",
			       "noTestRankRevealing",
			       &testRankRevealing,
			       defaultParams->getEntry("testRankRevealing").docString().c_str());
	cmdLineProc.setOption ("printFieldNames",
			       "noPrintFieldNames",
			       &printFieldNames,
			       defaultParams->getEntry("printFieldNames").docString().c_str());
	cmdLineProc.setOption ("printResults",
			       "noPrintResults",
			       &printResults,
			       defaultParams->getEntry("printResults").docString().c_str());
	cmdLineProc.setOption ("failIfInaccurate",
			       "noFailIfInaccurate",
			       &failIfInaccurate,
			       defaultParams->getEntry("failIfInaccurate").docString().c_str());
	cmdLineProc.setOption ("debug", 
			       "nodebug", 
			       &debug, 
			       defaultParams->getEntry("debug").docString().c_str());
	cmdLineProc.parse (argc, argv);
	cacheSizeHint = static_cast<size_t> (cacheSizeHintAsInt);
      } 
      catch (Teuchos::CommandLineProcessor::UnrecognizedOption& e) { 
	if (allowedToPrint)
	  cerr << "Unrecognized command-line option: " << e.what() << endl;
	throw e;
      }
      catch (Teuchos::CommandLineProcessor::HelpPrinted& e) { 
	return true;
      } 

      // Validate command-line options.  We provide default values
      // for unset options, so we don't have to validate those.
      TEST_FOR_EXCEPTION(numRowsLocal <= 0, std::invalid_argument,
			 "Number of rows per process must be positive.");
      TEST_FOR_EXCEPTION(numCols <= 0, std::invalid_argument,
			 "Number of columns must be positive.");
      return false; // Did not print help
    }
  };

  //
  // Given a default valid parameter list from FullTsqrVerifierCaller,
  // and the values of command-line options (that were read in from
  // the command line), return a parameter list describing the test.
  //
  Teuchos::RCP<Teuchos::ParameterList>
  testParameters (const Teuchos::RCP<const Teuchos::ParameterList>& validParams,
		  const CmdLineOptions& options)
  {
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;

    RCP<ParameterList> testParams = parameterList ("FullTsqrVerifier");
    testParams->set ("cacheSizeHint", options.cacheSizeHint);
    testParams->set ("numRowsLocal", options.numRowsLocal);
    testParams->set ("numCols", options.numCols);
    testParams->set ("testFactorExplicit", options.testFactorExplicit);
    testParams->set ("testRankRevealing", options.testRankRevealing);
    testParams->set ("contiguousCacheBlocks", options.contiguousCacheBlocks);
    testParams->set ("printFieldNames", options.printFieldNames);
    testParams->set ("printResults", options.printResults);
    testParams->set ("failIfInaccurate", options.failIfInaccurate);
    testParams->set ("debug", options.debug);

    testParams->validateParametersAndSetDefaults (*validParams);
    return testParams;
  }

  //
  // Return true if all tests were successful, else false.
  //
  bool
  test (int argc,
	char* argv[],
	const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
	const Teuchos::RCP<const Kokkos::SerialNode>& node,
	const bool allowedToPrint)
  {
    using TSQR::Test::NullCons;
    using TSQR::Test::Cons;
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    //
    // Get a default random seed, and set up the Caller (that iterates
    // the test over all Scalar types of interest).
    //
    typedef TSQR::Test::FullTsqrVerifierCaller caller_type;
    std::vector<int> randomSeed = caller_type::defaultRandomSeed ();
    caller_type caller (comm, node, randomSeed);

    //
    // Read command-line options
    //
    RCP<const ParameterList> defaultParams = caller.getValidParameterList();
    CmdLineOptions cmdLineOpts (defaultParams);
    const bool printedHelp = cmdLineOpts.read (argc, argv, defaultParams, allowedToPrint);
    // Don't run the tests (and do succeed) if help was printed.
    if (printedHelp)
      return true; 

    //
    // Use read-in command-line options to set up test parameters.
    //
    RCP<ParameterList> testParams = testParameters (defaultParams, cmdLineOpts);
    defaultParams = null; // save a little space 

    //
    // Define lists of Scalar types to test.  We keep separate lists
    // for real and complex types, since callers can control whether
    // each of these is tested independently on the command line.
    //
    typedef Cons<float, Cons<double, NullCons> > real_type_list;
#ifdef HAVE_TSQR_COMPLEX
    typedef Cons<std::complex<float>, Cons<std::complex<double>, NullCons> > complex_type_list;
#endif // HAVE_TSQR_COMPLEX

    //
    // Run the tests.  If the tests are set up to fail on
    // insufficiently inaccurate results, run() will throw an
    // exception in that case.  Otherwise, the tests return nothing,
    // and "succeed" if they don't crash or throw an exception.
    //
    // The testReal and testComplex options are read in at the command
    // line, but since they do not apply to all Scalar types, they
    // don't belong in testParams.
    //
    if (cmdLineOpts.testReal)
      caller.run<real_type_list> (testParams);
#ifdef HAVE_TSQR_COMPLEX
    if (cmdLineOpts.testComplex)
      caller.run<complex_type_list> (testParams);
#endif // HAVE_TSQR_COMPLEX

    return true; // for success
  }

} // namespace (anonymous)


int 
main (int argc, char* argv[])
{
  using TSQR::Test::NullCons;
  using TSQR::Test::Cons;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;

#ifdef HAVE_MPI
  typedef RCP< const Teuchos::Comm<int> > comm_ptr;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  comm_ptr comm = Teuchos::DefaultComm<int>::getComm();
  const int myRank = comm->getRank();
  // Only Rank 0 gets to write to cout and cerr.  The other MPI
  // process ranks send their output to a "black hole" (something that
  // acts like /dev/null, and may be /dev/null).
  const bool allowedToPrint = (myRank == 0);
  std::ostream& out = allowedToPrint ? std::cout : blackhole;
  std::ostream& err = allowedToPrint ? std::cerr : blackhole;
  // Make sure that err gets "used"
  (void) err;

#else // Don't HAVE_MPI: single-node test

  const bool allowedToPrint = true;
  std::ostream& out = std::cout;
  std::ostream& err = std::cerr;
#endif // HAVE_MPI

  RCP<ParameterList> nodeParams = 
    TSQR::Test::getValidNodeParameters<Kokkos::SerialNode> ();
  RCP<const Kokkos::SerialNode> node =
    TSQR::Test::getNode<Kokkos::SerialNode> (nodeParams);

  const bool success = test (argc, argv, comm, node, allowedToPrint);
  if (allowedToPrint && success)
    // The Trilinos test framework expects a message like this.
    out << "\nEnd Result: TEST PASSED" << endl;
  return 0;
}


