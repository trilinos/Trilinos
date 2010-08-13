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

#include "AnasaziConfigDefs.hpp"

#include "Teuchos_ConfigDefs.hpp" // HAVE_MPI
#ifdef HAVE_MPI
#  include "Teuchos_GlobalMPISession.hpp"
#  include "Teuchos_oblackholestream.hpp"
#endif // HAVE_MPI
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Time.hpp"

#include "Tsqr_CombineBenchmark.hpp"
#include "Tsqr_CombineTest.hpp"

#ifdef HAVE_TSQR_COMPLEX
#  include <complex>
#endif // HAVE_TSQR_COMPLEX

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR { 
  namespace Trilinos { 
    namespace Test {
      
      static const std::string& docStringFileName() {
	return TSQR_TEST_COMBINE_DOC_STRING_FILENAME;
      }

      // mfh 12 Aug 2010: docString should not be a constant string
      // literal if it's longer than some fixed, platform-dependent
      // amount (e.g., 2048 bytes).
      static const char docString[] = "This program tests TSQR::Combine.  "
	"Accuracy and performance tests are included.";

      using Teuchos::RCP;
      using Teuchos::Tuple;

      /// \class CombineTestParameters
      /// \brief Encapsulates values of command-line parameters
      ///
      struct CombineTestParameters {
	CombineTestParameters () :
	  verify (true),
	  benchmark (false),
	  numRows (1000),     // Number of rows in the test matrix
	  numCols (10),       // Number of columns in the test matrix
	  numTrials (10),     // Number of trials (action==Benchmark only)
	  debug (false),      // Whether to print debugging output to stderr
#ifdef HAVE_ANASAZI_COMPLEX
	  testComplex (true), // Whether to test complex-arithmetic routines
#endif // HAVE_ANASAZI_COMPLEX
	{}

	bool benchmark, verify;
	int numRows, numCols, numTrials;
	bool verbose, debug;
#ifdef HAVE_ANASAZI_COMPLEX
	bool testComplex;
#endif // HAVE_ANASAZI_COMPLEX
      };

      /// \brief Benchmark TSQR::Combine
      ///
      /// \param out [out] output stream for benchmark results.
      ///   It will only be used on rank 0.
      ///
      /// \param params [in] test parameter struct.  This method reads
      ///   the following field: numRows, numCols, numTrials,
      ///   testComplex.
      ///
      /// \warning Call only on (MPI) rank 0.  Otherwise, you'll run
      ///   the benchmark on every MPI rank simultaneously, but only
      ///   report results on rank 0.
      static void
      benchmarkCombine (std::ostream& out,
			const CombineTestParameters& params)
      {
	typedef Teuchos::Time timer_type;
	typedef int ordinal_type;

	const ordinal_type numRows = params.numRows;
	const ordinal_type numCols = params.numCols;
	const ordinal_type numTrials = params.numTrials;
#ifdef HAVE_ANASAZI_COMPLEX
	const bool testComplex = params.testComplex;
#else
	const bool testComplex = false;
#endif // HAVE_ANASAZI_COMPLEX

	std::vector<int> seed(4);
	const bool useSeedValues = false;
	TSQR::Test::benchmarkCombine< timer_type > (out, numRows, numCols, numTrials,
						    seed, useSeedValues, testComplex);
      }

      /// \brief Verify TSQR::Combine
      ///
      /// \param out [out] output stream for benchmark results.
      ///   It will only be used on rank 0.
      ///
      /// \param params [in] test parameter struct.  This method reads
      ///   the following field: numRows, numCols, numTrials,
      ///   testComplex.
      ///
      /// \warning Call only on (MPI) rank 0.  Otherwise, you'll run
      ///   the verification routine on every MPI rank simultaneously,
      ///   but only report results on rank 0.
      static void 
      verifyCombine (std::ostream& out,
		     const CombineTestParameters& params)
      {
	typedef Teuchos::Time timer_type;
	typedef int ordinal_type;

	const ordinal_type numRows = params.numRows;
	const ordinal_type numCols = params.numCols;
	const ordinal_type numTrials = params.numTrials;
#ifdef HAVE_ANASAZI_COMPLEX
	const bool testComplex = params.testComplex;
#else
	const bool testComplex = false;
#endif // HAVE_ANASAZI_COMPLEX
	const bool simulateSequentialTsqr = false;
	const bool debug = false;
	TSQR::Test::verifyCombine (numRows, numCols, testComplex, 
				   simulateSequentialTsqr, debug);
      }

      /// \brief Parse command-line options for this test
      ///
      /// \param argc [in] As usual in C(++)
      /// \param argv [in] As usual in C(++)
      /// \param allowedToPrint [in] Whether this (MPI) process is allowed
      ///   to print to stdout/stderr.  Different per (MPI) process.
      /// \param printedHelp [out] Whether this (MPI) process printed the
      ///   "help" display (summary of command-line options)
      ///
      /// \return Encapsulation of command-line options 
      static CombineTestParameters
      parseOptions (int argc, 
		    char* argv[], 
		    const bool allowedToPrint, 
		    bool& printedHelp)
      {
	using std::cerr;
	using std::endl;

	printedHelp = false;

	// Command-line parameters, set to their default values.
	CombineTestParameters params;
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
	  cmdLineProc.setOption ("debug", 
				 "nodebug", 
				 &params.debug, 
				 "Print debugging information");
	  cmdLineProc.setOption ("nrows", 
				 &params.numRows, 
				 "Number of rows in the test matrix");
	  cmdLineProc.setOption ("ncols", 
				 &params.numCols, 
				 "Number of columns in the test matrix");
	  cmdLineProc.setOption ("ntrials", 
				 &params.numTrials, 
				 "Number of trials (only used when \"--benchmark\"");
#ifdef HAVE_ANASAZI_COMPLEX
	  cmdLineProc.setOption ("complex", 
				 "nocomplex",
				 &params.testComplex,
				 "Test complex arithmetic, as well as real");
#endif // HAVE_ANASAZI_COMPLEX
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

	// Validate.  TODO (mfh 08 Jul 2010) Figure out how to do this with
	// ParameterList validators.
	if (params.numRows <= 0)
	  throw std::invalid_argument ("Number of rows must be positive");
	else if (params.numCols <= 0)
	  throw std::invalid_argument ("Number of columns must be positive");
	else if (params.numRows < params.numCols)
	  throw std::invalid_argument ("Number of rows must be >= number of columns");
	else if (params.benchmark && params.numTrials < 1)
	  throw std::invalid_argument ("Benchmark requires numTrials >= 1");

	return params;
      }

    } // namespace Test
  } // namespace Trilinos
} // namespace TSQR


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int 
main (int argc, char *argv[]) 
{
  using Teuchos::RCP;
  using TSQR::Trilinos::Test::CombineTestParameters;
  using TSQR::Trilinos::Test::parseOptions;

#ifdef HAVE_MPI
  typedef RCP< const Teuchos::Comm<int> > comm_ptr;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  comm_ptr comm = Teuchos::DefaultComm<int>::getComm();
  const int myRank = comm->getRank();
  // Only Rank 0 gets to write to stdout.  The other MPI process ranks
  // send their output to something that looks like /dev/null (and
  // likely is, on Unix-y operating systems).
  std::ostream& out = (myRank == 0) ? std::cout : blackhole;
  // Only Rank 0 performs the tests.
  const bool performingTests = (myRank==0);

#else // Don't HAVE_MPI: single-node test

  const bool performingTests = true;
  std::ostream& out = std::cout;
#endif // HAVE_MPI

  // Fetch command-line parameters.
  bool printedHelp = false;
  CombineTestParameters params = 
    parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp)
    return 0;

  if (performingTests)
    {
      using std::endl;

      if (params.benchmark)
	TSQR::Trilinos::Test::benchmarkCombine (out, params);
      // We allow the same run to do both benchmark and verify.
      if (params.verify)
	TSQR::Trilinos::Test::verifyCombine (out, params);

      // The Trilinos test framework expects a message like this.
      out << "\nEnd Result: TEST PASSED" << endl;
    }

  return 0;
}


