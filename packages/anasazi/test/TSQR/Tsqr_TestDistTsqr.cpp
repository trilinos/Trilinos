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
#ifdef HAVE_MPI
#  include "Teuchos_GlobalMPISession.hpp"
#  include "Teuchos_oblackholestream.hpp"
#endif // HAVE_MPI

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Time.hpp"

#include "Tsqr_Config.hpp"
#include "Tsqr_TeuchosMessenger.hpp"
#include "Tsqr_ParTest.hpp"

#ifdef HAVE_TSQR_COMPLEX
#  include <complex>
#endif // HAVE_TSQR_COMPLEX

#include <sstream>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR { 
  namespace Trilinos { 
    namespace Test {

      using Teuchos::RCP;
      using Teuchos::Tuple;

      const char docString[] = "This program tests TSQR::DistTsqr, which "
	"implements the internode-parallel part of TSQR (TSQR::Tsqr).  "
	"Accuracy and performance tests are included.";

      static void
      verifyDistTsqr (RCP< Teuchos::Comm<int> > comm,
		      const int numCols,
		      const bool testComplex,
		      const bool humanReadable,
		      const bool debug,
		      std::ostream& out,
		      std::ostream& err,
		      std::vector<int>& seed,
		      const bool useSeed)
      {
	if (! useSeed)
	  {
	    seed.resize (4);
	    seed[0] = 0;
	    seed[1] = 0;
	    seed[2] = 0;
	    seed[3] = 1;
	  }
	using TSQR::Trilinos::TeuchosMessenger;
	const bool testReal = true;

	if (testReal)
	  {
	    {
	      typedef float scalar_type;
	      std::string scalarTypeName ("float");
	      TeuchosMessenger< scalar_type > scalarComm (comm);
	      typedef ParTsqrVerifier< int, scalar_type > verifier_type;
	      verifier_type verifier (scalarComm, seed, scalarTypeName, 
				      humanReadable, debug, out, err);
	      verifier.verify (numCols);
	      verifier.getSeed (seed);
	    }
	    {
	      typedef double scalar_type;
	      std::string scalarTypeName ("double");
	      TeuchosMessenger< scalar_type > scalarComm (comm);
	      typedef ParTsqrVerifier< int, scalar_type > verifier_type;
	      verifier_type verifier (scalarComm, seed, scalarTypeName, 
				      humanReadable, debug, out, err);
	      verifier.verify (numCols);
	      verifier.getSeed (seed);
	    }
	  }
	
	if (testComplex)
	  {
#ifdef HAVE_TSQR_COMPLEX
	    using std::complex;
	    {
	      typedef complex<float> scalar_type;
	      std::string scalarTypeName ("complex<float>");
	      TeuchosMessenger< scalar_type > scalarComm (comm);
	      typedef ParTsqrVerifier< int, scalar_type > verifier_type;
	      verifier_type verifier (scalarComm, seed, scalarTypeName, 
				      humanReadable, debug, out, err);
	      verifier.verify (numCols);
	      verifier.getSeed (seed);
	    }
	    {
	      typedef complex<double> scalar_type;
	      std::string scalarTypeName ("complex<double>");
	      TeuchosMessenger< scalar_type > scalarComm (comm);
	      typedef ParTsqrVerifier< int, scalar_type > verifier_type;
	      verifier_type verifier (scalarComm, seed, scalarTypeName, 
				      humanReadable, debug, out, err);
	      verifier.verify (numCols);
	      verifier.getSeed (seed);
	    }
#else // Don't HAVE_TSQR_COMPLEX
	    throw std::logic_error("TSQR was not built with complex "
				   "arithmetic support");
#endif // HAVE_TSQR_COMPLEX
	  }
      }


      static void
      benchmarkDistTsqr (RCP< Teuchos::Comm<int> > comm,
			 const int numCols,
			 const int numTrials,
			 const bool testComplex,
			 const bool humanReadable,
			 const bool debug,
			 std::ostream& out,
			 std::ostream& err,
			 std::vector<int>& seed,
			 const bool useSeed)
      {
	if (! useSeed)
	  {
	    seed.resize (4);
	    seed[0] = 0;
	    seed[1] = 0;
	    seed[2] = 0;
	    seed[3] = 1;
	  }
	using TSQR::Test::ParTsqrBenchmarker;
	using TSQR::Trilinos::TeuchosMessenger;
	TeuchosMessenger< double > doubleComm (comm);
	const bool testReal = true;

	if (testReal)
	  {
	    {
	      typedef float scalar_type;
	      std::string scalarTypeName ("float");
	      TeuchosMessenger< scalar_type > scalarComm (comm);
	      typedef ParTsqrBenchmarker< int, scalar_type, timer_type > 
		benchmarker_type;
	      benchmarker_type benchmarker (scalarComm, doubleComm, seed, 
					    scalarTypeName, humanReadable,
					    debug, out, err);
	      benchmarker.benchmark (numTrials, numCols);
	      benchmarker.getSeed (seed);
	    }
	    {
	      typedef double scalar_type;
	      std::string scalarTypeName ("double");
	      TeuchosMessenger< scalar_type > scalarComm (comm);
	      typedef ParTsqrBenchmarker< int, scalar_type, timer_type > 
		benchmarker_type;
	      benchmarker_type benchmarker (scalarComm, doubleComm, seed, 
					    scalarTypeName, humanReadable,
					    debug, out, err);
	      benchmarker.benchmark (numTrials, numCols);
	      benchmarker.getSeed (seed);
	    }
	  }
	
	if (testComplex)
	  {
#ifdef HAVE_TSQR_COMPLEX
	    using std::complex;
	    {
	      typedef complex<float> scalar_type;
	      std::string scalarTypeName ("complex<float>");
	      TeuchosMessenger< scalar_type > scalarComm (comm);
	      typedef ParTsqrBenchmarker< int, scalar_type, timer_type > 
		benchmarker_type;
	      benchmarker_type benchmarker (scalarComm, doubleComm, seed, 
					    scalarTypeName, humanReadable,
					    debug, out, err);
	      benchmarker.benchmark (numTrials, numCols);
	      benchmarker.getSeed (seed);
	    }
	    {
	      typedef complex<double> scalar_type;
	      std::string scalarTypeName ("complex<double>");
	      TeuchosMessenger< scalar_type > scalarComm (comm);
	      typedef ParTsqrBenchmarker< int, scalar_type, timer_type > 
		benchmarker_type;
	      benchmarker_type benchmarker (scalarComm, doubleComm, seed, 
					    scalarTypeName, humanReadable,
					    debug, out, err);
	      benchmarker.benchmark (numTrials, numCols);
	      benchmarker.getSeed (seed);
	    }
#else // Don't HAVE_TSQR_COMPLEX
	    throw std::logic_error("TSQR was not built with complex "
				   "arithmetic support");
#endif // HAVE_TSQR_COMPLEX
	  }
      }


      /// \class DistTsqrTestParameters
      /// \brief Encapsulates values of command-line parameters
      ///
      struct DistTsqrTestParameters {
	DistTsqrTestParameters () :
	  numCols (10), 
	  numTrials (10), 
#ifdef HAVE_TSQR_COMPLEX
	  testComplex (true),
#endif // HAVE_TSQR_COMPLEX
	  verify (true),
	  benchmark (false),
	  debug (false),
	  humanReadable (false)
	{}

	int numCols, numTrials;
#ifdef HAVE_TSQR_COMPLEX
	bool testComplex;
#endif // HAVE_TSQR_COMPLEX
	bool verify, benchmark, debug, humanReadable;
      };

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
      static DistTsqrTestParameters
      parseOptions (int argc, 
		    char* argv[], 
		    const bool allowedToPrint, 
		    bool& printedHelp)
      {
	using std::cerr;
	using std::endl;

	printedHelp = false;

	// Command-line parameters, set to their default values.
	DistTsqrTestParameters params;
	try {
	  Teuchos::CommandLineProcessor cmdLineProc (/* throwExceptions=*/ true, 
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
	  cmdLineProc.setOption ("human-readable",
				 "machine-readable",
				 &params.humanReadable,
				 "If set, make output easy to read by humans "
				 "(but hard to parse)");
	  cmdLineProc.setOption ("ncols", 
				 &params.numCols, 
				 "Number of columns in the test matrix");
	  cmdLineProc.setOption ("ntrials", 
				 &params.numTrials, 
				 "Number of trials (only used when \"--benchmark\"");
#ifdef HAVE_TSQR_COMPLEX
	  cmdLineProc.setOption ("complex", 
				 "nocomplex",
				 &params.testComplex,
				 "Test complex arithmetic, as well as real");
#endif // HAVE_TSQR_COMPLEX
	  cmdLineProc.parse (argc, argv);
	} 
	catch (Teuchos::CommandLineProcessor::UnrecognizedOption& e) { 
	  if (allowedToPrint)
	    cerr << "Unrecognized command-line option: " << e.what() << endl;
	  throw e;
	}
	catch (Teuchos::CommandLineProcessor::HelpPrinted& e) { 
	  printedHelp = true;
	} 

	// Validate command-line options.  We provide default values
	// for unset options, so we don't have to validate those.
	if (params.numCols <= 0)
	  throw std::invalid_argument ("Number of columns must be positive");
	else if (params.benchmark && params.numTrials < 1)
	  throw std::invalid_argument ("\"--benchmark\" option requires numTrials >= 1");

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
  using TSQR::Trilinos::Test::DistTsqrTestParameters;
  using TSQR::Trilinos::Test::parseOptions;

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

#else // Don't HAVE_MPI: single-node test

  const bool allowedToPrint = true;
  std::ostream& out = std::cout;
  std::ostream& err = std::cerr;
#endif // HAVE_MPI

  // Fetch command-line parameters.
  bool printedHelp = false;
  DistTsqrTestParameters params = 
    parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp)
    return 0;

#ifdef HAVE_TSQR_COMPLEX
  const bool testComplex = params.testComplex;
#else // Don't HAVE_TSQR_COMPLEX
  const bool testComplex = false;
#endif // HAVE_TSQR_COMPLEX

  if (params.verify)
    {
      using TSQR::Trilinos::Test::verifyDistTsqr;
      std::vector<int> seed(4);
      const bool useSeed = false;
      verifyDistTsqr (comm, params.numCols, testComplex, params.humanReadable,
		      params.debug, out, err, seed, useSeed);
    }
  
  if (params.benchmark)
    {
      using TSQR::Trilinos::Test::benchmarkDistTsqr;
      std::vector<int> seed(4);
      const bool useSeed = false;
      benchmarkDistTsqr (comm, params.numCols, params.numTrials, testComplex,
			 params.humanReadable, params.debug, out, err, seed,
			 useSeed);
    }

  if (allowedToPrint)
    out << "\nEnd Result: TEST PASSED" << std::endl;
  return 0;
}


