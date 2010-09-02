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
#include "Teuchos_RCP.hpp"
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

using TSQR::MessengerBase;
using TSQR::Trilinos::TeuchosMessenger;
using TSQR::Test::DistTsqrVerifier;
using TSQR::Test::DistTsqrBenchmarker;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_implicit_cast;
using Teuchos::Tuple;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define TSQR_TEST_DIST_TSQR( theType, typeString )			\
do {							         	\
  typedef theType scalar_type;					        \
  typedef MessengerBase< scalar_type > base_messenger_type;	        \
  typedef RCP< base_messenger_type > base_messenger_ptr;		\
  typedef TeuchosMessenger< scalar_type > derived_messenger_type;       \
  typedef RCP< derived_messenger_type > derived_messenger_ptr;		\
  typedef DistTsqrVerifier< int, scalar_type > verifier_type;		\
									\
  std::string scalarTypeName (typeString);				\
  derived_messenger_ptr scalarCommDerived (new derived_messenger_type (comm)); \
  base_messenger_ptr scalarComm =					\
    rcp_implicit_cast< base_messenger_type > (scalarCommDerived);	\
  verifier_type verifier (scalarComm, seed, scalarTypeName,		\
			  out, err, testFactorExplicit,			\
			  testFactorImplicit, humanReadable, debug);	\
  verifier.verify (numCols);						\
  verifier.getSeed (seed);						\
} while(false)


#define TSQR_BENCHMARK_DIST_TSQR( theType, typeString )			\
do {									\
  typedef theType scalar_type;						\
  typedef MessengerBase< scalar_type > base_messenger_type;	        \
  typedef RCP< base_messenger_type > base_messenger_ptr;		\
  typedef TeuchosMessenger< scalar_type > derived_messenger_type;       \
  typedef RCP< derived_messenger_type > derived_messenger_ptr;		\
  typedef DistTsqrBenchmarker< int, scalar_type, timer_type >		\
    benchmarker_type;							\
									\
  std::string scalarTypeName (typeString);				\
  derived_messenger_ptr scalarCommDerived (new derived_messenger_type (comm)); \
  base_messenger_ptr scalarComm =					\
    rcp_implicit_cast< base_messenger_type > (scalarCommDerived);	\
  benchmarker_type benchmarker (scalarComm, doubleComm, seed,		\
				scalarTypeName, out, err,		\
				testFactorExplicit, testFactorImplicit, \
				humanReadable, debug);			\
  benchmarker.benchmark (numTrials, numCols);				\
  benchmarker.getSeed (seed);						\
} while(false)


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

const char docString[] = "This program tests TSQR::DistTsqr, which "
  "implements the internode-parallel part of TSQR (TSQR::Tsqr).  "
  "Accuracy and performance tests are included.";

/// \class DistTsqrTestParameters
/// \brief Encapsulates values of command-line parameters
///
struct DistTsqrTestParameters {
  int numCols, numTrials;
#ifdef HAVE_TSQR_COMPLEX
  bool testComplex;
#endif // HAVE_TSQR_COMPLEX
  bool verify, benchmark;
  bool testFactorExplicit, testFactorImplicit;
  bool humanReadable, debug;

  DistTsqrTestParameters () :
    numCols (10), 
    numTrials (10), 
#ifdef HAVE_TSQR_COMPLEX
    testComplex (true),
#endif // HAVE_TSQR_COMPLEX
    verify (true),
    benchmark (false),
    testFactorExplicit (true),
    testFactorImplicit (true),
    humanReadable (false),
    debug (false)
  {}
};

static void
verify (RCP< const Teuchos::Comm<int> > comm,
	const DistTsqrTestParameters& params,
	std::ostream& out,
	std::ostream& err,
	std::vector<int>& seed,
	const bool useSeed)
{
#ifdef HAVE_TSQR_COMPLEX
  const bool testComplex = params.testComplex;
#else // Don't HAVE_TSQR_COMPLEX
  const bool testComplex = false;
#endif // HAVE_TSQR_COMPLEX

  const int numCols = params.numCols;
  const bool testFactorExplicit = params.testFactorExplicit;
  const bool testFactorImplicit = params.testFactorImplicit;
  const bool humanReadable = params.humanReadable;
  const bool debug = params.debug;

  if (! useSeed)
    {
      seed.resize (4);
      seed[0] = 0;
      seed[1] = 0;
      seed[2] = 0;
      seed[3] = 1;
    }
  const bool testReal = true;

  if (testReal)
    {
      TSQR_TEST_DIST_TSQR( float, "float" );
      TSQR_TEST_DIST_TSQR( double, "double" );
    }
	
  if (testComplex)
    {
#ifdef HAVE_TSQR_COMPLEX
      using std::complex;

      TSQR_TEST_DIST_TSQR( complex<float>, "complex<float>" );
      TSQR_TEST_DIST_TSQR( complex<double>, "complex<double>" );

#else // Don't HAVE_TSQR_COMPLEX
      throw std::logic_error("TSQR was not built with complex "
			     "arithmetic support");
#endif // HAVE_TSQR_COMPLEX
    }
}


static void
benchmark (RCP< const Teuchos::Comm<int> > comm,
	   const DistTsqrTestParameters& params,
	   std::ostream& out,
	   std::ostream& err,
	   std::vector<int>& seed,
	   const bool useSeed)
{
  typedef Teuchos::Time timer_type;

#ifdef HAVE_TSQR_COMPLEX
  const bool testComplex = params.testComplex;
#else // Don't HAVE_TSQR_COMPLEX
  const bool testComplex = false;
#endif // HAVE_TSQR_COMPLEX

  const int numCols = params.numCols;
  const int numTrials = params.numTrials;
  const bool testFactorExplicit = params.testFactorExplicit;
  const bool testFactorImplicit = params.testFactorImplicit;
  const bool humanReadable = params.humanReadable;
  const bool debug = params.debug;

  if (! useSeed)
    {
      seed.resize (4);
      seed[0] = 0;
      seed[1] = 0;
      seed[2] = 0;
      seed[3] = 1;
    }
  RCP< MessengerBase< double > > doubleComm = 
    rcp_implicit_cast< MessengerBase< double > > (RCP< TeuchosMessenger< double > > (new TeuchosMessenger< double > (comm)));
  const bool testReal = true;

  if (testReal)
    {
      TSQR_BENCHMARK_DIST_TSQR( float, "float" );
      TSQR_BENCHMARK_DIST_TSQR( double, "double" );
    }
	
  if (testComplex)
    {
#ifdef HAVE_TSQR_COMPLEX
      using std::complex;

      TSQR_BENCHMARK_DIST_TSQR( complex<float>, "complex<float>" );
      TSQR_BENCHMARK_DIST_TSQR( complex<double>, "complex<double>" );

#else // Don't HAVE_TSQR_COMPLEX
      throw std::logic_error("TSQR was not built with complex "
			     "arithmetic support");
#endif // HAVE_TSQR_COMPLEX
    }
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
    cmdLineProc.setOption ("implicit",
			   "noimplicit",
			   &params.testFactorImplicit,
			   "Test DistTsqr\'s factor() and explicit_Q()");
    cmdLineProc.setOption ("explicit",
			   "noexplicit",
			   &params.testFactorExplicit,
			   "Test DistTsqr\'s factorExplicit()");
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int 
main (int argc, char *argv[]) 
{
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

  if (params.verify)
    {
      std::vector<int> seed(4);
      const bool useSeed = false;

      verify (comm, params, out, err, seed, useSeed);
    }
  
  if (params.benchmark)
    {
      std::vector<int> seed(4);
      const bool useSeed = false;

      benchmark (comm, params, out, err, seed, useSeed);
    }

  if (allowedToPrint)
    out << "\nEnd Result: TEST PASSED" << std::endl;
  return 0;
}


