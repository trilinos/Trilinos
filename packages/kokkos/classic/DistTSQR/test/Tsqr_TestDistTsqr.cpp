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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Tsqr_ConfigDefs.hpp>

#ifdef HAVE_MPI
#  include <Teuchos_GlobalMPISession.hpp>
#  include <Teuchos_oblackholestream.hpp>
#endif // HAVE_MPI

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>

#include <Tsqr_ParTest.hpp>
#include <Tsqr_TeuchosMessenger.hpp>

#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
#  include <complex>
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX

#include <sstream>
#include <stdexcept>
#include <vector>

using TSQR::MessengerBase;
using TSQR::TeuchosMessenger;
using TSQR::Test::DistTsqrVerifier;
using TSQR::Test::DistTsqrBenchmarker;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_implicit_cast;
using Teuchos::Tuple;


template< class Ordinal, class Scalar >
class MessengerPairMaker {
public:
  typedef int ordinal_type;
  typedef Scalar scalar_type;

  typedef std::pair<RCP<MessengerBase<ordinal_type> >, RCP<MessengerBase<scalar_type> > > pair_type;

  static pair_type
  makePair (const RCP< const Teuchos::Comm<int> >& comm)
  {
    RCP<TeuchosMessenger<ordinal_type> > derivedOrdinalComm = 
      rcp (new TeuchosMessenger<ordinal_type> (comm));
    RCP<MessengerBase<ordinal_type> > ordinalComm = 
      rcp_implicit_cast<MessengerBase<ordinal_type> > (derivedOrdinalComm);
    RCP<TeuchosMessenger<scalar_type> > derivedScalarComm =
      rcp (new TeuchosMessenger<scalar_type> (comm));
    RCP<MessengerBase<scalar_type> > scalarComm = 
      rcp_implicit_cast<MessengerBase<scalar_type> > (derivedScalarComm);

    return std::make_pair (ordinalComm, scalarComm);
  }
};


#define TSQR_TEST_DIST_TSQR( ScalarType, typeString )			\
do {							         	\
  typedef int ordinal_type;						\
  typedef ScalarType scalar_type;					\
  typedef MessengerPairMaker<ordinal_type, scalar_type>::pair_type pair_type; \
  typedef DistTsqrVerifier<int, scalar_type> verifier_type;		\
									\
  std::string scalarTypeName (typeString);				\
  pair_type messPair = MessengerPairMaker< ordinal_type, scalar_type >::makePair (comm); \
  verifier_type verifier (messPair.first, messPair.second, seed,	\
			  scalarTypeName, out, err, 			\
			  testFactorExplicit, testFactorImplicit,	\
			  humanReadable, printMatrices, debug);		\
  verifier.verify (numCols, params.additionalFieldNames,		\
		   params.additionalData, params.printFieldNames);	\
  verifier.getSeed (seed);						\
} while(false)


#define TSQR_BENCHMARK_DIST_TSQR( theType, typeString )			\
do {									\
  typedef theType scalar_type;						\
  typedef MessengerBase< scalar_type > base_messenger_type;	        \
  typedef RCP< base_messenger_type > base_messenger_ptr;		\
  typedef TeuchosMessenger< scalar_type > derived_messenger_type;       \
  typedef RCP< derived_messenger_type > derived_messenger_ptr;		\
  typedef DistTsqrBenchmarker<int, scalar_type, timer_type>		\
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
  benchmarker.benchmark (numTrials, numCols,				\
			 params.additionalFieldNames,			\
			 params.additionalData,				\
			 params.printFieldNames);			\
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
  DistTsqrTestParameters () :
    numCols (10), 
    numTrials (10), 
    verify (false),
    benchmark (false),
    testReal (true),
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
    testComplex (true),
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
    testFactorExplicit (true),
    testFactorImplicit (true),
    printFieldNames (true),
    printTrilinosTestStuff (true),
    humanReadable (false),
    printMatrices (false),
    debug (false)
  {}

  std::string additionalFieldNames, additionalData;
  int numCols, numTrials;
  bool verify, benchmark;
  bool testReal;
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
  bool testComplex;
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
  bool testFactorExplicit, testFactorImplicit;
  bool printFieldNames, printTrilinosTestStuff;
  bool humanReadable, printMatrices, debug;
};

static void
verify (RCP< const Teuchos::Comm<int> > comm,
	const DistTsqrTestParameters& params,
	std::ostream& out,
	std::ostream& err,
	std::vector<int>& seed,
	const bool useSeed)
{
  const bool testReal = params.testReal;
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
  const bool testComplex = params.testComplex;
#else // Don't HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
  const bool testComplex = false;
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX

  const int numCols = params.numCols;
  const bool testFactorExplicit = params.testFactorExplicit;
  const bool testFactorImplicit = params.testFactorImplicit;
  const bool humanReadable = params.humanReadable;
  const bool printMatrices = params.printMatrices;
  const bool debug = params.debug;

  if (! useSeed)
    {
      seed.resize (4);
      seed[0] = 0;
      seed[1] = 0;
      seed[2] = 0;
      seed[3] = 1;
    }
  if (testReal)
    {
      TSQR_TEST_DIST_TSQR( float, "float" );
      TSQR_TEST_DIST_TSQR( double, "double" );
    }
  if (testComplex)
    {
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
      using std::complex;

      TSQR_TEST_DIST_TSQR( complex<float>, "complex<float>" );
      TSQR_TEST_DIST_TSQR( complex<double>, "complex<double>" );

#else // Don't HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
      throw std::logic_error("TSQR was not built with complex "
			     "arithmetic support");
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
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

  const bool testReal = params.testReal;
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
  const bool testComplex = params.testComplex;
#else // Don't HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
  const bool testComplex = false;
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX

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

  if (testReal)
    {
      TSQR_BENCHMARK_DIST_TSQR( float, "float" );
      TSQR_BENCHMARK_DIST_TSQR( double, "double" );
    }
  if (testComplex)
    {
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
      using std::complex;

      TSQR_BENCHMARK_DIST_TSQR( complex<float>, "complex<float>" );
      TSQR_BENCHMARK_DIST_TSQR( complex<double>, "complex<double>" );

#else // Don't HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
      throw std::logic_error("TSQR was not built with complex "
			     "arithmetic support");
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
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
    cmdLineProc.setOption ("field-names", 
			   &params.additionalFieldNames,
			   "Any additional field name(s) (comma-delimited "
			   "string) to add to the benchmark output.  Empty "
			   "by default.  Good for things known when invoking "
			   "the benchmark executable, but not (easily) known "
			   "inside the benchmark -- e.g., environment "
			   "variables.");
    cmdLineProc.setOption ("output-data", 
			   &params.additionalData,
			   "Any additional data to add to the output, "
			   "corresponding to the above field name(s). "
			   "Empty by default.");
    cmdLineProc.setOption ("print-field-names",
			   "no-print-field-names",
			   &params.printFieldNames,
			   "Print field names (for machine-readable output only)");
    cmdLineProc.setOption ("print-trilinos-test-stuff", 
			   "no-print-trilinos-test-stuff", 
			   &params.printTrilinosTestStuff,
			   "Print output that makes the Trilinos test "
			   "framework happy (but makes benchmark results "
			   "parsing scripts unhappy)");
    cmdLineProc.setOption ("print-matrices", 
			   "no-print-matrices", 
			   &params.printMatrices, 
			   "Print global test matrices and computed results to stderr");
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
    cmdLineProc.setOption ("real", 
			   "noreal",
			   &params.testReal,
			   "Test real arithmetic routines");
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
    cmdLineProc.setOption ("complex", 
			   "nocomplex",
			   &params.testComplex,
			   "Test complex arithmetic routines");
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
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

  if (allowedToPrint && params.printTrilinosTestStuff)
    // The Trilinos test framework expects a message like this.
    out << "\nEnd Result: TEST PASSED" << std::endl;
  return 0;
}


