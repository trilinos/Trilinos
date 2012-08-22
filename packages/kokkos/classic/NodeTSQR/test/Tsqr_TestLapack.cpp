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

#include "Tsqr_ConfigDefs.hpp"
#include "Teuchos_ConfigDefs.hpp" // HAVE_MPI
#include "Teuchos_Tuple.hpp"
#ifdef HAVE_MPI
#  include "Teuchos_GlobalMPISession.hpp"
#  include "Teuchos_oblackholestream.hpp"
#endif // HAVE_MPI
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Tsqr_SeqTest.hpp"

#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
#  include <complex>
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX

#include <sstream>
#include <stdexcept>
#include <vector>


namespace TSQR { 
  namespace Trilinos { 
    namespace Test {

      const char docString[] = "This program compares LAPACK\'s QR factorization"
	" (with TSQR).  Accuracy and performance tests are included.";

      using Teuchos::RCP;
      using Teuchos::Tuple;

      /// \class LapackTestParameters
      /// \brief Encapsulates values of command-line parameters
      ///
      struct LapackTestParameters {
	LapackTestParameters () :
	  verify (false),
	  benchmark (false),
	  numRows (1000),
	  numCols (10),  
	  numTrials (10),
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	  testComplex (true),
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	  printFieldNames (true),
	  printTrilinosTestStuff (true),
	  humanReadable (false),
	  debug (false)
	{}

	bool verify, benchmark;
	int numRows, numCols, numTrials;
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	bool testComplex;
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	std::string additionalFieldNames, additionalData;
	bool printFieldNames, printTrilinosTestStuff, humanReadable, debug;
      };

      static void
      benchmark (std::ostream& out,
		 const LapackTestParameters& params)
      {
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	const bool testComplex = params.testComplex;
#else
	const bool testComplex = false;
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX

	using TSQR::Test::benchmarkLapack;
	benchmarkLapack (out, 
			 params.numRows,
			 params.numCols, 
			 params.numTrials, 
			 testComplex, 
			 params.additionalFieldNames,
			 params.additionalData,
			 params.printFieldNames,
			 params.humanReadable);
      }

      static void
      verify (std::ostream& out, 
	      const LapackTestParameters& params)
      {
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	const bool testComplex = params.testComplex;
#else
	const bool testComplex = false;
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX

	using TSQR::Test::verifyLapack;
	verifyLapack (out, 
		      params.numRows, 
		      params.numCols, 
		      testComplex,
		      params.additionalFieldNames,
		      params.additionalData,
		      params.printFieldNames, 
		      params.humanReadable, 
		      params.debug);
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
      static LapackTestParameters
      parseOptions (int argc, 
		    char* argv[], 
		    const bool allowedToPrint, 
		    bool& printedHelp)
      {
	using std::cerr;
	using std::endl;

	printedHelp = false;

	// Command-line parameters, set to their default values.
	LapackTestParameters params;

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
	  cmdLineProc.setOption ("nrows", 
				 &params.numRows, 
				 "Number of rows in the test matrix");
	  cmdLineProc.setOption ("ncols", 
				 &params.numCols, 
				 "Number of columns in the test matrix");
	  cmdLineProc.setOption ("ntrials", 
				 &params.numTrials, 
				 "Number of trials (only used when \"--benchmark\"");
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	  cmdLineProc.setOption ("complex", 
				 "nocomplex",
				 &params.testComplex,
				 "Test complex arithmetic, as well as real");
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
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
				 "Print field names for benchmark output (including "
				 "any arguments to --field-names).");
	  cmdLineProc.setOption ("print-trilinos-test-stuff", 
				 "no-print-trilinos-test-stuff", 
				 &params.printTrilinosTestStuff,
				 "Print output that makes the Trilinos test "
				 "framework happy (but makes benchmark results "
				 "parsing scripts unhappy)");
	  cmdLineProc.setOption ("human-readable",
				 "machine-readable",
				 &params.humanReadable,
				 "If set, make output easy to read by humans "
				 "(but hard to parse)");
	  cmdLineProc.setOption ("debug", 
				 "nodebug", 
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
	if (params.numRows <= 0)
	  throw std::invalid_argument ("Number of rows must be positive");
	else if (params.numCols <= 0)
	  throw std::invalid_argument ("Number of columns must be positive");
	else if (params.numRows < params.numCols)
	  throw std::invalid_argument ("Number of rows must be >= number of columns");
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
  using TSQR::Trilinos::Test::LapackTestParameters;
  using TSQR::Trilinos::Test::parseOptions;
  using std::endl;

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
  const bool performingTests = (myRank == 0);
  const bool allowedToPrint = (myRank == 0);

#else // Don't HAVE_MPI: single-node test

  const bool performingTests = true;
  const bool allowedToPrint = true;
  std::ostream& out = std::cout;
#endif // HAVE_MPI

  // Fetch command-line parameters.
  bool printedHelp = false;
  LapackTestParameters params = 
    parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp)
    return 0;

  if (performingTests)
    {
      if (params.benchmark)
	TSQR::Trilinos::Test::benchmark (out, params);

      // We allow the same run to do both benchmark and verify.
      if (params.verify)
	TSQR::Trilinos::Test::verify (out, params);

      if (params.printTrilinosTestStuff)
	// The Trilinos test framework expects a message like this.
	out << "\nEnd Result: TEST PASSED" << endl;
    }
  return 0;
}


