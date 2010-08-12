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
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// FIXME (mfh 04 Aug 2010) Should remove all Tpetra
// dependencies... Teuchos dependencies are necessary, of course.
#include "Tpetra_DefaultPlatform.hpp"
#include "Tsqr_TpetraMessenger.hpp"

#include "Tsqr_SeqTest.hpp"
#include "Tsqr_TsqrTest.hpp"
#include "Tsqr_CombineBenchmark.hpp"
#include "Teuchos_Time.hpp"

#ifdef HAVE_ANASAZI_COMPLEX
#  include <complex>
#endif // HAVE_ANASAZI_COMPLEX

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

      enum TsqrTestAction { Verify = 0, Benchmark, TsqrTestActionNumValues };
      TsqrTestAction TsqrTestActionValues[] = { Verify, Benchmark };
      const char* TsqrTestActionNames[] = {"verify", "benchmark"};

      const int numTsqrTestRoutines = 3;
      const char* TsqrTestRoutineNames[] = {"Combine", "MpiSeqTSQR", "SeqTSQR"};

      /// \class TsqrTestParameters
      /// \brief Encapsulates values of command-line parameters
      struct TsqrTestParameters {
	TsqrTestParameters () :
	  which ("Combine"),
	  action (Benchmark), 
	  nrows (1000), 
	  ncols (10), 
	  ntrials (10), 
	  cache_block_size (0),
	  verbose (true), 
	  debug (false),
#ifdef HAVE_ANASAZI_COMPLEX
	  test_complex_arithmetic (true),
#endif // HAVE_ANASAZI_COMPLEX
	  contiguous_cache_blocks (false), 
	  human_readable (false)
	{}

	std::string which;
	TsqrTestAction action;
	int nrows, ncols, ntrials, cache_block_size;
	bool verbose, debug;
#ifdef HAVE_ANASAZI_COMPLEX
	bool test_complex_arithmetic;
#endif // HAVE_ANASAZI_COMPLEX
	bool contiguous_cache_blocks, human_readable, tpetra;
      };

      static void
      benchmarkCombineAlone (RCP< const Teuchos::Comm<int> > comm,
			     std::ostream& out,
			     const TsqrTestParameters& params)
      {
	typedef Teuchos::Time timer_type;
	using TSQR::Test::benchmarkCombine;
	typedef int ordinal_type;

	if (comm->rank() == 0)
	  {
	    const ordinal_type numRows = params.nrows;
	    const ordinal_type numCols = params.ncols;
	    const ordinal_type numTrials = params.ntrials;
#ifdef HAVE_ANASAZI_COMPLEX
	    const bool testComplex = params.test_complex_arithmetic;
#else
	    const bool testComplex = false;
#endif // HAVE_ANASAZI_COMPLEX

	    std::vector<int> seed(4);
	    const bool useSeedValues = false;
	    benchmarkCombine< timer_type > (out, numRows, numCols, numTrials,
					    seed, useSeedValues, testComplex);
	  }
	comm->barrier();
      }

      template< class Scalar >
      static void
      verifyTsqrAlone (RCP< const Teuchos::Comm<int> > comm,
		       const TsqrTestParameters& params)
      {
	using TSQR::Random::NormalGenerator;

	typedef int ordinal_type;
	typedef Scalar scalar_type;
	typedef NormalGenerator< ordinal_type, scalar_type > generator_type;
	
	ordinal_type nrowsGlobal = params.nrows;
	ordinal_type ncols = params.ncols;
	size_t cacheBlockSize = static_cast<size_t> (params.cache_block_size);
	bool contiguousCacheBlocks = params.contiguous_cache_blocks;
#ifdef HAVE_ANASAZI_COMPLEX
	bool testComplexArithmetic = params.test_complex_arithmetic;
#else
	bool testComplexArithmetic = false;
#endif // HAVE_ANASAZI_COMPLEX
	bool humanReadable = params.human_readable;
	bool bDebug = params.debug;
	
	generator_type generator;
	//
	// FIXME (mfh 04 Aug 2010) Is this a Tpetra dependency?
	//
	TSQR::Trilinos::TpetraMessenger< ordinal_type > ordinalComm (comm);
	TSQR::Trilinos::TpetraMessenger< scalar_type > scalarComm (comm);
	
	if (params.which == "MpiSeqTSQR")
	  {
	    using TSQR::Test::verifyTsqr;
	    const int numCores = 1;

	    // 
	    // FIXME (mfh 04 Aug 2010) Remove the templating from this
	    // test.  It should just go through the usual scalar
	    // type(s), do its own random number generation, etc.
	    //
	    verifyTsqr< ordinal_type, scalar_type, generator_type > (params.which, 
								     generator, 
								     nrowsGlobal, 
								     ncols, 
								     &ordinalComm, 
								     &scalarComm, 
								     numCores, 
								     cacheBlockSize, 
								     contiguousCacheBlocks, 
								     humanReadable, 
								     bDebug);
	  }
	else if (params.which == "SeqTSQR")
	  {
	    using TSQR::Test::verifySeqTsqr;
	    if (ordinalComm.rank() == 0)
	      verifySeqTsqr (nrowsGlobal, 
			     ncols, 
			     cacheBlockSize, 
			     testComplexArithmetic,
			     contiguousCacheBlocks, 
			     humanReadable, bDebug);
	    ordinalComm.barrier ();
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
      static TsqrTestParameters
      parseOptions (int argc, 
		    char* argv[], 
		    const bool allowedToPrint, 
		    bool& printedHelp)
      {
	using std::cerr;
	using std::endl;

	printedHelp = false;

	// Command-line parameters, set to their default values.
	TsqrTestParameters params;
	try {
	  Teuchos::CommandLineProcessor cmdLineProc (/* throwExceptions= */ true, 
						     /* recognizeAllOptions=*/ true);
	  // mfh 08 Jul 2010
	  //
	  // Something I wish I could do here is represent the command-line
	  // options' values as a ParameterList, and automatically generate
	  // the CommandLineProcessor from that.  Unfortunately, C++ lacks
	  // sufficient introspection mechanisms for cleanly expressing this
	  // iteration over options of various types as a loop over
	  // previously set ParameterList entries.  Teuchos::ParameterList
	  // has a ConstIterator for iterating over all the parameters, but
	  // it stores each parameter as a Teuchos::any (wrapped in a
	  // ParameterEntry).  Teuchos::any only gives the type back as an
	  // std::type_info, and C++ can't treat that as a template
	  // parameter.  Thus, there are only two ways to express this as a
	  // loop over command-line arguments:
	  //
	  // 1. Map the std::type_info of the contents of the Teuchos::any
	  // object to a function that invokes the correct overload of
	  // setOption() 
	  //
	  // 2. Use a Smalltalk - style chain of "if type is this do that,
	  // else if type is this do that, else if...".  (Stroustrup will
	  // hate you for this.)
	  //
	  // #1 is made more difficult because std::type_info doesn't allow
	  // copy construction or assignment; it has a "const char* name()"
	  // method, though, so you could map from the type name to the
	  // appropriate handler.  #2 is bad in general, but in this case
	  // there are only three command-line argument types of interest
	  // (bool, int, std::string), so a Smalltalk-style loop would be
	  // acceptable.
	  cmdLineProc.setOption ("action", 
				 &params.action, 
				 (int) TsqrTestActionNumValues, 
				 TsqrTestActionValues,
				 TsqrTestActionNames,
				 "Which action to undertake");
	  cmdLineProc.setOption ("which", 
				 &params.which, 
				 "Which TSQR routine to test");
	  cmdLineProc.setOption ("verbose", 
				 "quiet", 
				 &params.verbose,
				 "Print messages and results");
	  cmdLineProc.setOption ("debug", 
				 "nodebug", 
				 &params.debug, 
				 "Print debugging information");
	  cmdLineProc.setOption ("nrows", 
				 &params.nrows, 
				 "Number of rows (globally) in the test matrix to factor");
	  cmdLineProc.setOption ("ncols", 
				 &params.ncols, 
				 "Number of columns in the test matrix to factor");
	  cmdLineProc.setOption ("ntrials", 
				 &params.ntrials, 
				 "Number of trials for the benchmark");
	  cmdLineProc.setOption ("cache-block-size", 
				 &params.cache_block_size, 
				 "Cache block size (0 means set a reasonable default)");
#ifdef HAVE_ANASAZI_COMPLEX
	  cmdLineProc.setOption ("complex", 
				 "nocomplex",
				 &params.test_complex_arithmetic,
				 "Test complex arithmetic");
#endif // HAVE_ANASAZI_COMPLEX
	  cmdLineProc.setOption ("contiguous-cache-blocks", 
				 "noncontiguous-cache-blocks", 
				 &params.contiguous_cache_blocks, 
				 "Reorganize cache blocks into contiguous storage");
	  cmdLineProc.setOption ("human-readable", 
				 "machine-parseable", 
				 &params.human_readable, 
				 "Make benchmark results human-readable (but hard to parse)");
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

	// Validate.  TODO (mfh 08 Jul 2010) Figure out how to do this with
	// ParameterList validators.
	if (params.nrows <= 0)
	  throw std::domain_error ("Number of rows must be positive");
	else if (params.ncols <= 0)
	  throw std::domain_error ("Number of columns must be positive");
	else if (params.nrows < params.ncols)
	  throw std::domain_error ("Number of rows must be >= number of columns");
	else if (params.cache_block_size < 0)
	  throw std::domain_error ("Cache block size must be nonnegative");

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
  using TSQR::Trilinos::Test::TsqrTestParameters;
  using TSQR::Trilinos::Test::parseOptions;
  typedef RCP< const Teuchos::Comm<int> > comm_ptr;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  comm_ptr comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  size_t myRank = comm->getRank();
  bool allowedToPrint = (myRank==0);
  bool printedHelp = false;
  
  TsqrTestParameters params = 
    parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp)
    return 0;

  if (allowedToPrint)
    TSQR::Trilinos::Test::benchmarkCombineAlone (comm, std::cout, params);
  else
    TSQR::Trilinos::Test::benchmarkCombineAlone (comm, blackhole, params);

  TSQR::Trilinos::Test::verifyTsqrAlone< double > (comm, params);
  TSQR::Trilinos::Test::verifyTsqrAlone< float > (comm, params);
#ifdef HAVE_ANASAZI_COMPLEX
  if (params.test_complex_arithmetic)
    {
      TSQR::Trilinos::Test::verifyTsqrAlone< std::complex<double> > (comm, params);
      TSQR::Trilinos::Test::verifyTsqrAlone< std::complex<float> > (comm, params);
    }
#endif // HAVE_ANASAZI_COMPLEX

  if (allowedToPrint) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}


