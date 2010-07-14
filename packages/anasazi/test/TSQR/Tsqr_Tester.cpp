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

#include "TsqrTpetraTest.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Tpetra_DefaultPlatform.hpp"

#include "Tsqr_SeqTest.hpp"
#include "Tsqr_TbbTest.hpp"
#include "Tsqr_TsqrTest.hpp"
#include "Teuchos_Time.hpp"

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

      const int numTsqrTestRoutines = 4;
      const char* TsqrTestRoutineNames[] = {"MpiSeqTSQR", "MpiTbbTSQR", "SeqTSQR", 
					    "TbbTSQR"};
      const int numSupportedKokkosNodeTypes = 2;
      const char* supportedKokkosNodeTypes[] = {"SerialNode", "TBBNode"};

      /// \class TsqrTestParameters
      /// \brief Encapsulates values of command-line parameters
      struct TsqrTestParameters {
	TsqrTestParameters () :
	  which ("MpiSeqTSQR"),
	  kokkos_node_type ("SerialNode"),
	  action (Verify), 
	  nrows (10000), ncols (10), ncores (1), ntrials (1), cache_block_size (0),
	  verbose (true), debug (false), contiguous_cache_blocks (false), 
	  human_readable (false), tpetra (false)
	{}
	std::string which, kokkos_node_type;
	TsqrTestAction action;
	int nrows, ncols, ncores, ntrials, cache_block_size;
	bool verbose, debug, contiguous_cache_blocks, human_readable, tpetra;
      };


      template< class Scalar >
      static void
      verifyTsqrAlone (RCP< const Teuchos::Comm<int> > comm,
		       const TsqrTestParameters& params)
      {
	typedef int ordinal_type;
	typedef Scalar scalar_type;
	typedef TSQR::Random::NormalGenerator< ordinal_type, scalar_type > generator_type;
	
	ordinal_type nrowsGlobal = params.nrows;
	ordinal_type ncols = params.ncols;
	int numCores = params.ncores;
	size_t cacheBlockSize = static_cast<size_t> (params.cache_block_size);
	bool contiguousCacheBlocks = params.contiguous_cache_blocks;
	bool humanReadable = params.human_readable;
	bool bDebug = params.debug;
	
	generator_type generator;
	TSQR::Trilinos::TrilinosMessenger< ordinal_type > ordinalComm (comm);
	TSQR::Trilinos::TrilinosMessenger< scalar_type > scalarComm (comm);
	
	if (params.which == "MpiSeqTSQR" || params.which == "MpiTbbTSQR")
	  {
	    using TSQR::Test::verifyTsqr;
	    verifyTsqr< ordinal_type, scalar_type, generator_type > (params.which, generator, nrowsGlobal, ncols, 
								     &ordinalComm, &scalarComm, numCores, 
								     cacheBlockSize, contiguousCacheBlocks, 
								     humanReadable, bDebug);
	  }
	else if (params.which == "TbbTSQR")
	  {
	    using TSQR::Test::verifyTbbTsqr;
	    if (ordinalComm.rank() == 0)
	      verifyTbbTsqr< ordinal_type, scalar_type, generator_type > (generator, nrowsGlobal, ncols, numCores,
									  cacheBlockSize, contiguousCacheBlocks, 
									  humanReadable, bDebug);
	    ordinalComm.barrier ();
	  }
	else if (params.which == "SeqTSQR")
	  {
	    using TSQR::Test::verifySeqTsqr;
	    if (ordinalComm.rank() == 0)
	      verifySeqTsqr< ordinal_type, scalar_type, generator_type > (generator, nrowsGlobal, ncols, cacheBlockSize, 
									  contiguousCacheBlocks, humanReadable, bDebug);
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
      parseOptions (int argc, char* argv[], const bool allowedToPrint, bool& printedHelp)
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
	  cmdLineProc.setOption ("node-type",
				 &params.kokkos_node_type,
				 "Kokkos node type to use (for Tpetra::MultiVector test)");
	  cmdLineProc.setOption ("tpetra",
				 "tsqr-alone",
				 &params.tpetra,
				 "Test Tpetra::MultiVector (default is TSQR alone)");
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
	  cmdLineProc.setOption ("ncores", 
				 &params.ncores, 
				 "Number of cores to use (per MPI process)");
	  cmdLineProc.setOption ("ntrials", 
				 &params.ntrials, 
				 "Number of trials for the benchmark");
	  cmdLineProc.setOption ("cache-block-size", 
				 &params.cache_block_size, 
				 "Cache block size (0 means set a reasonable default)");
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
	else if (params.ncores < 1)
	  throw std::domain_error ("Number of cores must be positive");
	else if (params.cache_block_size < 0)
	  throw std::domain_error ("Cache block size must be nonnegative");

	return params;
      }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

      template < class Node >
      RCP< Node > getNode (const Teuchos::ParameterList& plist) 
      {
	throw std::logic_error ("Node type not defined");
      }

      RCP< Kokkos::SerialNode > serialnode;
      template <>
      RCP< Kokkos::SerialNode > 
      getNode< Kokkos::SerialNode > (const Teuchos::ParameterList& plist) 
      {
	if (serialnode == Teuchos::null)
	  {
	    Teuchos::ParameterList plist2;
	    // Constructor requires a non-const plist, for some reason.
	    serialnode = Teuchos::rcp(new Kokkos::SerialNode(plist2));
	  }
	return serialnode;
      }

#ifdef HAVE_KOKKOS_TBB
      RCP< Kokkos::TBBNode > tbbnode;
      template <>
      RCP< Kokkos::TBBNode > 
      getNode< Kokkos::TBBNode >(const Teuchos::ParameterList& plist) 
      {
	using Teuchos::Exceptions::InvalidParameter;

	// TBBNode reads the "Num Threads" parameter, but we define
	// only the "numCores" parameter.  0, the default, tells TBB
	// to figure out on its own how many threads to use.
	std::string numCoresParamName ("numCores");
	int numCores = 0;

	// mfh 14 Jul 2010: I can never get the two-argument version
	// of plist.get() to work.  The C++ compiler is too stupid to
	// figure out which method to instantiate, even when I give it
	// all the type hints I possibly can.

	try {
	  numCores = plist.get< int > (numCoresParamName);
	} catch (InvalidParameter&) {
	  numCores = 0;
	}
	Teuchos::ParameterList tbbPlist;
	tbbPlist.set ("Num Threads", numCores);
	
	if (tbbnode == Teuchos::null)
	  tbbnode = Teuchos::rcp (new Kokkos::TBBNode (tbbPlist));
	return tbbnode;
      }
#endif

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
  using TSQR::Trilinos::Test::getNode;
  typedef RCP< const Teuchos::Comm<int> > comm_ptr;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  comm_ptr comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  size_t myRank = comm->getRank();
  //size_t numProc = comm->getSize();
  bool allowedToPrint = (myRank==0);
  bool printedHelp = false;
  
  TsqrTestParameters params = 
    parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp)
    return 0;

  if (params.tpetra)
    {
      typedef double S;
      typedef int LO;
      typedef int GO;
      typedef TSQR::ScalarTraits< S >::magnitude_type magnitude_type;

      Teuchos::ParameterList plist;
      plist.set ("cacheBlockSize", size_t(params.cache_block_size), 
		 "Cache block capacity per core (in bytes)");
      plist.set ("numCores", int(params.ncores), 
		 "If using a CPU-multithreaded intranode TSQR "
		 "implementation: number of CPU cores");
      using TSQR::Trilinos::Test::TpetraTsqrTest;

      std::pair< magnitude_type, magnitude_type > 
	results (magnitude_type(0), magnitude_type(0));
      if (params.kokkos_node_type == "SerialNode")
	{
	  typedef Kokkos::SerialNode Node;
	  RCP< Node > node = getNode< Node > (plist);
	  TpetraTsqrTest< S, LO, GO, Node > tester (params.nrows, params.ncols, 
						    comm, node, plist);
	  results = tester.getResults ();
	}
      // else if (params.kokkos_node_type == "TBBNode")
      // 	{
      // 	  typedef Kokkos::TBBNode Node;
      // 	  RCP< Node > node = getNode< Node > (plist);
      // 	  TpetraTsqrTest< S, LO, GO, Node > tester (params.nrows, params.ncols, 
      // 						    comm, node, plist);
      // 	  results = tester.getResults ();
      // 	}
      if (allowedToPrint)
	std::cout << results.first << ", " << results.second << std::endl;
    }
  else
    {
      TSQR::Trilinos::Test::verifyTsqrAlone< double > (comm, params);
      TSQR::Trilinos::Test::verifyTsqrAlone< float > (comm, params);
    }

  if (allowedToPrint) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}


