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

#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Tpetra_DefaultPlatform.hpp"

#include "Tsqr_TsqrTest.hpp"
#include "Teuchos_Time.hpp"
#include "Tsqr_Random_NormalGenerator.hpp"
#include "TsqrTrilinosMessenger.hpp"

#include <sstream>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// static void
// setCmdLineArgFromParamList (Teuchos::CommandLineProcessor& cmdLineProc,
// 			    const std::string& key,
// 			    Teuchos::ParameterEntry& entry)
// {
//   using Teuchos::any_cast;

//   Teuchos::any& value = entry.getAny();
//   if (value.type() == typeof(int))
//     cmdLineProc.setOption (key, &(any_cast(value)), entry.docString().c_str());
// }


static void
verifyTsqr (const std::string& which,
	    const int nrowsGlobal,
	    const int ncols,
	    Teuchos::RCP< const Teuchos::Comm<int> > comm,
	    const int numCores,
	    const size_t cacheBlockSize,
	    const bool contiguousCacheBlocks,
	    const humanReadable,
	    const bDebug)
{
  typedef int ordinal_type;
  typedef double scalar_type;

  TSQR::Random::NormalGenerator< ordinal_type, scalar_type > generator;
  TSQR::Trilinos::TrilinosMessenger< int > ordinalComm (comm);
  TSQR::Trilinos::TrilinosMessenger< double > scalarComm (comm);
  TSQR::Test::verifyTsqr (which, generator, nrowsGlobal, ncols, &ordinalComm,
			  &scalarComm, numCores, cacheBlockSize, 
			  contiguousCacheBlocks, humanReadable, bDebug);
}

enum TsqrTestAction {
  Verify = 0,
  Benchmark,
  TsqrTestActionNumValues
};
TsqrTestAction TsqrTestActionValues[] = {Verify, Benchmark};
const char* TsqrTestActionNames[] = {"verify", "benchmark"};

enum TsqrTestRoutine {
  FullTsqr = 0,
  FullMgs,
  IntraNodeOnly,
  InterNodeOnly,
  LapackOnly,
  TsqrTestRoutineNumValues
};
TsqrTestRoutine TsqrTestRoutineValues[] = 
  {FullTsqr, FullMgs, IntraNodeOnly, InterNodeOnly, LapackOnly};
const char* TsqrTestRoutineNames[] = {"full-tsqr", "full-mgs", "intranode-only", 
				      "internode-only", "lapack-only"};

static Teuchos::ParameterList
parseOptions (int argc, char* argv[], const bool allowedToPrint, bool& printedHelp)
{
  using std::cerr;
  using std::endl;
  using std::string;

  printedHelp = false;

  // Command-line parameters, set to their default values.
  TsqrTestAction action = Verify;
  TsqrTestRoutine routine = FullTsqr; 
  bool verbose = true;
  bool debug = false;
  int nrows = 100000; 
  int ncols = 10;
  int ncores = 1;
  int ntrials = 10;
  int cache_block_size = 0;
  bool contiguous_cache_blocks = false;
  bool human_readable = false;

  Teuchos::ParameterList plist;
  plist.set ("action", action, 
	     "Which action to undertake");
  plist.set ("routine", routine,
	     "Which TSQR routine to test");
  plist.set ("verbose", verbose, 
	     "Print messages and results");
  plist.set ("debug", debug, 
	     "Print debugging information");
  plist.set ("nrows", nrows, 
	     "Number of rows (globally) in the test matrix to factor");
  plist.set ("ncols", ncols, 
	     "Number of columns in the test matrix to factor");
  plist.set ("ncores", ncores, 
	     "Number of cores to use (per MPI process)");
  plist.set ("ntrials", ntrials, 
	     "Number of trials for the benchmark");
  plist.set ("cache-block-size", cache_block_size, 
	     "Cache block size (0 means we set a reasonable default)");
  plist.set ("contiguous-cache-blocks", contiguous_cache_blocks, 
	     "Whether TSQR will reorganize cache blocks into contiguous storage");
  plist.set ("human-readable", human_readable, 
	     "Make benchmark results human-readable (but hard to parse)");
  try {
    Teuchos::CommandLineProcessor cmdLineProc (/* throwExceptions= */ true, /* recognizeAllOptions=*/ true);

    // mfh 08 Jul 2010
    //
    // C++ lacks sufficient introspection mechanisms for cleanly
    // expressing this iteration over options of various types as a
    // loop.  Teuchos::ParameterList has a ConstIterator for iterating
    // over all the parameters, but it stores each parameter as a
    // Teuchos::any (wrapped in a ParameterEntry).  Teuchos::any only
    // gives the type back as an std::type_info, and C++ can't treat
    // that as a template parameter.  Thus, there are only two ways to
    // express this as a loop over command-line arguments:
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
    cmdLineProc.setOption ("action", &action, 
			   (int) TsqrTestActionNumValues, 
			   TsqrTestActionValues,
			   TsqrTestActionNames,
			   plist.getEntry("action").docString().c_str());
    cmdLineProc.setOption ("routine", &routine, (int) TsqrTestRoutineNumValues, 
			   TsqrTestRoutineValues, TsqrTestRoutineNames,
			   plist.getEntry("routine").docString().c_str());
    cmdLineProc.setOption ("verbose", "quiet", &verbose, 
			   plist.getEntry("verbose").docString().c_str());
    cmdLineProc.setOption ("debug", "nodebug", &debug, 
			   plist.getEntry("debug").docString().c_str());
    cmdLineProc.setOption ("nrows", &nrows, 
			   plist.getEntry("nrows").docString().c_str());
    cmdLineProc.setOption ("ncols", &ncols, 
			   plist.getEntry("ncols").docString().c_str());
    cmdLineProc.setOption ("ncores", &ncores, 
			   plist.getEntry("ncores").docString().c_str());
    cmdLineProc.setOption ("ntrials", &ntrials, 
			   plist.getEntry("ntrials").docString().c_str());
    cmdLineProc.setOption ("cache-block-size", &cache_block_size, 
			   plist.getEntry("cache-block-size").docString().c_str());
    cmdLineProc.setOption ("contiguous-cache-blocks", "noncontiguous-cache-blocks", 
			   &contiguous_cache_blocks, 
			   plist.getEntry("contiguous-cache-blocks").docString().c_str());
    cmdLineProc.setOption ("human-readable", "machine-parseable", &human_readable, 
			   plist.getEntry("human-readable").docString().c_str());
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
  if (nrows <= 0)
    throw std::domain_error ("Number of rows must be positive");
  else if (ncols <= 0)
    throw std::domain_error ("Number of columns must be positive");
  else if (nrows < ncols)
    throw std::domain_error ("Number of rows must be >= number of columns");
  else if (ncores < 1)
    throw std::domain_error ("Number of cores must be positive");
  else if (cache_block_size < 0)
    throw std::domain_error ("Cache block size must be nonnegative");

  if (! printedHelp)
    {
      // Fetch the (possibly altered) values of the command-line options.
      plist.set ("action", action, 
		 "Which action to undertake");
      plist.set ("routine", routine,
		 "Which TSQR routine to test");
      plist.set ("verbose", verbose, 
		 "Print messages and results");
      plist.set ("debug", debug, 
		 "Print debugging information");
      plist.set ("nrows", nrows, 
		 "Number of rows (globally) in the test matrix to factor");
      plist.set ("ncols", ncols, 
		 "Number of columns in the test matrix to factor");
      plist.set ("ncores", ncores, 
		 "Number of cores to use (per MPI process)");
      plist.set ("ntrials", ntrials, 
		 "Number of trials for the benchmark");
      plist.set ("cache-block-size", cache_block_size, 
		 "Cache block size (0 means we set a reasonable default)");
      plist.set ("contiguous-cache-blocks", contiguous_cache_blocks, 
		 "Whether TSQR will reorganize cache blocks into contiguous storage");
      plist.set ("human-readable", human_readable, 
		 "Make benchmark results human-readable (but hard to parse)");
    }
  return plist;
}


int 
main (int argc, char *argv[]) 
{
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  Teuchos::RCP< const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  size_t myRank = comm->getRank();
  //size_t numProc = comm->getSize();
  bool allowedToPrint = (myRank==0);
  bool printedHelp = false;
  
  Teuchos::ParameterList plist = parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp)
    return 0;

  verifyTsqr ("MpiSeqTSQR", 10000, 10, comm, 1, size_t(0), false, false, false);

  if (allowedToPrint) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}


