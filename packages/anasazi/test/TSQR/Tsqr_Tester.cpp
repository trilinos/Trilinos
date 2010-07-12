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

#include "Tsqr_SeqTest.hpp"
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


enum TsqrTestAction {
  Verify = 0,
  Benchmark,
  TsqrTestActionNumValues
};
TsqrTestAction TsqrTestActionValues[] = {Verify, Benchmark};
const char* TsqrTestActionNames[] = {"verify", "benchmark"};

// enum TsqrTestRoutine {
//   FullTsqr = 0,
//   FullMgs,
//   IntraNodeOnly,
//   InterNodeOnly,
//   LapackOnly,
//   TsqrTestRoutineNumValues
// };
// TsqrTestRoutine TsqrTestRoutineValues[] = 
//   {FullTsqr, FullMgs, IntraNodeOnly, InterNodeOnly, LapackOnly};
// const char* TsqrTestRoutineNames[] = {"full-tsqr", "full-mgs", "intranode-only", 
// 				      "internode-only", "lapack-only"};

const int numTsqrTestRoutines = 4;
const char* TsqrTestRoutineNames[] = {"MpiSeqTSQR", "MpiTbbTSQR", "SeqTSQR", "TbbTSQR"};

struct TsqrTestParameters {
  TsqrTestParameters () :
    which ("MpiSeqTSQR"),
    action (Verify), 
    nrows (10000), ncols (10), ncores (1), ntrials (1), cache_block_size (0),
    verbose (true), debug (false), contiguous_cache_blocks (false), human_readable (false) 
  {}

  std::string which;
  TsqrTestAction action;
  int nrows, ncols, ncores, ntrials, cache_block_size;
  bool verbose, debug, contiguous_cache_blocks, human_readable;
};


template< class Scalar >
static void
verifyTsqr (Teuchos::RCP< const Teuchos::Comm<int> > comm,
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
      TSQR::Test::verifyTsqr (params.which, generator, nrowsGlobal, ncols, 
			      &ordinalComm, &scalarComm, numCores, 
			      cacheBlockSize, contiguousCacheBlocks, 
			      humanReadable, bDebug);
    }
  else if (params.which == "SeqTSQR")
    {
      using TSQR::Test::verifySeqTsqr;
      verifySeqTsqr< ordinal_type, scalar_type, generator_type > (generator, nrowsGlobal, ncols, cacheBlockSize, 
								  contiguousCacheBlocks, humanReadable, bDebug);
    }
}




static TsqrTestParameters
parseOptions (int argc, char* argv[], const bool allowedToPrint, bool& printedHelp)
{
  using std::cerr;
  using std::endl;
  using std::string;

  printedHelp = false;

  // Command-line parameters, set to their default values.
  TsqrTestParameters params;
  try {
    Teuchos::CommandLineProcessor cmdLineProc (/* throwExceptions= */ true, /* recognizeAllOptions=*/ true);

    // mfh 08 Jul 2010
    //
    // C++ lacks sufficient introspection mechanisms for cleanly
    // expressing this iteration over options of various types as a
    // loop over previously set ParameterList entries.
    // Teuchos::ParameterList has a ConstIterator for iterating over
    // all the parameters, but it stores each parameter as a
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
			   "Whether to reorganize cache blocks into contiguous storage");
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
  
  TsqrTestParameters params = parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp)
    return 0;
  std::cerr << "FOO BAR" << std::endl;

  verifyTsqr< double > (comm, params);
  //  verifyTsqr< float > (comm, params);

  if (allowedToPrint) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}


