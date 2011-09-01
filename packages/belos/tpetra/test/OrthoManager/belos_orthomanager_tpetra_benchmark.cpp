//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2011 Sandia Corporation
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

/// \file belos_orthomanager_tpetra_benchmark.cpp
/// \brief Benchmark (Mat)OrthoManager subclass(es) with Tpetra
///
/// Benchmark various subclasses of (Mat)OrthoManager, using
/// Tpetra::MultiVector as the multivector implementation, and
/// Tpetra::Operator as the operator implementation.
///
#include "belos_orthomanager_tpetra_util.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <algorithm>

using std::endl;

//
// These typedefs make main() as generic as possible.
//
typedef double scalar_type;
typedef int local_ordinal_type;
typedef int global_ordinal_type;

#ifdef HAVE_KOKKOS_TBB
typedef Kokkos::TBBNode node_type;
#else
typedef Kokkos::SerialNode node_type;
#endif // HAVE_KOKKOS_TBB

typedef Teuchos::ScalarTraits<scalar_type> SCT;
typedef SCT::magnitudeType magnitude_type;
typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;
typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> OP;
typedef Belos::MultiVecTraits<scalar_type, MV> MVT;
typedef Belos::OperatorTraits<scalar_type, MV, OP> OPT;
typedef Teuchos::SerialDenseMatrix<int, scalar_type> serial_matrix_type;
typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> sparse_matrix_type;

/* ******************************************************************* */

/// \fn main
/// \brief Benchmark driver for (Mat)OrthoManager subclasses
int 
main (int argc, char *argv[]) 
{
  using Belos::OrthoManager;
  using Belos::OrthoManagerFactory;
  using Belos::OutputManager;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &std::cout);
  RCP<const Teuchos::Comm<int> > pComm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  // This factory object knows how to make a (Mat)OrthoManager
  // subclass, given a name for the subclass.  The name is not the
  // same as the class' syntactic name: e.g., "TSQR" is the name of
  // TsqrOrthoManager.
  OrthoManagerFactory<scalar_type, MV, OP> factory;

  // The name of the (Mat)OrthoManager subclass to instantiate.
  std::string orthoManName (factory.defaultName());

  // For SimpleOrthoManager: the normalization method to use.  Valid
  // values: "MGS", "CGS".
  std::string normalization ("CGS");

  // Name of the Harwell-Boeing sparse matrix file from which to read
  // the inner product operator matrix.  If name is "" or not provided
  // at the command line, use the standard Euclidean inner product.
  std::string filename;

  bool verbose = false; // Verbosity of output
  bool debug = false;   // Whether to print debugging-level output
  // Whether or not to run the benchmark.  If false, we let this
  // "test" pass trivially.
  bool benchmark = false; 

  // Whether to display benchmark results compactly (in a CSV format),
  // or in a human-readable table.
  bool displayResultsCompactly = false;

  // Default _local_ (per MPI process) number of rows.  This will
  // change if a sparse matrix is loaded in as an inner product
  // operator.  Regardless, the number of rows per MPI process must be
  // no less than numCols*numBlocks in order for TSQR to work.  To
  // ensure that the test always passes with default parameters, we
  // scale by the number of processes.  The default value below may be
  // changed by a command-line parameter with a corresponding name.
  int numRowsPerProcess = 100;

  // The OrthoManager is benchmarked with numBlocks multivectors of
  // width numCols each, for numTrials trials.  The values below are
  // defaults and may be changed by the corresponding command-line
  // arguments.
  int numCols = 10;
  int numBlocks = 5;
  int numTrials = 3;

  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("benchmark", "nobenchmark", &benchmark, 
		  "Whether to run the benchmark.  If not, this \"test\" "
		  "passes trivially.");
  cmdp.setOption ("verbose", "quiet", &verbose,
		  "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug,
		  "Print debugging information.");
  cmdp.setOption ("compact", "human", &displayResultsCompactly,
		  "Whether to display benchmark results compactly (in a "
		  "CSV format), or in a human-readable table.");
  cmdp.setOption ("filename", &filename,
		  "Filename of a Harwell-Boeing sparse matrix, used as the "
		  "inner product operator by the orthogonalization manager."
		  "  If not provided, no matrix is read and the Euclidean "
		  "inner product is used.");
  {
    std::ostringstream os;
    const int numValid = factory.numOrthoManagers();
    const bool plural = numValid > 1 || numValid == 0;

    os << "OrthoManager subclass to benchmark.  There ";
    os << (plural ? "are " : "is ") << numValid << (plural ? "s: " : ": ");
    factory.printValidNames (os);
    os << ".  If none is provided, the test trivially passes.";
    cmdp.setOption ("ortho", &orthoManName, os.str().c_str());
  }
  cmdp.setOption ("normalization", &normalization, 
		  "For SimpleOrthoManager (--ortho=Simple): the normalization "
		  "method to use.  Valid values: \"MGS\", \"CGS\".");
  cmdp.setOption ("numRowsPerProcess", &numRowsPerProcess, 
		  "Number of rows per MPI process in the test multivectors.  "
		  "If an input matrix is given, this value is ignored, since "
		  "the vectors must be commensurate with the dimensions of "
		  "the matrix.");
  cmdp.setOption ("numCols", &numCols, 
		  "Number of columns in the input multivector (>= 1).");
  cmdp.setOption ("numBlocks", &numBlocks, 
		  "Number of block(s) to benchmark (>= 1).");
  cmdp.setOption ("numTrials", &numTrials,
		  "Number of trial(s) per timing run (>= 1).");
  
  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult = cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // explicitly say to run the benchmark, we let this "test" pass
    // trivially.
    if (! benchmark || parseResult == CommandLineProcessor::PARSE_HELP_PRINTED)
      {
	if (Teuchos::rank(*pComm) == 0)
	  std::cout << "End Result: TEST PASSED" << endl;
	return EXIT_SUCCESS;
      }
    TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
  }
  // Total number of rows in the test vector(s).
  // This may be changed if we load in a sparse matrix.
  int numRows = numRowsPerProcess * pComm->getSize();
  //
  // Validate command-line arguments
  //
  TEST_FOR_EXCEPTION(numRowsPerProcess <= 0, std::invalid_argument, 
		     "numRowsPerProcess <= 0 is not allowed");
  TEST_FOR_EXCEPTION(numCols <= 0, std::invalid_argument, 
		     "numCols <= 0 is not allowed");
  TEST_FOR_EXCEPTION(numBlocks <= 0, std::invalid_argument, 
		     "numBlocks <= 0 is not allowed");  
    
  // Declare an output manager for handling local output.  Initialize,
  // using the caller's desired verbosity level.
  RCP<OutputManager<scalar_type> > outMan = 
    Belos::Test::makeOutputManager<scalar_type> (verbose, debug);

  // Stream for debug output.  If debug output is not enabled, then
  // this stream doesn't print anything sent to it (it's a "black
  // hole" stream).
  std::ostream& debugOut = outMan->stream(Belos::Debug);
  Belos::Test::printVersionInfo (debugOut);

  // Load the inner product operator matrix from the given filename.
  // If filename == "", use the identity matrix as the inner product
  // operator (the Euclidean inner product), and leave M as
  // Teuchos::null.  Also return an appropriate Map (which will
  // always be initialized; it should never be Teuchos::null).
  RCP<map_type> map;
  RCP<sparse_matrix_type> M; 
  {
    using Belos::Test::loadSparseMatrix;
    // If the sparse matrix is loaded successfully, this call will
    // modify numRows to be the total number of rows in the sparse
    // matrix.  Otherwise, it will leave numRows alone.
    std::pair<RCP<map_type>, RCP<sparse_matrix_type> > results = 
      loadSparseMatrix<local_ordinal_type, global_ordinal_type, node_type> (pComm, filename, numRows, debugOut);
    map = results.first;
    M = results.second;
  }
  TEST_FOR_EXCEPTION(map.is_null(), std::logic_error,
		     "Error: (Mat)OrthoManager test code failed to "
		     "initialize the Map");
  if (M.is_null())
    {
      // Number of rows per process has to be >= number of rows.
      TEST_FOR_EXCEPTION(numRowsPerProcess <= numCols,
			 std::invalid_argument,
			 "numRowsPerProcess <= numCols is not allowed");
    }
  // Loading the sparse matrix may have changed numRows, so check
  // again that the number of rows per process is >= numCols.
  // getNodeNumElements() returns a size_t, which is unsigned, and you
  // shouldn't compare signed and unsigned values.  
  if (map->getNodeNumElements() < static_cast<size_t>(numCols))
    {
      std::ostringstream os;
      os << "The number of elements on this process " << pComm->getRank() 
	 << " is too small for the number of columns that you want to test."
	 << "  There are " << map->getNodeNumElements() << " elements on "
	"this process, but the normalize() method of the MatOrthoManager "
	"subclass will need to process a multivector with " << numCols 
	 << " columns.  Not all MatOrthoManager subclasses can handle a "
	"local row block with fewer rows than columns.";
      // QUESTION (mfh 26 Jan 2011) Should this be a logic error
      // instead?  It's really TSQR's fault that it can't handle a
      // local number of elements less than the number of columns.
      throw std::invalid_argument(os.str());
    }

  // Using the factory object, instantiate the specified OrthoManager
  // subclass to be tested.  Specify "fast" parameters for a fair
  // benchmark comparison, but override the fast parameters to get the
  // desired normalization method for SimpleOrthoManaager.
  RCP<OrthoManager<scalar_type, MV> > orthoMan;
  {
    std::string label (orthoManName);
    RCP<const ParameterList> params = factory.getFastParameters (orthoManName);
    if (orthoManName == "Simple")
      {
	RCP<ParameterList> paramsCopy (new ParameterList (*params));
	paramsCopy->set ("Normalization", normalization);
	params = paramsCopy;
	label = label + " (" + normalization + " normalization)";
      }
    orthoMan = factory.makeOrthoManager (orthoManName, M, outMan, label, params);
  }

  // "Prototype" multivector.  The test code will use this (via
  // Belos::MultiVecTraits) to clone other multivectors as necessary.
  // (This means the test code doesn't need the Map, and it also makes
  // the test code independent of the idea of a Map.)  We only have to
  // allocate one column, because the entries are S are not even read.
  // (We could allocate zero columns, if the MV object allows it.  We
  // play it safe and allocate 1 column instead.)
  RCP<MV> X = rcp (new MV (map, 1));

  // "Compact" mode means that we have to override
  // TimeMonitor::summarize(), which both handles multiple MPI
  // processes correctly (only Rank 0 prints to std::cout), and prints
  // verbosely in a table form.  We deal with the former by making an
  // ostream which is std::cout on Rank 0, and prints nothing (is a
  // "bit bucket") elsewhere.  We deal with the latter inside the
  // benchmark itself.
  Teuchos::oblackholestream bitBucket;
  std::ostream& resultStream = 
    (displayResultsCompactly && Teuchos::rank(*pComm) != 0) ? bitBucket : std::cout;

  // Benchmark the OrthoManager subclass.
  typedef Belos::Test::OrthoManagerBenchmarker<scalar_type, MV> benchmarker_type;
  benchmarker_type::benchmark (orthoMan, orthoManName, normalization, X, 
			       numCols, numBlocks, numTrials, 
			       outMan, resultStream, displayResultsCompactly);
  // Only Rank 0 gets to write to cout.
  if (Teuchos::rank(*pComm) == 0)
    std::cout << "End Result: TEST PASSED" << endl;
  return EXIT_SUCCESS;
}



