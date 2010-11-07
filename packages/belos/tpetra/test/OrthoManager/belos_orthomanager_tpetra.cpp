//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2010 Sandia Corporation
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

/// \file belos_orthomanager_tpetra.cpp
/// \brief Test (Mat)OrthoManager subclass(es) with Tpetra
///
/// Test various subclasses of (Mat)OrthoManager, using
/// Tpetra::MultiVector as the multivector implementation,
/// and Tpetra::Operator as the operator implementation.
///
#include "BelosConfigDefs.hpp"
#include "BelosOutputManager.hpp"
#include "BelosOrthoManagerFactory.hpp"
#include "BelosOrthoManagerTest.hpp"
#include "BelosTpetraAdapter.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_DefaultNode.hpp>

// I/O for Harwell-Boeing files
#include <iohb.h>

#include <iostream>
#include <sstream>
#include <stdexcept>

using std::cout;
using std::endl;
using std::vector;

//
// These typedefs make main() as generic as possible.
//
typedef double scalar_type;
typedef int local_ordinal_type;
typedef int global_ordinal_type;
//typedef Kokkos::DefaultNode::DefaultNodeType node_type;
typedef Kokkos::SerialNode node_type;

typedef Teuchos::ScalarTraits< scalar_type > SCT;
typedef SCT::magnitudeType magnitude_type;
typedef Tpetra::MultiVector< scalar_type, local_ordinal_type, global_ordinal_type, node_type > MV;
typedef Tpetra::Operator< scalar_type, local_ordinal_type, global_ordinal_type, node_type > OP;
typedef Belos::MultiVecTraits< scalar_type, MV > MVT;
typedef Belos::OperatorTraits< scalar_type, MV, OP > OPT;
typedef Teuchos::SerialDenseMatrix< int, scalar_type > serial_matrix_type;
typedef Tpetra::Map< local_ordinal_type, global_ordinal_type, node_type > map_type;
typedef Tpetra::CrsMatrix< scalar_type, local_ordinal_type, global_ordinal_type, node_type > sparse_matrix_type;


////////////////////////////////////////////////////////////////////////////////

// The accepted way to restrict the scope of functions to their source
// file, is to use an anonymous namespace, rather than to declare the
// functions "static."  Besides, "static" is a confusingly overloaded
// term in C++.
namespace {

  void
  printVersionInfo (std::ostream& debugOut)
  {
    using std::endl;

    debugOut << "Belos version information:" << endl 
	     << Belos::Belos_Version() << endl << endl;
  }

  int
  selectVerbosity (const bool verbose, const bool debug)
  {
    // FIXME Calling this a "MsgType" (its correct type) or even an
    // "enum MsgType" confuses the compiler.
    int theType = Belos::Errors; // default (always print errors)
    if (verbose) 
      {
	// "Verbose" also means printing out Debug messages (as well
	// as everything else).
	theType = theType | 
	  Belos::Warnings | 
	  Belos::IterationDetails |
	  Belos::OrthoDetails | 
	  Belos::FinalSummary | 
	  Belos::TimingDetails |
	  Belos::StatusTestDetails | 
	  Belos::Debug;
      }
    if (debug)
      // "Debug" doesn't necessarily mean the same thing as
      // "Verbose".  We interpret "Debug" to mean printing out
      // messages marked as Debug (as well as Error messages).
      theType = theType | Belos::Debug;
    return theType;
  }

  /// \fn getNode
  /// \brief Return an RCP to a Kokkos Node
  ///
  template< class Node >
  Teuchos::RCP< Node >
  getNode() {
    throw std::runtime_error ("This Kokkos Node type not supported (compile-time error)");
  }

  template< >
  Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType >
  getNode() {
    return Kokkos::DefaultNode::getDefaultNode();
  }

  // TPINode is the default Node type, if Kokkos was built with
  // ThreadPool support.  Otherwise, if Kokkos was built with Intel
  // Threading Building Blocks support, TBBNode is the default Node
  // type.  Otherwise, SerialNode is the default.
#if defined(HAVE_KOKKOS_THREADPOOL) || defined(HAVE_KOKKOS_TBB)
  template< >
  Teuchos::RCP< Kokkos::SerialNode >
  getNode() {
    Teuchos::ParameterList defaultParams;
    return Teuchos::rcp (new Kokkos::SerialNode (defaultParams));
  }
#endif // defined(HAVE_KOKKOS_THREADPOOL) || defined(HAVE_KOKKOS_TBB)

  /// \fn loadSparseMatrix
  /// \brief Load a sparse matrix from a Harwell-Boeing file
  ///
  /// Load a sparse matrix from a Harwell-Boeing file, distribute
  /// it, and return RCPs to a map_type (the map object describing the
  /// distribution of the sparse matrix: we distribute in a way such
  /// that the domain, range, and row maps are the same) and a
  /// sparse_matrix_type (the sparse matrix itself).
  ///
  std::pair< Teuchos::RCP< map_type >, 
	     Teuchos::RCP< sparse_matrix_type > >
  loadSparseMatrix (const Teuchos::RCP< const Teuchos::Comm<int> > comm,
		    const std::string& filename,
		    int& numRows,
		    std::ostream& debugOut)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::vector;

    typedef Tpetra::Map< local_ordinal_type, global_ordinal_type, node_type > map_type;
    typedef Tpetra::CrsMatrix< scalar_type, local_ordinal_type, global_ordinal_type, node_type > sparse_matrix_type;

    const int myRank = Teuchos::rank (*comm);
    RCP< map_type > pMap;
    RCP< sparse_matrix_type > pMatrix;

    if (filename != "") 
      {
	debugOut << "Loading sparse matrix file \"" << filename << "\"" << endl;

	int loadedNumRows = 0;
	int numCols = 0;
	int nnz = -1;
	int rnnzmax = 0;
	double *dvals = NULL;
	int *colptr = NULL;
	int *rowind = NULL;

	// The Harwell-Boeing routines use info == 0 to signal failure.
	int info = 0;

	if (myRank == 0) 
	  {
	    // Proc 0 reads the sparse matrix (stored in Harwell-Boeing
	    // format) from the file into the tuple (loadedNumRows, numCols, nnz,
	    // colptr, rowind, dvals).  The routine allocates memory for
	    // colptr, rowind, and dvals using malloc().
	    info = readHB_newmat_double (filename.c_str(), &loadedNumRows, 
					 &numCols, &nnz, &colptr, &rowind, 
					 &dvals);
	    // Make sure that loadedNumRows has a sensible value,
	    // since we'll need to allocate an std::vector with that
	    // many elements.
	    TEST_FOR_EXCEPTION(loadedNumRows < 0, std::runtime_error,
			       "Harwell-Boeing sparse matrix file reports that "
			       "the matrix has # rows = " << loadedNumRows 
			       << " < 0.");

	    // The Harwell-Boeing routines use info == 0 to signal failure.
	    if (info != 0)
	      {
		// rnnzmax := maximum number of nonzeros per row, over all
		// rows of the sparse matrix.
		std::vector<int> rnnz (loadedNumRows, 0);
		for (int *ri = rowind; ri < rowind + nnz; ++ri) {
		  ++rnnz[*ri-1];
		}
		// This business with the iterator ensures that results
		// are sensible even if the sequence is empty.
		std::vector<int>::const_iterator iter = 
		  std::max_element (rnnz.begin(),rnnz.end());
		if (iter != rnnz.end())
		  rnnzmax = *iter;
		else
		  // The matrix has zero rows, so the max number of
		  // nonzeros per row is trivially zero.
		  rnnzmax = 0;
	      }
	  }

	// Proc 0 now broadcasts the sparse matrix data to the other
	// process(es).  First things broadcast are info and nnz, which
	// tell the other process(es) whether reading the sparse matrix
	// succeeded.  (info should be nonzero if so.  The
	// Harwell-Boeing routines return "C boolean true" rather than
	// the POSIX-standard "zero for success.")
	Teuchos::broadcast (*comm, 0, &info);
	Teuchos::broadcast (*comm, 0, &nnz);

	TEST_FOR_EXCEPTION(info == 0, std::runtime_error,
			   "Error reading Harwell-Boeing sparse matrix file \"" 
			   << filename << "\"" << std::endl);
	
	TEST_FOR_EXCEPTION(nnz < 0, std::runtime_error,
			   "Harwell-Boeing sparse matrix file \"" 
			   << filename << "\" reports having negative nnz "
			   << "(= " << nnz << ")"
			   << std::endl);
	
	TEST_FOR_EXCEPTION(nnz == 0, std::runtime_error,
			   "Test matrix in Harwell-Boeing sparse matrix file '" 
			   << filename << "' " << "has zero nonzero values, which "
			   << "means it does not define a valid inner product." 
			   << std::endl);

	Teuchos::broadcast (*comm, 0, &loadedNumRows);
	Teuchos::broadcast (*comm, 0, &numCols);
	Teuchos::broadcast (*comm, 0, &rnnzmax);

	TEST_FOR_EXCEPTION(loadedNumRows != numCols, std::runtime_error,
			   "Test matrix in Harwell-Boeing sparse matrix file '" 
			   << filename << "' " << "is not square: it is " 
			   << loadedNumRows << " by " << numCols << std::endl);
	// We've fully validated the number of rows, so set the
	// appropriate output parameter.
	numRows = loadedNumRows;

	// Create Tpetra::Map to represent multivectors in the range of
	// the sparse matrix.
	pMap = rcp (new map_type (numRows, 0, comm, 
				  Tpetra::GloballyDistributed,
				  getNode< node_type >()));
	// Second argument: max number of nonzero entries per row.
	pMatrix = rcp (new sparse_matrix_type (pMap, rnnzmax));

	if (myRank == 0) 
	  {
	    // Convert from Harwell-Boeing format (compressed sparse
	    // column, one-indexed) to CrsMatrix format (compressed
	    // sparse row, zero-index).  We do this by iterating over
	    // all the columns of the matrix.
	    int curNonzeroIndex = 0;
	    for (int c = 0; c < numCols; ++c) 
	      {
		for (int colnnz = 0; colnnz < colptr[c+1] - colptr[c]; ++colnnz) 
		  {
		    // Row index: *rptr - 1 (1-based -> 0-based indexing)
		    // Column index: c
		    // Value to insert there: *dptr
		    const int curGlobalRowIndex = rowind[curNonzeroIndex] - 1;
		    const scalar_type curValue = dvals[curNonzeroIndex];
		    pMatrix->insertGlobalValues (curGlobalRowIndex, 
						 Teuchos::tuple(c), 
						 Teuchos::tuple(curValue));
		    curNonzeroIndex++;
		  }
	      }
	  }
	if (myRank == 0) 
	  {
	    // Free memory allocated by the Harwell-Boeing input routine.
	    if (dvals != NULL)
	      {
		free (dvals);
		dvals = NULL;
	      }
	    if (colptr != NULL)
	      {
		free (colptr);
		colptr = NULL;
	      }
	    if (rowind != NULL)
	      {
		free (rowind);
		rowind = NULL;
	      }
	  }
	// We're done reading in the sparse matrix.  Now distribute it
	// among the processes.  The domain, range, and row maps are
	// the same (the matrix must be square).
	pMatrix->fillComplete();
	debugOut << "Completed loading and distributing sparse matrix" << endl;
      } // else M == null
    else 
      {
	debugOut << "Testing with Euclidean inner product" << endl;

	// Let M remain null, and allocate map using the number of rows
	// (numRows) specified on the command line.
	pMap = rcp (new map_type (numRows, 0, comm, 
				  Tpetra::GloballyDistributed, 
				  getNode< node_type >()));
      }
    return std::make_pair (pMap, pMatrix);
  }
}


/// \fn main
/// \brief Test driver for (Mat)OrthoManager subclasses
int 
main (int argc, char *argv[]) 
{
  using Teuchos::CommandLineProcessor;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpisess(&argc,&argv,&std::cout);
  RCP< const Teuchos::Comm<int> > comm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  // This factory object knows how to make a (Mat)OrthoManager
  // subclass, given a name for the subclass.  The name is not the
  // same as the class' syntactic name: e.g., "DKGS" is the name of
  // DkgsOrthoManager.
  Belos::OrthoManagerFactory< scalar_type, MV, OP > factory;
  // The name of the (Mat)OrthoManager subclass to instantiate.
  std::string ortho (factory.defaultName());

  // Name of the Harwell-Boeing sparse matrix file from which to read
  // the inner product operator matrix.  If name is "" or not provided
  // at the command line, use the standard Euclidean inner product.
  std::string filename;

  bool verbose = false;
  bool debug = false;
  int numRows = 100;
  int sizeS  = 5;
  int sizeX1 = 11; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1
  int sizeX2 = 13; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption ("verbose", "quiet", &verbose,
		  "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug,
		  "Print debugging information.");
  cmdp.setOption ("filename", &filename,
		  "Filename of a Harwell-Boeing sparse matrix, used as the "
		  "inner product operator by the orthogonalization manager."
		  "  If not provided, no matrix is read and the Euclidean "
		  "inner product is used.");
  {
    std::ostringstream os;
    const int numValid = factory.numOrthoManagers();
    const bool plural = numValid > 1 || numValid == 0;

    os << "OrthoManager subclass to test.  There ";
    os << (plural ? "are " : "is ") << numValid << (plural ? "s: " : ": ");
    factory.printValidNames (os);
    os << ".";
    cmdp.setOption ("ortho", &ortho, os.str().c_str());
  }
  cmdp.setOption ("numRows", &numRows, 
		  "Controls the number of rows of the test "
		  "multivectors.  If an input matrix is given, this "
		  "parameter\'s value is ignored, since the vectors must "
		  "be commensurate with the dimensions of the matrix.");
  cmdp.setOption ("sizeS", &sizeS, "Controls the number of columns of the "
		  "input multivector.");
  cmdp.setOption ("sizeX1", &sizeX1, "Controls the number of columns of the "
		  "first basis.");
  cmdp.setOption ("sizeX2", &sizeX2, "Controls the number of columns of the "
		  "second basis.  We require for simplicity of testing (the "
		  "routines do not require it) that sizeX1 >= sizeX2.");
  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult = cmdp.parse (argc,argv);
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED)
      {
	if (Teuchos::rank(*comm) == 0)
	  std::cout << "End Result: TEST PASSED" << endl;
	return EXIT_SUCCESS;
      }
    TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
  }
  //
  // Validate command-line arguments
  //
  TEST_FOR_EXCEPTION(numRows <= 0, std::invalid_argument, "numRows <= 0 is not allowed");
  TEST_FOR_EXCEPTION(numRows <= sizeS + sizeX1 + sizeX2, std::invalid_argument, 
		     "numRows <= sizeS + sizeX1 + sizeX2 is not allowed");
    
  // Declare an output manager for handling local output.  Initialize,
  // using the caller's desired verbosity level.
  //
  // FIXME In Anasazi, this class is called BasicOutputManager.  In
  // Belos, this class is called OutputManager.  The difference is
  // annoying and should be factored out.
  RCP< Belos::OutputManager< scalar_type > > MyOM (new Belos::OutputManager< scalar_type > (selectVerbosity (verbose, debug)));

  // Stream for debug output.  If debug output is not enabled, then
  // this stream doesn't print anything sent to it (it's a "black
  // hole" stream).
  std::ostream& debugOut = MyOM->stream(Belos::Debug);
  printVersionInfo (debugOut);

  // Load the inner product operator matrix from the given filename.
  // If filename == "", use the identity matrix as the inner product
  // operator (the Euclidean inner product), and leave M as
  // Teuchos::null.  Also return an appropriate Map (which will
  // always be initialized; it should never be Teuchos::null).
  RCP< map_type > map;
  RCP< sparse_matrix_type > M; 
  {
    // If the sparse matrix is loaded successfully, this call will
    // modify numRows to be the number of rows in the sparse matrix.
    // Otherwise, it will leave numRows alone.
    std::pair< RCP< map_type >, RCP< sparse_matrix_type > > results = 
      loadSparseMatrix (comm, filename, numRows, debugOut);
    map = results.first;
    M = results.second;
  }
  TEST_FOR_EXCEPTION(Teuchos::is_null(map), std::logic_error,
		     "Error: (Mat)OrthoManager test code failed to "
		     "initialize the Map");

  // Using the factory object, instantiate the specified
  // OrthoManager subclass to be tested.
  //
  // TODO Read in a ParameterList for configuring the specific
  // OrthoManager subclass.

  RCP< Belos::OrthoManager< scalar_type, MV > > OM = 
    factory.makeOrthoManager (ortho, M, Teuchos::null);

  // Whether the specific OrthoManager subclass promises to compute
  // rank-revealing orthogonalizations.  If yes, then test it on
  // rank-deficient multivectors, otherwise only test it on full-rank
  // multivectors.
  const bool isRankRevealing = factory.isRankRevealing (ortho);

  // "Prototype" multivector.  The test code will use this (via
  // Belos::MultiVecTraits) to clone other multivectors as
  // necessary.  (This means the test code doesn't need the Map, and
  // it also makes the test code independent of the idea of a Map.)
  RCP< MV > S = rcp (new MV (map, sizeS));

  // Test the OrthoManager subclass.  Return the number of tests
  // that failed.  None of the tests should fail (this function
  // should return zero).
  int numFailed = 0;
  {
    typedef Belos::Test::OrthoManagerTester< scalar_type, MV > tester_type;
    numFailed = tester_type::runTests (OM, isRankRevealing, S, 
				       sizeX1, sizeX2, MyOM);
  }

  // Only Rank 0 gets to write to cout.  The other processes dump
  // output to a black hole.
  //std::ostream& finalOut = (Teuchos::rank(*comm) == 0) ? std::cout : Teuchos::oblackholestream;

  if (numFailed != 0)
    {
      MyOM->stream(Belos::Errors) << numFailed << " errors." << endl;

      // The Trilinos test framework depends on seeing this message,
      // so don't rely on the OutputManager to report it correctly.
      if (Teuchos::rank(*comm) == 0)
	std::cout << "End Result: TEST FAILED" << endl;	
      return EXIT_FAILURE;
    }
  else 
    {
      if (Teuchos::rank(*comm) == 0)
	std::cout << "End Result: TEST PASSED" << endl;
      return EXIT_SUCCESS;
    }
}



