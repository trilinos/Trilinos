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

/// \file cxx_main_tsqr.cpp 
/// \brief Test TsqrOrthoManager and TsqrMatOrthoManager
///
/// Test the OrthoManager interface to TsqrOrthoManager and
/// TsqrMatOrthoManager, using Tpetra::MultiVector as the multivector
/// implementation.

#include "BelosConfigDefs.hpp"
#include "BelosOutputManager.hpp"
#include "BelosOrthoManagerFactory.hpp"
#include "BelosTpetraAdapter.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <iohb.h>

#include <complex>
#include <stdexcept>

#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif // HAVE_KOKKOS_TBB

using namespace Belos;
using namespace Teuchos;

using std::cout;
using std::endl;
using std::vector;

typedef double scalar_type;
typedef int local_ordinal_type;
typedef int global_ordinal_type;
typedef Kokkos::DefaultNode::DefaultNodeType node_type;

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

template< class Node >
Teuchos::RCP< Node > 
getNode() {
  TEST_FOR_EXCEPTION(true, std::logic_error, "Node type not defined.");
}

template <>
Teuchos::RCP< Kokkos::SerialNode > 
getNode< Kokkos::SerialNode >() {
  static Teuchos::RCP< Kokkos::SerialNode > serialnode;
  if (Teuchos::is_null (serialnode)) 
    {
      ParameterList pl;
      serialnode = Teuchos::rcp (new Kokkos::SerialNode (pl));
    }
  return serialnode;
}

#ifdef HAVE_KOKKOS_TBB
template <>
Teuchos::RCP< Kokkos::TBBNode > 
getNode< Kokkos::TBBNode >() {
  static Teuchos::RCP< Kokkos::TBBNode > tbbnode;
  if (Teuchos::is_null (tbbnode))
    {
      ParameterList pl;
      tbbnode = Teuchos::rcp (new Kokkos::TBBNode(pl));
    }
  return tbbnode;
}
#endif

// this is the tolerance that all tests are performed against
const magnitude_type TOL = 1.0e-12;
const magnitude_type ATOL = 10;


std::pair< Teuchos::RCP< map_type >, Teuchos::RCP< sparse_matrix_type > >
loadSparseMatrix (const Teuchos::RCP< const Teuchos::Comm<int> > comm,
		  const std::string& filename,
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
	  // format) from the file into the tuple (numRows, numCols, nnz,
	  // colptr, rowind, dvals).  The routine allocates memory for
	  // colptr, rowind, and dvals using malloc().
	  info = readHB_newmat_double (filename.c_str(), &numRows, &numCols,
				       &nnz, &colptr, &rowind, &dvals);
	  // The Harwell-Boeing routines use info == 0 to signal failure.
	  if (info != 0)
	    {
	      // rnnzmax := maximum number of nonzeros per row, over all
	      // rows of the sparse matrix.
	      std::vector<int> rnnz (numRows, 0);
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

      Teuchos::broadcast (*comm, 0, &numRows);
      Teuchos::broadcast (*comm, 0, &numCols);
      Teuchos::broadcast (*comm, 0, &rnnzmax);

      TEST_FOR_EXCEPTION(numRows != numCols, std::runtime_error,
			 "Test matrix in Harwell-Boeing sparse matrix file '" 
			 << filename << "' " << "is not square: it is " 
			 << numRows << " by " << numCols << std::endl);

      // Create Tpetra::Map to represent multivectors in the range of
      // the sparse matrix.
      pMap = Teuchos::rcp (new map_type (numRows, 0, comm, 
					 Tpetra::GloballyDistributed, 
					 getNode< node_type >()));
      // Second argument: max number of nonzero entries per row.
      pMatrix = Teuchos::rcp (new sparse_matrix_type (pMap, rnnzmax));

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
      // among the processes.
      pMatrix->fillComplete();
      debugOut << "Completed loading and distributing sparse matrix" << endl;
    } // else M == null
  else 
    {
      debugOut << "Testing with Euclidean inner product" << endl;

      // Let M remain null, and allocate map using the number of rows
      // (numRows) specified on the command line.
      pMap = Teuchos::rcp (new map_type (numRows, 0, comm, 
					 Tpetra::GloballyDistributed, 
					 getNode< node_type >()));
    }
  return std::make_pair (pMap, pMatrix);
}

////////////////////////////////////////////////////////////////////////////////
// Some forward declarations
////////////////////////////////////////////////////////////////////////////////



int 
main (int argc, char *argv[]) 
{
  using Teuchos::RCP;

  Teuchos::GlobalMPISession mpisess(&argc,&argv,&std::cout);
  RCP< const Teuchos::Comm<int> > comm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  std::ostream& finalOut = (Teuchos::rank(*comm) == 0) ? std::cout : Teuchos::oblackholestream;

  Belos::OrthoManagerFactory< scalar_type, MV, OP > factory;


  int numFailed = 0;
  bool verbose = false;
  bool debug = false;
  std::string filename;
  std::string ortho (factory.defaultName());
  int numRows = 100;
  int sizeS  = 5;
  int sizeX1 = 11; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1
  int sizeX2 = 13; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1
  bool success = true;

  // Declare an output manager for handling local output.
  // In Anasazi, this class is called BasicOutputManager.
  // In Belos, this class is called OutputManager.
  RCP< Belos::OutputManager<scalar_type> > MyOM;
  try {

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption ("verbose", "quiet", &verbose,
		    "Print messages and results.");
    cmdp.setOption ("debug", "nodebug", &debug,
		    "Print debugging information.");
    cmdp.setOption ("filename", &filename,
		    "Filename for Harwell-Boeing sparse matrix, used as the M "
		    "operator for the inner product.  If not provided, no "
		    "matrix is read and the Euclidean inner product is used.");
    {
      std::ostringstream os;
      os << "OrthoManager subclass to test.  There ";
      if (numValidOrthoManagers > 1)
	os << "are " << numValidOrthoManagers << "options: ";
      else
	os << "is " << numValidOrthoManagers << "option: ";

      os << printValidOrthoManagerList() << ".";
      cmdp.setOption ("ortho", &ortho, os.str().c_str());
    }
    cmdp.setOption ("numRows", &numRows, 
		    "Controls the number of rows of the test "
		    "multivectors.  If an input matrix is given, this "
		    "parameter\'s value is ignored.");
    cmdp.setOption ("sizeS", &sizeS, "Controls the number of columns of the "
		    "input multivector.");
    cmdp.setOption ("sizeX1", &sizeX1, "Controls the number of columns of the "
		    "first basis.");
    cmdp.setOption ("sizeX2", &sizeX2, "Controls the number of columns of the "
		    "second basis.  We require for simplicity of testing (the "
		    "routines do not require it) that sizeX1 >= sizeX2.");
    if (cmdp.parse (argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) 
      return -1;
    if (debug) 
      verbose = true;
    
    MyOM = rcp( new Belos::OutputManager<scalar_type>() );

    // Select which type(s) of messages to print
    {
      // FIXME: Calling this a "MsgType" or even an "enum MsgType"
      // confuses the compiler.
      int theType = Errors; // default (always print errors)
      if (verbose) 
        {
	  // "Verbose" also means printing out Debug messages.
	  theType = theType | Warnings | IterationDetails |
	    OrthoDetails | FinalSummary | TimingDetails |
	    StatusTestDetails | Debug;
        }
      if (debug)
        theType = theType | Debug;

      MyOM->setVerbosity (theType);
    }
    // Stream for debug output.  If debug output is not enabled, then
    // this stream doesn't print anything sent to it (it's a "black
    // hole" stream).
    std::ostream& debugOut = MyOM->stream(Debug);

    debugOut << "Belos version information:" << endl 
	     << Belos_Version() << endl << endl;

    // Map describing layout of (multi)vectors
    RCP< map_type > map;
    // Inner product operator (to be loaded)
    RCP< sparse_matrix_type > M;
    {
      // Load the inner product operator, or make it the identity
      // matrix if filename=="".  Also return an appropriate Map.
      std::pair< RCP< map_type >, RCP< sparse_matrix_type > > results = 
	loadSparseMatrix (comm, filename, debugOut);
      map = results.first;
      M = results.second;
    }

    // Instantiate the specified OrthoManager subclass for testing.
    RCP< Belos::OrthoManager< scalar_type, MV > > OM = 
      factory.makeOrthoManager (ortho, M, Teuchos::null);

    // "Prototype" multivector, from which to clone other multivectors.
    RCP< MV > S = rcp (new MV (map, sizeS));

    const int numFailed = 
      Belos::Test::runTests< scalar_type, MV > (OM, S, sizeX1, sizeX2, MyOM);

  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,finalOut,success);

  if (numFailed != 0 || ! success) 
    {
      if (numFailed != 0) {
	MyOM->stream(Errors) << numFailed << " errors." << endl;
      }
      // The Trilinos test framework depends on seeing this message,
      // so don't rely on the OutputManager to report it correctly.
      finalOut << "End Result: TEST FAILED" << endl;	
      return -1;
    }
  else 
    {
      finalOut << "End Result: TEST PASSED" << endl;
      return 0;
    }
}



