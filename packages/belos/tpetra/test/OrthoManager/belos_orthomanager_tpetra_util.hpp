//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
#ifndef __belos_orthomanager_tpetra_util_hpp
#define __belos_orthomanager_tpetra_util_hpp

#include <BelosTypes.hpp> // includes BelosConfigDefs.hpp
#include <BelosOrthoManagerFactory.hpp>
#include <BelosOrthoManagerTest.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_DefaultNode.hpp>

// I/O for Harwell-Boeing files
#include <iohb.h>

namespace Belos {
  namespace Test {

    //! Print Belos version information
    void
    printVersionInfo (std::ostream& debugOut)
    {
      using std::endl;
      
      debugOut << "Belos version information:" << endl 
	       << Belos::Belos_Version() << endl << endl;
    }

    //! Return a MsgType enum to specify Belos::OutputManager verbosity
    int
    selectVerbosity (const bool verbose, const bool debug)
    {
      // NOTE Calling this a "MsgType" (its correct type) or even an
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

    template<class Scalar>
    Teuchos::RCP<Belos::OutputManager<Scalar> > 
    makeOutputManager (const bool verbose, const bool debug)
    {
      return Teuchos::rcp (new Belos::OutputManager<Scalar> (selectVerbosity (verbose, debug)));
    }

    /// \fn getNode
    /// \brief Return an RCP to a Kokkos Node
    ///
    template<class NodeType>
    Teuchos::RCP<NodeType>
    getNode() {
      throw std::runtime_error ("This Kokkos Node type not supported (compile-time error)");
    }

    template<>
    Teuchos::RCP<Kokkos::SerialNode>
    getNode() {
      Teuchos::ParameterList defaultParams;
      return Teuchos::rcp (new Kokkos::SerialNode (defaultParams));
    }

#if defined(HAVE_KOKKOSCLASSIC_TBB)
    template<>
    Teuchos::RCP<Kokkos::TBBNode>
    getNode() {
      // "Num Threads" specifies the number of threads.  Defaults to an
      // automatically chosen value.
      Teuchos::ParameterList defaultParams;
      return Teuchos::rcp (new Kokkos::TBBNode (defaultParams));
    }
#endif // defined(HAVE_KOKKOSCLASSIC_TBB)

    /// \fn loadSparseMatrix
    /// \brief Load a sparse matrix from a Harwell-Boeing file
    ///
    /// Load a sparse matrix from a Harwell-Boeing file, distribute
    /// it, and return RCPs to a map_type (the map object describing the
    /// distribution of the sparse matrix: we distribute in a way such
    /// that the domain, range, and row maps are the same) and a
    /// sparse_matrix_type (the sparse matrix itself).
    ///
    template<class LO, class GO, class NodeType>
    std::pair<Teuchos::RCP<Tpetra::Map<LO, GO, NodeType> >, Teuchos::RCP<Tpetra::CrsMatrix<double, LO, GO, NodeType> > >
    loadSparseMatrix (const Teuchos::RCP< const Teuchos::Comm<int> > pComm,
		      const std::string& filename,
		      int& numRows,
		      std::ostream& debugOut)
    {
      typedef double scalar_type;
      typedef LO local_ordinal_type;
      typedef GO global_ordinal_type;
      typedef NodeType node_type;

      typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
      typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> sparse_matrix_type;

      using Teuchos::RCP;
      using Teuchos::rcp;
      using std::vector;

      const int myRank = Teuchos::rank (*pComm);
      RCP<map_type> pMap;
      RCP<sparse_matrix_type> pMatrix;

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
	      TEUCHOS_TEST_FOR_EXCEPTION(loadedNumRows < 0, std::runtime_error,
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
	  Teuchos::broadcast (*pComm, 0, &info);
	  Teuchos::broadcast (*pComm, 0, &nnz);

	  TEUCHOS_TEST_FOR_EXCEPTION(info == 0, std::runtime_error,
			     "Error reading Harwell-Boeing sparse matrix file \"" 
			     << filename << "\"" << std::endl);
	
	  TEUCHOS_TEST_FOR_EXCEPTION(nnz < 0, std::runtime_error,
			     "Harwell-Boeing sparse matrix file \"" 
			     << filename << "\" reports having negative nnz "
			     << "(= " << nnz << ")"
			     << std::endl);
	
	  TEUCHOS_TEST_FOR_EXCEPTION(nnz == 0, std::runtime_error,
			     "Test matrix in Harwell-Boeing sparse matrix file '" 
			     << filename << "' " << "has zero nonzero values, which "
			     << "means it does not define a valid inner product." 
			     << std::endl);

	  Teuchos::broadcast (*pComm, 0, &loadedNumRows);
	  Teuchos::broadcast (*pComm, 0, &numCols);
	  Teuchos::broadcast (*pComm, 0, &rnnzmax);

	  TEUCHOS_TEST_FOR_EXCEPTION(loadedNumRows != numCols, std::runtime_error,
			     "Test matrix in Harwell-Boeing sparse matrix file '" 
			     << filename << "' " << "is not square: it is " 
			     << loadedNumRows << " by " << numCols << std::endl);
	  // We've fully validated the number of rows, so set the
	  // appropriate output parameter.
	  numRows = loadedNumRows;

	  // Create Tpetra::Map to represent multivectors in the range of
	  // the sparse matrix.
	  pMap = rcp (new map_type (numRows, 0, pComm, 
				    Tpetra::GloballyDistributed,
				    getNode<node_type>()));
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
	  pMap = rcp (new map_type (numRows, 0, pComm, 
				    Tpetra::GloballyDistributed, 
				    getNode<node_type>()));
	}
      return std::make_pair (pMap, pMatrix);
    }

  } // namespace Test
} // namespace Belos

#endif // __belos_orthomanager_tpetra_util_hpp
