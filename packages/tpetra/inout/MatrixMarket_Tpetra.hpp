//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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

#ifndef __MatrixMarket_Tpetra_hpp
#define __MatrixMarket_Tpetra_hpp

#include "MatrixMarket_Banner.hpp"
#include "MatrixMarket_CoordDataReader.hpp"
#include "MatrixMarket_util.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdexcept>

namespace MatrixMarket {
  namespace Tpetra {

    /// \class Adder
    /// \author Mark Hoemmen
    /// \brief Add entries read from file to Tpetra::CrsMatrix.
    ///
    /// This Adder assumes that entries are read in on MPI Rank 0.
    /// fillComplete() will be called later (after this Adder is done)
    /// to distribute the entries.  A more memory-scalable technique
    /// would be to buffer and distribute (i,j,Aij) triples in
    /// round-robin fashion, for each MPI rank to add.  We don't do
    /// this here.
    ///
    /// Wrap this with a SymmetrizingAdder to get the desired
    /// symmetrizing behavior.
    template<class SparseMatrixType>
    class Adder {
    public:
      typedef SparseMatrixType sparse_matrix_type;
      typedef typename SparseMatrixType::scalar_type value_type;
      /// \typedef index_type
      ///
      /// Indices of the sparse matrix are read in as global ordinals,
      /// since Matrix Market files represent the whole matrix and
      /// don't have a notion of distribution.
      typedef typename SparseMatrixType::global_ordinal_type index_type;

      /// Constructor
      /// 
      /// \param pMatrix [in/out] The sparse matrix to which to add
      ///   entries.  Currently this only needs to be nonnull and
      ///   valid on MPI Rank 0.  More advanced implementations may
      ///   choose to do more clever things for memory scalability,
      ///   though, so this should be nonnull and valid on all MPI
      ///   ranks that partake of the sparse matrix.  (It normally
      ///   would be, for a Tpetra::CrsMatrix.)
      Adder (const Teuchos::RCP<sparse_matrix_type>& pMatrix) :
	pMatrix_ (pMatrix), 
	myRank_ (Teuchos::rank (pMatrix->getComm()))
      {}

      //! Set entry (i,j) of the sparse matrix to Aij
      void 
      operator() (const index_type i, const index_type j, const value_type& Aij)
      {
	using Teuchos::tuple;
	if (myTurnToAdd())
	  {
	    TEST_FOR_EXCEPTION(pMatrix_.is_null(), std::logic_error, 
			       "The sparse matrix has not been set up yet.  "
			       "Call prepare() with the global number of rows "
			       "in the sparse matrix, before calling operator().");
	    // Matrix Market indices are 1-based, but
	    // Tpetra::CrsMatrix wants zero-based indices.
	    pMatrix_->insertGlobalValues (i-1, tuple(j-1), tuple(Aij));
	  }
      }

    private:
      //! Whether this MPI process should add the current matrix entry
      bool myTurnToAdd () const { return myRank_ == 0; }

      //! The sparse matrix to which to add entries
      Teuchos::RCP<sparse_matrix_type> pMatrix_;
      //! This MPI process' MPI rank
      int myRank_;
    };


    /// \class Reader
    /// \author Mark Hoemmen
    /// \brief Tpetra::CrsMatrix Matrix Market sparse matrix reader.
    ///
    /// The readFile() and read() class methods read in a Matrix
    /// Market "coordinate" format sparse matrix file resp. input
    /// stream (valid on MPI Rank 0), and return a Tpetra::CrsMatrix
    /// sparse matrix (SparseMatrixType) distributed in block row
    /// fashion over all the MPI ranks represented in the given
    /// communicator.
    template<class SparseMatrixType>
    class Reader {
    public:
      typedef SparseMatrixType sparse_matrix_type;
      typedef Teuchos::RCP<sparse_matrix_type> sparse_matrix_ptr;

      typedef typename SparseMatrixType::scalar_type scalar_type;
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      /// \typedef global_ordinal_type
      ///
      /// Indices of the sparse matrix are read in as global ordinals,
      /// since Matrix Market files represent the whole matrix and
      /// don't have a notion of distribution.
      typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
      typedef typename SparseMatrixType::node_type node_type;

      typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
      typedef Teuchos::RCP<map_type> map_ptr;

    private:
      /// Read in the Banner line from the given input stream
      ///
      /// \param in [in/out, valid only on Rank 0]
      /// \param lineNumber [in/out, valid only on Rank 0]
      /// \param pComm [in, global]
      ///
      /// \return Banner [non-null and valid only on Rank 0]
      static Teuchos::RCP<const Banner>
      readBanner (std::istream& in,
		  size_t& lineNumber,
		  const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
		  const bool tolerant=false)
      {
	using Teuchos::RCP;
	using Teuchos::rcp;
	const int myRank = Teuchos::rank (*pComm);

	// Non-null only on MPI Rank 0.  Using a pointer lets the data
	// persist outside the "myRank==0" scopes.
	RCP<Banner> pBanner;
	if (myRank == 0)
	  {
	    typedef Teuchos::ScalarTraits<Scalar> STS;

	    std::string line;
	    if (! in.getline(line))
	      throw std::invalid_argument ("Failed to get first (banner) line");
	    else
	      lineNumber++;

	    pBanner = rcp (new Banner (line, tolerant));
	    if (pBanner->matrixType() != "coordinate")
	      throw std::invalid_argument ("Matrix Market input file must contain a "
					   "\"coordinate\"-format sparse matrix in "
					   "order to create a Tpetra::CrsMatrix "
					   "object from it.");
	    else if (! STS::isComplex && pBanner->dataType() == "complex")
	      throw std::invalid_argument ("Matrix Market file contains complex-"
					   "valued data, but your chosen Scalar "
					   "type is real.");
	    else if (pBanner->dataType() != "real" && pBanner->dataType() != "complex")
	      throw std::invalid_argument ("Only real or complex data types (no "
					   "pattern or integer matrices) are "
					   "currently supported");
	  }
	return pBanner;
      }

      /// Call on all MPI ranks.  MPI Rank 0 attempts to read in the
      /// coordinate dimensions from the input stream.  If it
      /// succeeds, it broadcasts them to all the other MPI ranks.
      /// (All ranks need to know the matrix dimensions in order to
      /// create domain, range, and column Maps.)
      ///
      /// \return (numRows, numCols, numNonzeros)
      static Teuchos::Tuple<global_ordinal_type, 3>
      readCoordDims (std::istream& in,
		     size_t& lineNumber,
		     const Teuchos::RCP<const Banner>& pBanner,
		     const Teuchos::RCP<const comm_type>& pComm,
		     const bool tolerant = false) // ignored except on Rank 0
      {
	// Packed coordinate matrix dimensions (numRows, numCols,
	// numNonzeros); computed on Rank 0 and broadcasted to all MPI
	// ranks.
	Teuchos::Tuple<global_ordinal_type, 3> dims;

	// Read in the coordinate matrix dimensions from the input
	// stream.  "success" tells us whether reading in the
	// coordinate matrix dimensions succeeded ("Guilty unless
	// proven innocent").
	bool success = false;
	if (Teuchos::rank(*pComm) == 0)
	  { 
	    TEST_FOR_EXCEPTION(pBanner->matrixType() != "coordinate", 
			       std::invalid_argument,
			       "The Tpetra::CrsMatrix Matrix Market reader "
			       "only accepts \"coordinate\" (sparse) matrix "
			       "data.");
	    // Unpacked coordinate matrix dimensions
	    global_ordinal_type numRows, numCols, numNonzeros;
	    // Only MPI Rank 0 reads from the input stream
	    success = readCoordinateDimensions (in, numRows, numCols, 
						numNonzeros, lineNumber, 
						tolerant);
	  }
	lineNumber++;
	// Only Rank 0 did the reading, so it decides success.
	Teuchos::broadcast (*pComm, 0, &success);
	if (success)
	  {
	    // Broadcast (numRows, numCols, numNonzeros) from Rank 0
	    // to all the other MPI ranks.  We gather them up into an
	    // Tuple so we can send them with one broadcast instead of
	    // three.
	    if (Teuchos::rank(*pComm) == 0)
	      {
		dims[0] = numRows;
		dims[1] = numCols;
		dims[2] = numNonzeros;
	      }
	    Teuchos::broadcast (*pComm, 0, dims);
	    if (Teuchos::rank(*pComm) != 0)
	      {
		numRows = dims[0];
		numCols = dims[1];
		numNonzeros = dims[2];
	      }
	  }
	else
	  {
	    // Perhaps in tolerant mode, we could set all the
	    // dimensions to zero for now, and deduce correct
	    // dimensions by reading all of the file's entries and
	    // computing the max(row index) and max(column index).
	    // However, for now we just error out in that case.
	    throw std::invalid_argument ("Error reading Matrix Market sparse "
					 "matrix: failed to read coordinate "
					 "matrix dimensions.");
	  }
	return dims;
      }

      /// \class MatrixTriple
      /// \brief (domain map, range map, sparse matrix) 
      ///
      /// First two elements are the domain map and range map which
      /// are the first two arguments of CrsMatrix::fillComplete().
      /// Last element is the CrsMatrix itself, on which
      /// fillComplete() should eventually be called.  After
      /// fillComplete() is called, you can discard the maps.
      ///
      /// \note One of those annoying little definitions you have to
      /// make because you're stuck using ancient compilers that don't
      /// support useful widgets like tuple (a.k.a. tr1::tuple).
      struct MatrixTriple {
	MatrixTriple (const map_ptr& __pDomainMap, 
		      const map_ptr& __pRangeMap, 
		      const sparse_matrix_ptr& __pMatrix) : 
	  pDomainMap (__pDomainMap), 
	  pRangeMap (__pRangeMap), 
	  pMatrix (__pMatrix) {}

	map_ptr pDomainMap, pRangeMap;
	sparse_matrix_ptr pMatrix;
      };

      /// \brief Set up domain and range map and the sparse matrix.
      ///
      /// The sparse matrix is set up with the usual contiguous block
      /// row distribution.  We need the two maps when the number of
      /// rows is not equal to the number of columns.  When done
      /// adding entries to the sparse matrix, you should invoke
      /// fillComplete() with both maps (as in \c
      /// A->fillComplete(pDomainMap, pRangeMap) ).  This will do the
      /// same thing as the usual fillComplete() call with no
      /// arguments when numRows==numCols, and it will be correct
      /// otherwise.
      ///
      /// \param pComm [in] Communicator [all MPI ranks]
      /// \param pNode [in] Kokkos Node object
      /// \param dims [in, same on all MPI ranks] 
      ///   (# rows in matrix, # columns in matrix, # entries in matrix
      /// \return (domain map, range map, sparse matrix) 
      ///   [valid on all MPI ranks]
      static MatrixTriple
      setUpMatrix (const Teuchos::RCP<const Teuchos::Comm<int> >& pComm, 
		   const Teuchos::RCP<node_type>& pNode,
		   const Teuchos::Tuple<global_ordinal_type, 3>& dims)
      {
	using Teuchos::RCP;
	using Teuchos::rcp;
	typedef typename sparse_matrix_type::local_ordinal_type local_ordinal_type;
	typedef typename sparse_matrix_type::global_ordinal_type global_ordinal_type;
	typedef typename sparse_matrix_type::node_type node_type;	
	typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

	const global_ordinal_type numRows = dims[0];
	const global_ordinal_type numCols = dims[1];

	// Map to represent the range of the sparse matrix.  In our
	// case, it also represents the distribution of rows of the
	// matrix: the row Map is the same as the domain Map.  This is
	// because we are using a block row distribution, regardless
	// of whether or not the matrix is square.  The row Map is
	// always 1-1.
	RCP<map_type> pRangeMap (new map_type (numRows, 0, pComm, 
					       Tpetra::GloballyDistributed, 
					       pNode));
	// Map to represent the domain of the sparse matrix.  We only
	// need this if the matrix is not square (numRows != numCols).
	RCP<map_type> pDomainMap (new map_type (numCols, 0, pComm, 
						Tpetra::GloballyDistributed, 
						pNode));
	// (Column) Map to represent the set of columns owned by this
	// MPI rank.  This stays null if the matrix is square, since
	// we can use the default in that case.  The column Map is
	// _not_ 1-1 unless there is only 1 MPI rank.
	RCP<map_type> pColMap;
	if (numRows != numCols)
	  // Since we're setting up a block row distribution, each MPI
	  // rank owns all the columns: the local number of elements
	  // is the same as the global number of elements.
	  //
	  // FIXME (mfh 01 Feb 2011) Should check that casting numCols
	  // to size_t (the second argument) doesn't overflow.  It
	  // probably doesn't for the usual global_ordinal_type types.
	  pColMap = Tpetra::createContigMapWithNode (numCols, numCols, pComm, pNode);

	// Set up the sparse matrix.  The range Map is the same as the
	// row Map in this case, but the column Map may be different
	// (in case numRows != numCols).
	//
	// Third argument: (approx) number of nonzero entries per
	// row.  We don't know it until we read the whole Matrix
	// Market file, so we just set it to zero here and rely on the
	// (slow, but correct) default "DynamicProfile" behavior.
	sparse_matrix_ptr pMatrix = (numRows == numCols) ? 
	  rcp (new sparse_matrix_type (pRangeMap, 0)) : 
	  rcp (new sparse_matrix_type (pRangeMap, pColMap, 0));

	return MatrixTriple (pDomainMap, pRangeMap, pMatrix);
      }

      /// Return an Adder object that optionally symmetrizes the
      /// entries of the sparse matrix.
      /// 
      /// \param banner [in, nonnull and valid on Rank 0 only]
      /// \param triple [in/out, all MPI ranks]
      ///
      /// \return adder object [nonnull and valid on Rank 0 only]
      ///
      /// \note To be called only on MPI Rank 0
      typedef adder_type SymmetrizingAdder<Adder<sparse_matrix_type> >;
      static Teuchos::RCP<adder_type>
      makeAdder (const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
		 Teuchos::RCP<const Banner>& pBanner, 
		 const MatrixTriple& triple)
      {
	if (Teuchos::rank (*pComm) == 0)
	  {
	    using Teuchos::RCP;
	    using Teuchos::rcp;
	    typedef raw_adder_type Adder<sparse_matrix_type>;
	    RCP<raw_adder_type> pRaw (new raw_adder_type (triple.pMatrix));
	    return rcp (new adder_type (pRaw, pBanner->symmType()));
	  }
	else
	  return Teuchos::null;
      }

      /// Call fillComplete() on the sparse matrix using the domain
      /// and range maps in the triple.  Discard the maps and return
      /// the sparse matrix.
      ///
      /// \note To be called globally (on all MPI ranks)
      static sparse_matrix_ptr
      completeMatrix (const MatrixTriple& triple)
      {
	triple.pMatrix->fillComplete (triple.pDomainMap, triple.pRangeMap);
	return pMatrix;
      }

    public:

      static Teuchos::RCP<sparse_matrix_type>
      readFile (const std::string& filename,
		const Teuchos::RCP<const Teuchos::Comm<int> >& pComm, 
		const Teuchos::RCP<node_type>& pNode,
		const bool tolerant=false)
      {
	std::ifstream in (filename.c_str());
	return read (in, pComm, tolerant);
      }
      
      static Teuchos::RCP<sparse_matrix_type>
      read (std::istream& in,	
	    const Teuchos::RCP<const Teuchos::Comm<int> >& pComm, 
	    const Teuchos::RCP<node_type>& pNode,
	    const bool tolerant=false)
      {
	using Teuchos::RCP;
	using Teuchos::rcp;
	using Teuchos::Tuple;
	typedef Teuchos::ScalarTraits<Scalar> STS;
	const int myRank = Teuchos::rank (*pComm);

	// I _could_ write this as a giant one-liner pipe, to please
	// all you APL programmers out there.  Instead, though, I'll
	// name the output variables.
	size_t lineNumber = 1;
	RCP<const Banner> pBanner = readBanner (in, lineNumber, pComm, tolerant);
	Tuple<global_ordinal_type, 3> dims = 
	  readCoordDims (in, lineNumber, pBanner, pComm, tolerant);
	MatrixTriple triple = setUpMatrix (pComm, pNode, dims);
	RCP<adder_type> pAdder = makeAdder (pComm, pBanner, triple);

	// Read the sparse matrix entries.
	//
	// "Guilty until proven innocent."
	bool readSuccess = false;
	if (myRank == 0)
	  {
	    std::pair<bool, std::vector<size_t> > results;	    
	    // Reader for "coordinate" format sparse matrix data.
	    typedef CoordDataReader<adder_type, global_ordinal_type, scalar_type, STS::isComplex> reader_type;
	    reader_type reader (*pAdder);

	    // Read the sparse matrix entries
	    std::pair<bool, std::vector<size_t> > results = 
	      reader.read (in, lineNumber, tolerant, debug);

	    readSuccess = results.first;
	  }
	// The broadcast of success serves as a barrier.  Note that
	// Tpetra::CrsMatrix::fillComplete() only starts with a
	// barrier in debug mode, so we need some kind of barrier or
	// synchronization beforehand.
	Teuchos::broadcast (*pComm, 0, &readSuccess);

	// TODO (mfh 01 Feb 2011)
	//
	// In tolerant mode, given a "verbose" flag, report / log any
	// bad line number(s) on MPI Rank 0.
	//
	// Should we run fillComplete() in tolerant mode, if the read
	// did not succeed?
	TEST_FOR_EXCEPTION(! readSuccess, std::invalid_argument, 
			   "Failed to read in the Matrix Market sparse matrix.");
	
	// Call fillComplete() with the appropriate domain and range
	// maps, and return the sparse matrix.
	return completeMatrix (triple);
      }
    };
    
  } // namespace Tpetra
} // namespace MatrixMarket

