// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
#ifndef __Tpetra_MultiVectorFiller_hpp
#define __Tpetra_MultiVectorFiller_hpp

#include <Tpetra_MultiVector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <iterator>

namespace {
  /// \class MultiVectorFillerData
  /// \brief Implementation of fill and local assembly for \c MultiVectorFiller.
  /// \author Mark Hoemmen
  ///
  /// \tparam MV Specialization of \c Tpetra::MultiVector.
  template<class MV>
  class MultiVectorFillerData {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;

    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    /// \brief Default constructor (sets number of columns to zero).
    ///
    /// \param map [in] Map, which this object may or may not use as a
    ///   hint to separate local from nonlocal data.  This need not be
    ///   the same Map as that of the multivector output of \c
    ///   globalAssemble().
    ///
    /// Before using this object, you should call \c setNumColumns()
    /// to set the number of columns in the output multivector.
    /// Otherwise, the two-argument version of \c
    /// sumIntoGlobalValues() won't actually do anything.
    MultiVectorFillerData (const Teuchos::RCP<const map_type>& map) : 
      map_ (map),
      numCols_ (0)
    {}

    /// \brief Constructor.
    ///
    /// \param map [in] Map, which this object may or may not use as a
    ///   hint to separate local from nonlocal data.  This need not be
    ///   the same Map as that of the multivector output of \c
    ///   globalAssemble().
    ///
    /// \param numColumns [in] The (expected) number of columns in the
    ///   output multivector.  You can always change this later by
    ///   calling \c setNumColumns().
    ///
    /// \note If the number of columns given here is not the same as
    ///   the number of columns in the output multivector, you should
    ///   call \c setNumColumns() first before inserting any data.
    ///   Otherwise, the two-argument version of \c
    ///   sumIntoGlobalValues() won't do the right thing.
    MultiVectorFillerData (const Teuchos::RCP<const map_type>& map,
			   const size_t numColumns) : 
      map_ (map),
      numCols_ (numColumns),
      sourceIndices_ (numColumns),
      sourceValues_ (numColumns)
    {}

    //! Set the number of columns in the output multivector.
    void
    setNumColumns (const size_t newNumColumns) 
    {
      const size_t oldNumColumns = getNumColumns();
      if (newNumColumns >= oldNumColumns) {
	for (size_t j = oldNumColumns; j < newNumColumns; ++j) {
	  sourceIndices_.push_back (Teuchos::Array<global_ordinal_type> ());
	  sourceValues_.push_back (Teuchos::Array<scalar_type> ());
	}
      } 
      else {
	// This may not necessarily deallocate any data, but that's OK.
	sourceIndices_.resize (newNumColumns);
	sourceValues_.resize (newNumColumns);
      }
      numCols_ = oldNumColumns;
    }

    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 size_t column,
			 Teuchos::ArrayView<const scalar_type> values)
    {
      if (column >= getNumColumns()) {
	for (size_t j = column; j < getNumColumns(); ++j) {
	  sourceIndices_.push_back (Teuchos::Array<global_ordinal_type> ());
	  sourceValues_.push_back (Teuchos::Array<scalar_type> ());
	}
      }
      std::copy (rows.begin(), rows.end(), std::back_inserter (sourceIndices_[column]));
      std::copy (values.begin(), values.end(), std::back_inserter (sourceValues_[column]));
    }

    /// Data for each column are stored contiguously in rows and in
    /// values.  Thus, rows and values are in rowwise order, even
    /// though they may be stored in columnwise order in the
    /// multivector.
    ///
    /// Be sure that the number of columns is set correctly before
    /// calling this.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 Teuchos::ArrayView<const scalar_type> values)
    {
      typedef global_ordinal_type GO;
      typedef scalar_type ST;
      typedef typename Teuchos::ArrayView<const GO>::const_iterator GoIter;
      typedef typename Teuchos::ArrayView<const ST>::const_iterator StIter;

      const size_t numColumns = getNumColumns();
      GoIter rowIter = rows.begin();
      StIter valIter = values.begin();
      for (size_t j = 0; j < numColumns; ++j) {
	GoIter rowIterNext = rowIter + numColumns;
	StIter valIterNext = valIter + numColumns;
	std::copy (rowIter, rowIterNext, std::back_inserter (sourceIndices_[j]));
	std::copy (valIter, valIterNext, std::back_inserter (sourceValues_[j]));
	rowIter = rowIterNext;
	valIter = valIterNext;
      }
    }

    /// \brief Locally assemble into X.
    ///
    /// \param X [in/out] Multivector (overlapping source distribution).
    ///
    /// \param f [in/out] Binary function that defines the combine
    ///   mode.  It must define scalar_type operator (const
    ///   scalar_type&, const scalar_type&).  It need not necessarily
    ///   be commutative or even associative, but it should be
    ///   thread-safe in case we decide to parallelize local assembly.
    ///   We call it via X(i,j) = f(X(i,j), Y(i,j)), so write your
    ///   possibly nonassociative or noncommutative operation
    ///   accordingly.
    ///
    /// X is distributed by the source Map (with possible overlap) of
    /// the Export operation.  The source Map of the Export includes
    /// both the elements owned by this object's constructor's input
    /// Map, and the indices inserted by \c sumIntoGlobalValues().
    ///
    /// Precondition: The set of global indices in X's Map equals the
    /// the union of the entries of nonlocalIndices_[j] for all valid
    /// columns j.
    ///
    /// \note You can get the usual ADD combine mode by supplying f =
    ///   std::plus<scalar_type>.
    template<class BinaryFunction>
    void
    locallyAssemble (MV& X, BinaryFunction& f)
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;
      typedef node_type NT;
      typedef Tpetra::Map<LO, GO, NT> map_type;

      RCP<const map_type> map = X.getMap();
      Array<LO> localIndices;
      const size_t numColumns = getNumColumns();
      for (size_t j = 0; j < numColumns; ++j) {
	const typename Array<const GO>::size_type numIndices = sourceIndices_[j].size();
	// Precompute all the local indices before saving to the
	// vector.  Hopefully this gives us a little bit of locality
	// in the global->local conversion, at the expense of a little
	// more storage.
	if (localIndices.size() < numIndices) {
	  localIndices.resize (numIndices);
	}
	ArrayView<LO> localIndicesView = localIndices.view (0, numIndices);
	ArrayView<const GO> globalIndicesView = sourceIndices_[j].view (0, numIndices);
	for (typename ArrayView<const GO>::size_type i = 0; i < numIndices; ++i) {
	  localIndices[i] = map->getLocalElement (globalIndicesView[i]);
	}

	ArrayRCP<ST> X_j = X.getDataNonConst (j);
	ArrayView<const ST> localValues = sourceValues_[j].view (0, numIndices);
	for (typename ArrayView<const GO>::size_type i = 0; i < numIndices; ++i) {
	  const LO localInd = localIndices[i];
	  X_j[localInd] = f (X_j[localInd], localValues[i]);
	}
      }
    }

    //! \c locallyAssemble() for the usual ADD combine mode.
    void 
    locallyAssemble (MV& X)
    {
      std::plus<double> f;
      locallyAssemble<std::plus<scalar_type> > (X, f);
    }

    //! Clear the contents of the vector, making it implicitly a vector of zeros.
    void clear() {
      Teuchos::Array<Teuchos::Array<global_ordinal_type> > newSourceIndices;
      Teuchos::Array<Teuchos::Array<scalar_type> > newSourceValues;
      // The standard STL idiom for clearing the contents of a vector completely.
      std::swap (sourceIndices_, newSourceIndices);
      std::swap (sourceValues_, newSourceValues);
    }

    //! All source indices (local and nonlocal) of the source Map, sorted and unique.
    Teuchos::Array<global_ordinal_type> 
    getSourceIndices () const 
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::as;
      typedef global_ordinal_type GO;
      typedef typename Array<GO>::iterator iter_type;
      typedef typename Array<GO>::size_type size_type;

      Array<GO> allInds (0); // will resize below
      const size_t numCols = getNumColumns();

      if (numCols == 1) {
	// Special case for 1 column avoids copying indices.  Pick the
	// size of the array exactly the first time so there are at
	// most two allocations (the constructor may choose to
	// allocate).
	allInds.resize (sourceIndices_[0].size());
	std::copy (sourceIndices_[0].begin(), sourceIndices_[0].end(),
		   allInds.begin());
	std::sort (allInds.begin(), allInds.end());
	iter_type it = std::unique (allInds.begin(), allInds.end());
	const size_type curSize = as<size_type> (it - allInds.begin());
	allInds.resize (curSize);
      }
      else {
	// Carefully collect all the row indices one column at a time.
	// This ensures that the total allocation size in this routine
	// is independent of the number of columns.  Also, only sort
	// the current column's indices.  Use merges to ensure sorted
	// order in the collected final result.
	Array<GO> curInds;
	for (size_t j = 0; j < numCols; ++j) {
	  // Collect the current column's source indices into curInds.
	  // Sort it and make the entries unique.
	  curInds.resize (sourceIndices_[j].size());
	  std::copy (sourceIndices_[j].begin(), sourceIndices_[j].end(),
		     curInds.begin());
	  std::sort (curInds.begin(), curInds.end());
	  iter_type it = std::unique (curInds.begin(), curInds.end());
	  const size_type curSize = as<size_type> (it - curInds.begin());
	  curInds.resize (curSize);

	  if (curSize > 0) {
	    // Assume inductively that the entries of allInds are
	    // sorted and unique.  Make more space in allInds for the
	    // new entries (since std::inplace_merge can only work in
	    // one array), copy them in, merge, and make the results
	    // unique.
	    const size_type oldSize = allInds.size();
	    allInds.resize (oldSize + curSize);

	    //iter_type middle = &allInds[oldSize];
	    iter_type middle;
	    {
	      ArrayView<GO> oldAllInds = allInds.view (oldSize, curSize);
	      middle = oldAllInds.begin();
	    }
	    std::copy (curInds.begin(), curInds.end(), middle);
	    std::inplace_merge (allInds.begin(), middle, allInds.end());
	    iter_type it2 = std::unique (allInds.begin(), allInds.end());
	    const size_type newSize = as<size_type> (it2 - allInds.begin());
	    allInds.resize (newSize);
	  }
	}
      }
      return allInds;
    }

  private:
    Teuchos::RCP<const map_type> map_;
    size_t numCols_;
    Teuchos::Array<Teuchos::Array<global_ordinal_type> > sourceIndices_;
    Teuchos::Array<Teuchos::Array<scalar_type> > sourceValues_;

    size_t getNumColumns() const { return numCols_; }
  };

  /// \class MultiVectorFillerData2
  /// \brief Second implementation of fill and local assembly for \c MultiVectorFiller.
  /// \author Mark Hoemmen
  ///
  /// \tparam MV Specialization of \c Tpetra::MultiVector.
  template<class MV>
  class MultiVectorFillerData2 {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;

    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    /// \brief Default constructor (sets number of columns to zero).
    ///
    /// \param map [in] Map over which to distribute the initial fill.
    ///   This need not be the same Map as that of the multivector
    ///   output of \c globalAssemble().
    ///
    /// Before using this object, you should call \c setNumColumns()
    /// to set the number of columns in the output multivector.
    /// Otherwise, the two-argument version of \c
    /// sumIntoGlobalValues() won't actually do anything.
    MultiVectorFillerData2 (const Teuchos::RCP<const map_type>& map) : 
      map_ (map),
      numCols_ (0)
    {}

    /// \brief Constructor.
    ///
    /// \param map [in] Map over which to distribute the initial fill.
    ///   This need not be the same Map as that of the multivector
    ///   output of \c globalAssemble().
    ///
    /// \param numColumns [in] The (expected) number of columns in the
    ///   output multivector.  You can always change this later by
    ///   calling \c setNumColumns().
    ///
    /// \note If the number of columns given here is not the same as
    ///   the number of columns in the output multivector, you should
    ///   call \c setNumColumns() first before inserting any data.
    ///   Otherwise, the two-argument version of \c
    ///   sumIntoGlobalValues() won't do the right thing.
    MultiVectorFillerData2 (const Teuchos::RCP<const map_type>& map,
			    const size_t numColumns) : 
      map_ (map),
      numCols_ (numColumns),
      localVec_ (new MV (map, numColumns)),
      nonlocalIndices_ (numColumns),
      nonlocalValues_ (numColumns)
    {}

    //! All source indices (local and nonlocal) of the source Map, sorted and unique.
    Teuchos::Array<global_ordinal_type> 
    getSourceIndices () const 
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Tpetra::global_size_t;
      typedef global_ordinal_type GO;
      typedef typename Array<GO>::iterator iter_type;

      // Get the nonlocal row indices, sorted and made unique.
      // It's fair to assume that these are not contiguous.
      Array<GO> nonlocalIndices = getSortedUniqueNonlocalIndices();

      // Get the local row indices, not necessarily sorted or unique.
      ArrayView<const GO> localIndices = getLocalIndices ();

      Array<GO> indices (localIndices.size() + nonlocalIndices.size());
      std::copy (localIndices.begin(), localIndices.end(), indices.begin());
      std::sort (indices.begin(), indices.end());

      // Merge the local and nonlocal indices.
      if (nonlocalIndices.size() > 0) {
	//iter_type middle = &indices[nonlocalIndices.size()];
	iter_type middle;
	{
	  ArrayView<GO> indView = indices.view (0, nonlocalIndices.size());
	  middle = indView.end();
	}
	std::copy (nonlocalIndices.begin(), nonlocalIndices.end(), middle);
	std::inplace_merge (indices.begin(), middle, indices.end());
      }
      return indices;
    }

    /// \brief Set the number of columns in the output multivector.
    ///
    /// Setting the number of columns to zero effectively clears out
    /// all local storage, but may not necessarily deallocate nonlocal
    /// storage.  Call \c clear() to clear out all nonlocal storage.
    void
    setNumColumns (const size_t newNumColumns) 
    {
      using Teuchos::Array;
      using Teuchos::Range1D;
      using Teuchos::RCP;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;

      const size_t oldNumColumns = numCols_;
      if (newNumColumns == oldNumColumns) {
	return; // No side effects if no change.
      }

      RCP<MV> newLocalVec;
      if (newNumColumns > oldNumColumns) {
	newLocalVec = rcp (new MV (map_, newNumColumns));
	// Assign the contents of the old local multivector to the
	// first oldNumColumns columns of the new local multivector,
	// then get rid of the old local multivector.
	RCP<MV> newLocalVecView = 
	  newLocalVec->subViewNonConst (Range1D (0, oldNumColumns-1));
	*newLocalVecView = *localVec_;
      } 
      else {
	if (newNumColumns == 0) {
	  // Tpetra::MultiVector doesn't let you construct a
	  // multivector with zero columns.
	  newLocalVec = Teuchos::null;
	}
	else {
	  newLocalVec = 
	    localVec_->subViewNonConst (Range1D (0, newNumColumns-1));
	}
      }

      // Leave most side effects until the end, for exception safety.
      nonlocalIndices_.resize (newNumColumns);
      nonlocalValues_.resize (newNumColumns);
      localVec_ = newLocalVec;
      numCols_ = newNumColumns;
    }

    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 size_t columnIndex,
			 Teuchos::ArrayView<const scalar_type> values)
    {
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;

      if (columnIndex >= getNumColumns()) {
	// Automatically expand the number of columns.  This
	// implicitly ensures that localVec_ is not null.
	setNumColumns (columnIndex + 1);
      }
      
      typename ArrayView<const GO>::const_iterator rowIter = rows.begin();
      typename ArrayView<const ST>::const_iterator valIter = values.begin();
      for ( ; rowIter != rows.end() && valIter != values.end(); ++rowIter, ++valIter) {
	const GO globalRowIndex = *rowIter;
	// Converting from global to local index could be logarithmic
	// in the number of global indices that this process owns,
	// depending on the Map implementation.  However, the lookup
	// allows us to store data in the local multivector, rather
	// than in a separate data structure.
	const LO localRowIndex = map_->getLocalElement (globalRowIndex);
	if (localRowIndex == Teuchos::OrdinalTraits<LO>::invalid()) {
	  nonlocalIndices_[columnIndex].push_back (globalRowIndex);
	  nonlocalValues_[columnIndex].push_back (*valIter);
	}
	else {
	  // FIXME (mfh 27 Feb 2012) This will be very slow for GPU
	  // Node types.  In that case, we should hold on to the view
	  // of localVec_ as long as the number of columns doesn't
	  // change, and make modifications to the view until
	  // localAssemble() is called.
	  ArrayRCP<ST> X_j = localVec_->getDataNonConst (columnIndex);
	  // FIXME (mfh 27 Feb 2012) Allow different combine modes.
	  // The current combine mode just adds to the current value
	  // at that location.
	  X_j[localRowIndex] += *valIter;
	}
      }
    }

    /// Data for each column are stored contiguously in rows and in
    /// values.  Thus, rows and values are in rowwise order, even
    /// though they may be stored in columnwise order in the
    /// multivector.
    ///
    /// Be sure that the number of columns is set correctly before
    /// calling this.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 Teuchos::ArrayView<const scalar_type> values)
    {
      using Teuchos::ArrayView;
      typedef typename ArrayView<const global_ordinal_type>::size_type size_type;

      const size_t numCols = getNumColumns();
      for (size_t j = 0; j < numCols; ++j) {
	const size_type offset = numCols*j;
	const size_type len = numCols;
	sumIntoGlobalValues (rows.view (offset, len), j, values.view (offset, len));
      }
    }

    /// \brief Locally assemble into X, with user-specified combine mode.
    ///
    /// \param X [in/out] Multivector (overlapping source distribution).
    ///
    /// \param f [in/out] Binary function that defines the combine
    ///   mode.  It must define scalar_type operator (const
    ///   scalar_type&, const scalar_type&).  It need not necessarily
    ///   be commutative or even associative, but it should be
    ///   thread-safe in case we decide to parallelize local assembly.
    ///   We call it via X(i,j) = f(X(i,j), Y(i,j)), so write your
    ///   possibly nonassociative or noncommutative operation
    ///   accordingly.
    ///
    /// X is distributed by the source Map (with possible overlap) of
    /// the Export operation.  The source Map of the Export includes
    /// both the elements owned by this object's constructor's input
    /// Map, and the indices inserted by \c sumIntoGlobalValues().
    ///
    /// Precondition: The set of global indices in X's Map equals the
    /// union of the global indices in map_ and (the union of the
    /// entries of nonlocalIndices_[j] for all valid columns j).
    ///
    /// \note You can get the usual Tpetra::ADD combine mode by
    ///   supplying f = std::plus<scalar_type>.
    template<class BinaryFunction>
    void 
    locallyAssemble (MV& X, BinaryFunction& f)
    {
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;

      RCP<const map_type> srcMap = X.getMap();
      ArrayView<const GO> localIndices = map_->getNodeElementList ();

      for (size_t j = 0; j < X.getNumVectors(); ++j) {
	ArrayRCP<ST> X_j = X.getDataNonConst (j);

	// First add all the local data into X_j.
	ArrayRCP<const ST> local_j = localVec_->getDataNonConst (j);
	for (typename ArrayView<const GO>::const_iterator it = localIndices.begin(); 
	     it != localIndices.end(); ++it) {
	  const LO rowIndLocal = map_->getLocalElement (*it);
	  const LO rowIndX = srcMap->getLocalElement (*it);

	  TEUCHOS_TEST_FOR_EXCEPTION(rowIndX == Teuchos::OrdinalTraits<LO>::invalid(), 
            std::invalid_argument, "locallyAssemble(): Input multivector X does "
            "not own the global index " << *it << ".  This probably means that "
            "X was not constructed with the right Map.");
	  // FIXME (mfh 27 Feb 2012) We hard-code the ADD combine mode
	  // for now.  Later, accept other combine modes.
	  X_j[rowIndX] = f (X_j[rowIndX], local_j[rowIndLocal]);
	}

	// Now add the nonlocal data into X_j.
	ArrayView<const GO> nonlocalIndices = nonlocalIndices_[j]();
	typename ArrayView<const GO>::const_iterator indexIter = nonlocalIndices.begin(); 
	ArrayView<const ST> nonlocalValues = nonlocalValues_[j]();
	typename ArrayView<const ST>::const_iterator valueIter = nonlocalValues.begin();
	for ( ; indexIter != nonlocalIndices.end() && valueIter != nonlocalValues.end();
	      ++indexIter, ++valueIter) {
	  const LO rowIndX = srcMap->getLocalElement (*indexIter);
	  X_j[rowIndX] = f (X_j[rowIndX], *valueIter);
	}
      }
    }

    //! \c locallyAssemble() for the usual Tpetra::ADD combine mode.
    void 
    locallyAssemble (MV& X)
    {
      std::plus<double> f;
      locallyAssemble<std::plus<scalar_type> > (X, f);
    }

    /// \brief Clear the contents of the multivector.
    /// 
    /// This fills the vector with zeros, and also removes nonlocal
    /// data.  It does <i>not</i> deallocate all storage.  For that,
    /// you need to set the number of columns to zero.
    void clear() {
      Teuchos::Array<Teuchos::Array<global_ordinal_type> > newNonlocalIndices;
      Teuchos::Array<Teuchos::Array<scalar_type> > newNonlocalValues;
      // The standard STL idiom for clearing the contents of a vector
      // completely.  Setting the size to zero may not necessarily
      // deallocate data.
      std::swap (nonlocalIndices_, newNonlocalIndices);
      std::swap (nonlocalValues_, newNonlocalValues);

      // Don't actually deallocate the multivector of local entries.
      // Just fill it with zero.  This is because the caller hasn't
      // reset the number of columns.
      if (! localVec_.is_null()) {
	localVec_->putScalar (Teuchos::ScalarTraits<scalar_type>::zero());
      }
    }

  private:
    //! Map that tells us which entries are locally owned.
    Teuchos::RCP<const map_type> map_;

    //! The number of columns in the (output) multivector.
    size_t numCols_;

    //! Multivector of locally stored (i.e., owned by \c map_) entries.
    Teuchos::RCP<MV> localVec_;

    //! Nonlocal global indices: one array for each column of the multivector.
    Teuchos::Array<Teuchos::Array<global_ordinal_type> > nonlocalIndices_;

    //! Nonlocal global values: one array for each column of the multivector.
    Teuchos::Array<Teuchos::Array<scalar_type> > nonlocalValues_;

    //! The number of columns in the (output) multivector.
    size_t getNumColumns() const { return numCols_; }

    //! The locally owned row indices, not necessarily sorted or unique.
    Teuchos::ArrayView<const global_ordinal_type> 
    getLocalIndices() const
    {
      return map_->getNodeElementList ();
    }

    //! The inserted nonlocal row indices, sorted and made unique.
    Teuchos::Array<global_ordinal_type> 
    getSortedUniqueNonlocalIndices() const
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::as;
      typedef global_ordinal_type GO;
      typedef typename Array<GO>::iterator iter_type;
      typedef typename Array<GO>::size_type size_type;

      Array<GO> allNonlocals (0); // will resize below
      const size_t numCols = getNumColumns();

      if (numCols == 1) {
	// Special case for 1 column avoids copying indices.  Pick the
	// size of the array exactly the first time so there are at
	// most two allocations (the constructor may choose to
	// allocate).
	allNonlocals.resize (nonlocalIndices_[0].size());
	std::copy (nonlocalIndices_[0].begin(), nonlocalIndices_[0].end(),
		   allNonlocals.begin());
	std::sort (allNonlocals.begin(), allNonlocals.end());
	iter_type it = std::unique (allNonlocals.begin(), allNonlocals.end());
	const size_type curSize = as<size_type> (it - allNonlocals.begin());
	allNonlocals.resize (curSize);
      }
      else {
	// Carefully collect all the row indices one column at a time.
	// This ensures that the total allocation size in this routine
	// is independent of the number of columns.  Also, only sort
	// the current column's indices.  Use merges to ensure sorted
	// order in the collected final result.
	Array<GO> curNonlocals;
	for (size_t j = 0; j < numCols; ++j) {
	  // Collect the current column's nonlocal indices into
	  // curNonlocals.  Sort it and make the entries unique.
	  curNonlocals.resize (nonlocalIndices_[j].size());
	  std::copy (nonlocalIndices_[j].begin(), nonlocalIndices_[j].end(),
		     curNonlocals.begin());
	  std::sort (curNonlocals.begin(), curNonlocals.end());
	  iter_type it = std::unique (curNonlocals.begin(), curNonlocals.end());
	  const size_type curSize = as<size_type> (it - curNonlocals.begin());
	  curNonlocals.resize (curSize);

	  if (curSize > 0) {
	    // Assume inductively that the entries of allNonlocals are
	    // sorted and unique.  Make more space in allNonlocals for
	    // the new entries (since std::inplace_merge can only work
	    // in one array), copy them in, merge, and make the results
	    // unique.
	    const size_type oldSize = allNonlocals.size();
	    allNonlocals.resize (oldSize + curSize);

	    //iter_type middle = &allInds[oldSize];
	    iter_type middle;
	    {
	      ArrayView<GO> oldAllNonlocals = allNonlocals.view (oldSize, curSize);
	      middle = oldAllNonlocals.begin();
	    }
	    std::copy (curNonlocals.begin(), curNonlocals.end(), middle);
	    std::inplace_merge (allNonlocals.begin(), middle, allNonlocals.end());
	    iter_type it2 = std::unique (allNonlocals.begin(), allNonlocals.end());
	    const size_type newSize = as<size_type> (it2 - allNonlocals.begin());
	    allNonlocals.resize (newSize);
	  }
	}
      }
      return allNonlocals;
    }
  };
} // namespace (anonymous)

namespace Tpetra {

  /// \class MultiVectorFiller
  /// \brief Adds nonlocal sum-into functionality to Tpetra::MultiVector.
  /// \author Mark Hoemmen
  ///
  /// \tparam MV Specialization of Tpetra::MultiVector.
  template<class MV>
  class MultiVectorFiller {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;
    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    /// \brief Constructor.
    ///
    /// The constructor takes the same arguments as MV's constructor.
    /// We did this on purpose, so that you can construct
    /// MultiVectorFiller as you would a standard MV.
    ///
    /// \param map [in] A Map with the same communicator and Kokkos
    ///   Node as the output multivector of \c globalAssemble().
    ///
    /// \param numCols [in] Expected number of columns in the output
    ///   multivector of \c globalAssemble().
    ///
    /// \note If the number of columns given here is not the same as
    ///   the number of columns in the output multivector, the
    ///   two-argument version of \c sumIntoGlobalValues() won't do
    ///   the right thing.  Use the three-argument version of that
    ///   method if you don't know how many columns there will be in
    ///   advance.
    ///
    /// \note Not providing the output multivector in the constructor
    ///   gives \c globalAssemble() more flexibility.  For example,
    ///   its output multivector may have any distribution with the
    ///   same global number of rows.  Furthermore, decoupling
    ///   insertion from the output multivector lets this object store
    ///   preassembled data in whatever format it likes.  It doesn't
    ///   have to insert directly into the output multivector until
    ///   assembly time.  (This may improve performance by amortizing
    ///   the cost of global->local index conversions, for example.)
    MultiVectorFiller (const Teuchos::RCP<const map_type>& map, 
		       const size_t numCols);

    /// \brief Assemble all the data (local and nonlocal) into X_out.
    ///
    /// You can call this method multiple times with different
    /// multivector output arguments.  If those arguments have the
    /// same Map, this method will attempt to reuse the Export object
    /// each time.  It will only reuse if the new target Map is the
    /// same as (in the sense of \c isSameAs()) the previous target
    /// Map, unless you force reuse with the second argument (that
    /// saves a few global reductions for the check).
    ///
    /// \param X_out [in/out] MultiVector with a nonoverlapping Map.
    ///   That Map need not be the same as the one with which this
    ///   object was created.
    ///
    /// \param forceReuseMap [in] Assume that X_out has the same Map
    ///   (in the sense of \c isSameAs()) as the target Map in the
    ///   previous call to this method.  If this method was not called
    ///   before, then assume that X_out has the same Map as the
    ///   argument to the constructor of this object.
    void globalAssemble (MV& X_out, const bool forceReuseMap = false);

    /// \brief Sum data into the multivector.
    ///
    /// \param rows [in] Array of global rows for which to insert
    ///   values.  Must have the same length as values.
    ///
    /// \param column [in] Index of the column in which to insert.
    ///
    /// \param values [in] Array of values to insert.  Must have the
    ///   same length as rows.  rows[i] is the row in which values[i]
    ///   is to be inserted.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 size_t column,
			 Teuchos::ArrayView<const scalar_type> values)
    {
      data_.sumIntoGlobalValues (rows, column, values);
    }

    /// \brief Sum data into the multivector.
    ///
    /// In rows and values, the data for each column are stored
    /// contiguously.  Thus, each array is really a matrix in rowwise
    /// order.  (The multivector may use columnwise order internally.)
    ///
    /// Be sure that the number of columns is set correctly before
    /// calling this.
    ///
    /// \param rows [in] Array of global rows for which to insert
    ///   values.  Must have the same length as values.
    ///
    /// \param values [in] Array of values to insert.  Must have the
    ///   same length as rows.  rows[i] is the row in which values[i]
    ///   is to be inserted.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 Teuchos::ArrayView<const scalar_type> values)
    {
      data_.sumIntoGlobalValues (rows, values);
    }

  private:
    //! Map with which this object was created ("ctor" == "constructor").
    Teuchos::RCP<const map_type> ctorMap_;

    /// \brief Source Map of the last call to \c globalAssemble().
    ///
    /// A possibly overlapping Map that describes the distribution of
    /// input data, and is the source of the Export.
    Teuchos::RCP<const map_type> sourceMap_;

    /// \brief The target Map of the last call to \c globalAssemble().
    ///
    /// A nonoverlapping Map which is the target of the Export, and
    /// describes the distribution of the output multivector of \c
    /// globalAssemble().
    Teuchos::RCP<const map_type> targetMap_;

    /// \brief The source MV of the Export operation in \c globalAssemble().
    /// 
    /// The \c globalAssemble() method uses this multivector as the
    /// source of the Export operation.  The Export redistributes the
    /// data from a possibly overlapping distribution (reflecting how
    /// elements were inserted) to a nonoverlapping distribution (that
    /// of the output multivector \c X_out of \c globalAssemble()).
    /// This is the simplest way to implement redistribution
    ///
    /// We avoid resizing and reallocating \c sourceVec_ by using a
    /// contiguous subview as the source of the Export, if \c X_out
    /// has fewer columns than \c sourceVec_.
    Teuchos::RCP<MV> sourceVec_;

    /// \brief Storage for inserted indices and values.
    ///
    /// We separate this to facilitate experimentation with different
    /// storage formats.
    MultiVectorFillerData2<MV> data_;

    typedef Tpetra::Export<local_ordinal_type, global_ordinal_type, node_type> export_type;
    Teuchos::RCP<export_type> exporter_;

    /// \brief Assemble the local data into \c X_in.
    ///
    /// This method is called by \c globalAssemble(), in which \c X_in
    /// is the multivector with the (possibly overlapping) source
    /// distribution.
    void locallyAssemble (MV& X_in) {
      data_.locallyAssemble (X_in);
    }

    //! All source indices (local and nonlocal) of the source Map.
    Teuchos::Array<global_ordinal_type> getSourceIndices () const {
      return data_.getSourceIndices();
    }

    /// \brief Compute the source Map (for the Export) from the source indices.
    ///
    /// indexBase, comm, and node are the same as the input arguments
    /// of the Map constructors.
    ///
    /// \note This is a collective operation.
    ///
    /// \return The (possibly overlapping) source Map.
    Teuchos::RCP<const map_type>
    computeSourceMap (const global_ordinal_type indexBase,
		      const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
		      const Teuchos::RCP<node_type>& node);
  };

  template<class MV>
  MultiVectorFiller<MV>::MultiVectorFiller (const Teuchos::RCP<const typename MultiVectorFiller<MV>::map_type>& map, const size_t numCols) 
    : ctorMap_ (map), 
      sourceMap_ (Teuchos::null), 
      targetMap_ (Teuchos::null), 
      data_ (map, numCols), 
      exporter_ (Teuchos::null)
  {}

  template<class MV>
  Teuchos::RCP<const typename MultiVectorFiller<MV>::map_type>
  MultiVectorFiller<MV>::
  computeSourceMap (const global_ordinal_type indexBase,
		    const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
		    const Teuchos::RCP<node_type>& node)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::rcp;
    typedef global_ordinal_type GO;

    Array<GO> indices = getSourceIndices ();

    // Passing "invalid" for the numGlobalElements argument of the Map
    // constructor tells the Map to compute the global number of
    // elements itself.
    const Tpetra::global_size_t invalid = 
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    return rcp (new map_type (invalid, indices, indexBase, comm, node));
  }

  template<class MV>
  void 
  MultiVectorFiller<MV>::globalAssemble (MV& X_out, const bool forceReuseMap)
  {
    using Teuchos::ArrayView;
    using Teuchos::Array;
    using Teuchos::Range1D;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef global_ordinal_type GO;

    const size_t numVecs = X_out.getNumVectors();

    if (numVecs == 0) {
      // Nothing to do!  Of course, this does not check for whether
      // X_out has the right number of rows.  That's OK, though.
      // Compare to the definition of the BLAS' _AXPY for an input
      // vector containing NaNs, but multiplied by alpha=0.
      return;
    }
    //
    // Get the target Map of the Export.  If X_out's Map is the same
    // as the target Map from last time, then we may be able to
    // recycle the Export from last time, if the source Map is also
    // the same.
    //
    RCP<const map_type> targetMap;
    bool assumeSameTargetMap = false;
    if (targetMap_.is_null()) {
      targetMap_ = X_out.getMap();
      targetMap = targetMap_;
      assumeSameTargetMap = false;
    }
    else {
      if (forceReuseMap) {
	targetMap = targetMap_;
	assumeSameTargetMap = true;
      }
      else  {
	// If X_out's Map is the same as targetMap_, we may be able to
	// reuse the Export.  Constructing the Export may be more
	// expensive than calling isSameAs() (which requires just a
	// few reductions and reading through the lists of owned
	// global indices), so it's worth checking.
	if (targetMap_->isSameAs (*(X_out.getMap()))) {
	  assumeSameTargetMap = true;
	  targetMap = targetMap_;
	}
      }
    }
    //
    // Get the source Map of the Export.  If the source Map of the
    // Export is the same as last time, then we may be able to recycle
    // the Export from last time, if the target Map is also the same.
    //
    RCP<const map_type> sourceMap;
    bool computedSourceMap = false;
    {
      if (forceReuseMap && ! sourceMap_.is_null()) {
	sourceMap = sourceMap_;
      }
      else {
	sourceMap = computeSourceMap (ctorMap_->getIndexBase(), 
				      ctorMap_->getComm(), 
				      ctorMap_->getNode());
        computedSourceMap = true;
      }
    }
    if (computedSourceMap) {
      sourceMap_ = sourceMap;
    }
    // 
    // Now that we have the source and target Maps of the Export, we
    // can check whether we can recycle the Export from last time.
    //
    const bool mustComputeExport = 
      (exporter_.is_null() || (assumeSameTargetMap && ! computedSourceMap));
    if (mustComputeExport) {
      exporter_ = rcp (new export_type (sourceMap_, targetMap_));
    }

    // Source multivector for the Export.
    RCP<MV> X_in;
    const bool mustReallocateVec = sourceVec_.is_null() || 
      sourceVec_->getNumVectors() < numVecs || ! assumeSameTargetMap;

    if (mustReallocateVec) {
      // Reallocating nonlocalVec_ ensures that it has the right Map.
      sourceVec_ = rcp (new MV (sourceMap_, numVecs));
      X_in = sourceVec_;
    } else {
      if (sourceVec_->getNumVectors() == numVecs) {
	X_in = sourceVec_;
      } else { // sourceVec_ has more vectors than needed.
	X_in = sourceVec_->subViewNonConst (Range1D (0, numVecs-1));
      }
    }

    // "Locally assemble" the data into X_in by summing together
    // entries with the same indices.
    locallyAssemble (*X_in);
    
    // Do the Export.
    const Tpetra::CombineMode combineMode = Tpetra::ADD;
    X_in->doExport (X_out, *exporter_, combineMode);
  }

  namespace Test {

    /// \class MultiVectorFillerTester
    /// \brief Tests for \c MultiVectorFiller
    /// \author Mark Hoemmen
    ///
    /// \tparam MV A specialization of \c Tpetra::MultiVector.
    template<class MV>
    class MultiVectorFillerTester {
    public:
      typedef typename MV::scalar_type scalar_type;
      typedef typename MV::local_ordinal_type local_ordinal_type;
      typedef typename MV::global_ordinal_type global_ordinal_type;
      typedef typename MV::node_type node_type;
      typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

      /// \brief Test global assembly when constructor Map = target Map.
      ///
      /// Constructor Map = target Map is a common case for finite
      /// element assembly.  This method current only tests the version
      /// of \c sumIntoGlobalValues() that works on one column at a
      /// time.
      ///
      /// If any test fails, this method throws an exception.
      static void 
      testSameMap (const Teuchos::RCP<const map_type>& targetMap, 
		   const global_ordinal_type eltSize, // Must be odd
		   const size_t numCols)
      {
	using Teuchos::Array;
	using Teuchos::ArrayRCP;
	using Teuchos::ArrayView;
	using Teuchos::as;
	using Teuchos::Comm;
	using Teuchos::ptr;
	using Teuchos::RCP;
	using Teuchos::rcp;
	using Teuchos::REDUCE_SUM;
	using Teuchos::reduceAll;
	using std::cerr;
	using std::endl;

	typedef local_ordinal_type LO;
	typedef global_ordinal_type GO;
	typedef scalar_type ST;
	typedef Teuchos::ScalarTraits<ST> STS;

	TEUCHOS_TEST_FOR_EXCEPTION(eltSize % 2 == 0, std::invalid_argument,
				   "Element size (eltSize) argument must be odd.");
	TEUCHOS_TEST_FOR_EXCEPTION(numCols == 0, std::invalid_argument,
				   "Number of columns (numCols) argument must be nonzero.");

	//RCP<MV> X = rcp (new MV (targetMap, numCols));

	Array<GO> rows (eltSize);
	Array<ST> values (eltSize);
	std::fill (values.begin(), values.end(), STS::one());

	// Make this a pointer so we can free its contents, in case
	// those contents depend on the input to globalAssemble().
	RCP<MultiVectorFiller<MV> > filler = 
	  rcp (new MultiVectorFiller<MV> (targetMap, numCols));

	TEUCHOS_TEST_FOR_EXCEPTION(! targetMap->isContiguous(), 
          std::invalid_argument, "MultiVectorFiller test currently only works "
          "for contiguous Maps.");

	const GO minGlobalIndex = targetMap->getMinGlobalIndex();
	const GO maxGlobalIndex = targetMap->getMaxGlobalIndex();
	const GO minAllGlobalIndex = targetMap->getMinAllGlobalIndex();
	const GO maxAllGlobalIndex = targetMap->getMaxAllGlobalIndex();
	for (size_t j = 0; j < numCols; ++j) {
	  for (GO i = minGlobalIndex; i <= maxGlobalIndex; ++i) {
	    // Overlap over processes, without running out of bounds.
	    const GO start = std::max (i - eltSize/2, minAllGlobalIndex);
	    const GO end = std::min (i + eltSize/2, maxAllGlobalIndex);
	    const GO len = end - start + 1;

	    TEUCHOS_TEST_FOR_EXCEPTION(len > eltSize, std::logic_error,
              "At start,end = " << start << "," << end << ", len = " << len 
              << " > eltSize = " << eltSize << ".");

	    for (GO k = 0; k < len; ++k) {
	      rows[k] = start + k;
	    }
	    filler->sumIntoGlobalValues (rows.view(0, len), j, values.view(0, len));
	  }
	}

	MV X_out (targetMap, numCols);
	filler->globalAssemble (X_out);
	filler = Teuchos::null;

	const GO indexBase = targetMap->getIndexBase();
	Array<GO> errorLocations;
	for (size_t j = 0; j < numCols; ++j) {
	  ArrayRCP<const ST> X_j = X_out.getData (j);

	  // Each entry of the column should have the value eltSize,
	  // except for the first and last few entries in the whole
	  // column (globally, not locally).
	  for (GO i = minGlobalIndex; i <= maxGlobalIndex; ++i) {
	    const LO localIndex = targetMap->getLocalElement (i);
	    TEUCHOS_TEST_FOR_EXCEPTION(i == Teuchos::OrdinalTraits<LO>::invalid(),
	      std::logic_error, "Global index " << i << " is not in the "
	      "multivector's Map.");

	    if (i <= minAllGlobalIndex + eltSize/2) {
	      if (X_j[localIndex] != STS::one() + as<ST>(i) - as<ST>(indexBase)) {
		errorLocations.push_back (i);
	      }
	    } 
	    else if (i >= maxAllGlobalIndex - eltSize/2) {
	      if (X_j[localIndex] != STS::one() + as<ST>(maxAllGlobalIndex) - as<ST>(i)) {
		errorLocations.push_back (i);
	      }
	    }
	    else {
	      if (X_j[localIndex] != as<ST>(eltSize)) {
		errorLocations.push_back (i);
	      }
	    }
	  } // for each global index which my process owns

	  const typename Array<GO>::size_type localNumErrors = errorLocations.size();
	  typename Array<GO>::size_type globalNumErrors = 0;
	  RCP<const Comm<int> > comm = targetMap->getComm();
	  reduceAll (*comm, REDUCE_SUM, localNumErrors, ptr (&globalNumErrors));

	  if (globalNumErrors != 0) {
	    std::ostringstream os;
	    os << "Proc " << comm->getRank() << ": " << localNumErrors 
	       << " incorrect value" << (localNumErrors != 1 ? "s" : "") 
	       << ".  Error locations: [ ";
	    std::copy (errorLocations.begin(), errorLocations.end(),
		       std::ostream_iterator<GO> (os, " "));
	    os << "].";
	    // Iterate through all processes in the communicator,
	    // printing out each process' local errors.
	    for (int p = 0; p < comm->getSize(); ++p) {
	      if (p == comm->getRank()) {
		cerr << os.str() << endl;
	      }
	      // Barriers to let output finish.
	      comm->barrier();
	      comm->barrier();
	      comm->barrier();
	    }
	    TEUCHOS_TEST_FOR_EXCEPTION(globalNumErrors != 0, std::logic_error, 
	      "Over all procs: " << globalNumErrors << " total error" 
	      << (globalNumErrors != 1 ? "s" : "") << ".");
	  } // if there were any errors in column j
	} // for each column j
      }
    };

    //! Instantiate a \c MultiVectorFillerTester and run the test.
    template<class ScalarType, 
	     class LocalOrdinalType, 
	     class GlobalOrdinalType, 
	     class NodeType>
    void 
    testMultiVectorFiller (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
			   const Teuchos::RCP<NodeType>& node,
			   const size_t unknownsPerNode,
			   const GlobalOrdinalType unknownsPerElt,
			   const size_t numCols)
    {
      using Tpetra::createContigMapWithNode;
      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using std::cerr;
      using std::endl;

      typedef ScalarType ST;
      typedef LocalOrdinalType LO;
      typedef GlobalOrdinalType GO;
      typedef NodeType NT;
      typedef Tpetra::Map<LO, GO, NT> MT;
      typedef Tpetra::MultiVector<ST, LO, GO, NT> MV;

      const Tpetra::global_size_t invalid = 
	Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      RCP<const MT> targetMap = 
	createContigMapWithNode<LO, GO, NT> (invalid, unknownsPerNode, comm, node);

      std::ostringstream os;
      int success = 1;
      try {
	MultiVectorFillerTester<MV>::testSameMap (targetMap, unknownsPerElt, numCols);
	success = 1;
      } catch (std::exception& e) {
	success = 0;
	os << e.what();
      }

      for (int p = 0; p < comm->getSize(); ++p) {
	if (p == comm->getRank()) {
	  cerr << "On process " << comm->getRank() << ": " << os.str() << endl;
	}
	comm->barrier();
	comm->barrier();
	comm->barrier();
      }
    }
  } // namespace Test
} // namespace Tpetra


#endif // __Tpetra_MultiVectorFiller_hpp
