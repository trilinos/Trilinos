/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_LOCALFILTER_DEF_HPP
#define IFPACK2_LOCALFILTER_DEF_HPP

#include <Ifpack2_LocalFilter_decl.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>

#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif

namespace Ifpack2 {


template<class MatrixType>
bool
LocalFilter<MatrixType>::
mapPairsAreFitted (const row_matrix_type& A)
{
  const map_type& rangeMap = * (A.getRangeMap ());
  const map_type& rowMap = * (A.getRowMap ());
  const bool rangeAndRowFitted = mapPairIsFitted (rangeMap, rowMap);

  const map_type& domainMap = * (A.getDomainMap ());
  const map_type& columnMap = * (A.getColMap ());
  const bool domainAndColumnFitted = mapPairIsFitted (domainMap, columnMap);

  return rangeAndRowFitted && domainAndColumnFitted;
}


template<class MatrixType>
bool
LocalFilter<MatrixType>::
mapPairIsFitted (const map_type& map1, const map_type& map2)
{
  using Teuchos::ArrayView;
  typedef global_ordinal_type GO; // a handy abbreviation
  typedef typename ArrayView<const GO>::size_type size_type;

  bool fitted = true;
  if (&map1 == &map2) {
    fitted = true;
  }
  else if (map1.isContiguous () && map2.isContiguous () &&
           map1.getMinGlobalIndex () == map2.getMinGlobalIndex () &&
           map1.getMaxGlobalIndex () <= map2.getMaxGlobalIndex ()) {
    // Special case where both Maps are contiguous.
    fitted = true;
  }
  else {
    ArrayView<const GO> inds_map2 = map2.getNodeElementList ();
    const size_type numInds_map1 = static_cast<size_type> (map1.getNodeNumElements ());

    if (map1.isContiguous ()) {
      // Avoid calling getNodeElementList() on the always one-to-one
      // Map, if it is contiguous (a common case).  When called on a
      // contiguous Map, getNodeElementList() causes allocation of an
      // array that sticks around, even though the array isn't needed.
      // (The Map is contiguous, so you can compute the entries; you
      // don't need to store them.)
      if (numInds_map1 > inds_map2.size ()) {
        // There are fewer indices in map1 on this process than in
        // map2.  This case might be impossible.
        fitted = false;
      }
      else {
        // Do all the map1 indices match the initial map2 indices?
        const GO minInd_map1 = map1.getMinGlobalIndex ();
        for (size_type k = 0; k < numInds_map1; ++k) {
          const GO inds_map1_k = static_cast<GO> (k) + minInd_map1;
          if (inds_map1_k != inds_map2[k]) {
            fitted = false;
            break;
          }
        }
      }
    }
    else { // map1 is not contiguous.
      // Get index lists from both Maps, and compare their indices.
      ArrayView<const GO> inds_map1 = map1.getNodeElementList ();
      if (numInds_map1 > inds_map2.size ()) {
        // There are fewer indices in map1 on this process than in
        // map2.  This case might be impossible.
        fitted = false;
      }
      else {
        // Do all the map1 indices match those in map2?
        for (size_type k = 0; k < numInds_map1; ++k) {
          if (inds_map1[k] != inds_map2[k]) {
            fitted = false;
            break;
          }
        }
      }
    }
  }
  return fitted;
}


template<class MatrixType>
LocalFilter<MatrixType>::
LocalFilter (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
  NumNonzeros_ (0),
  MaxNumEntries_ (0),
  MaxNumEntriesA_ (0)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

#ifdef HAVE_IFPACK2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! mapPairsAreFitted (*A), std::invalid_argument, "Ifpack2::LocalFilter: "
    "A's Map pairs are not fitted to each other on Process "
    << A_->getRowMap ()->getComm ()->getRank () << " of the input matrix's "
    "communicator.  "
    "This means that LocalFilter does not currently know how to work with A.  "
    "This will change soon.  Please see discussion of Bug 5992.");
#endif // HAVE_IFPACK2_DEBUG

  // Build the local communicator (containing this process only).
  RCP<const Teuchos::Comm<int> > localComm;
#ifdef HAVE_MPI
  localComm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));
#else
  localComm = rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_MPI


  // // FIXME (mfh 21 Nov 2013) Currently, the implementation implicitly
  // // assumes that the range Map is fitted to the row Map.  Otherwise,
  // // it probably won't work at all.
  // TEUCHOS_TEST_FOR_EXCEPTION(
  //   mapPairIsFitted (* (A_->getRangeMap ()), * (A_->getRowMap ())),
  //   std::logic_error, "Ifpack2::LocalFilter: Range Map of the input matrix "
  //   "is not fitted to its row Map.  LocalFilter does not currently work in "
  //   "this case.  See Bug 5992.");

  // Build the local row, domain, and range Maps.  They both use the
  // local communicator built above.  The global indices of each are
  // different than those of the corresponding original Map; they
  // actually are the same as the local indices of the original Map.
  //
  // It's even OK if signedness of local_ordinal_type and
  // global_ordinal_type are different.  (That would be a BAD IDEA but
  // some users seem to enjoy making trouble for themselves and us.)
  //
  // Both the local row and local range Maps must have the same number
  // of entries, namely, that of the local number of entries of A's
  // range Map.

  // FIXME (mfh 21 Nov 2013) For some reason, we have to use this,
  // even if it differs from the number of entries in the range Map.
  // Otherwise, AdditiveSchwarz Test1 fails, down in the local solve,
  // where the matrix has 8 columns but the local part of the vector
  // only has five rows.
  const size_t numRows = A_->getNodeNumRows ();

  // using std::cerr;
  // using std::endl;
  // cerr << "A_ has " << A_->getNodeNumRows () << " rows." << endl
  //      << "Range Map has " << A_->getRangeMap ()->getNodeNumElements () << " entries." << endl
  //      << "Row Map has " << A_->getRowMap ()->getNodeNumElements () << " entries." << endl;

  const global_ordinal_type indexBase = static_cast<global_ordinal_type> (0);

  localRowMap_ =
    rcp (new map_type (numRows, indexBase, localComm,
                       Tpetra::GloballyDistributed, A_->getNode ()));
  // If the original matrix's range Map is not fitted to its row Map,
  // we'll have to do an Export when applying the matrix.
  localRangeMap_ = localRowMap_;

  // If the original matrix's domain Map is not fitted to its column
  // Map, we'll have to do an Import when applying the matrix.
  const size_t numCols = A_->getDomainMap ()->getNodeNumElements ();
  if (A_->getRangeMap ().getRawPtr () == A_->getDomainMap ().getRawPtr ()) {
    // The input matrix's domain and range Maps are identical, so the
    // locally filtered matrix's domain and range Maps can be also.
    localDomainMap_ = localRangeMap_;
  }
  else {
    localDomainMap_ =
      rcp (new map_type (numCols, indexBase, localComm,
                         Tpetra::GloballyDistributed, A_->getNode ()));
  }

  // NodeNumEntries_ will contain the actual number of nonzeros for
  // each localized row (that is, without external nodes, and always
  // with the diagonal entry)
  NumEntries_.resize (numRows);

  // tentative value for MaxNumEntries. This is the number of
  // nonzeros in the local matrix
  MaxNumEntries_  = A_->getNodeMaxNumRowEntries ();
  MaxNumEntriesA_ = A_->getNodeMaxNumRowEntries ();

  // Allocate temporary arrays for getLocalRowCopy().
  localIndices_.resize (MaxNumEntries_);
  Values_.resize (MaxNumEntries_);

  // now compute:
  // - the number of nonzero per row
  // - the total number of nonzeros
  // - the diagonal entries

  // compute nonzeros (total and per-row), and store the
  // diagonal entries (already modified)
  size_t ActualMaxNumEntries = 0;

  for (size_t i = 0; i < numRows; ++i) {
    NumEntries_[i] = 0;
    size_t Nnz, NewNnz = 0;
    A_->getLocalRowCopy (i, localIndices_, Values_, Nnz);
    for (size_t j = 0; j < Nnz; ++j) {
      // FIXME (mfh 03 Apr 2013) This assumes the following:
      //
      // 1. Row Map, range Map, and domain Map are all the same.
      //
      // 2. The column Map's list of GIDs on this process is the
      //    domain Map's list of GIDs, followed by remote GIDs.  Thus,
      //    for any GID in the domain Map on this process, its LID in
      //    the domain Map (and therefore in the row Map, by (1)) is
      //    the same as its LID in the column Map.  (Hence the
      //    less-than test, which if true, means that localIndices_[j]
      //    belongs to the row Map.)
      if (static_cast<size_t> (localIndices_[j]) < numRows) {
        ++NewNnz;
      }
    }

    if (NewNnz > ActualMaxNumEntries) {
      ActualMaxNumEntries = NewNnz;
    }

    NumNonzeros_ += NewNnz;
    NumEntries_[i] = NewNnz;
  }

  MaxNumEntries_ = ActualMaxNumEntries;
}


template<class MatrixType>
LocalFilter<MatrixType>::~LocalFilter()
{}


template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
LocalFilter<MatrixType>::getComm () const
{
  return localRowMap_->getComm ();
}


template<class MatrixType>
Teuchos::RCP<typename MatrixType::node_type>
LocalFilter<MatrixType>::getNode () const
{
  return A_->getNode ();
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
LocalFilter<MatrixType>::getRowMap () const
{
  return localRowMap_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
LocalFilter<MatrixType>::getColMap() const
{
  return localRowMap_; // FIXME (mfh 20 Nov 2013)
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
LocalFilter<MatrixType>::getDomainMap() const
{
  return localDomainMap_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
LocalFilter<MatrixType>::getRangeMap() const
{
  return localRangeMap_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
LocalFilter<MatrixType>::getGraph () const
{
  // FIXME (mfh 20 Nov 2013) This is not what the documentation says
  // this method should do!  It should return the graph of the locally
  // filtered matrix, not the original matrix's graph.
  return A_->getGraph ();
}


template<class MatrixType>
global_size_t LocalFilter<MatrixType>::getGlobalNumRows() const
{
  return static_cast<global_size_t> (localRangeMap_->getNodeNumElements ());
}


template<class MatrixType>
global_size_t LocalFilter<MatrixType>::getGlobalNumCols() const
{
  return static_cast<global_size_t> (localDomainMap_->getNodeNumElements ());
}


template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeNumRows() const
{
  return static_cast<size_t> (localRangeMap_->getNodeNumElements ());
}


template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeNumCols() const
{
  return static_cast<size_t> (localDomainMap_->getNodeNumElements ());
}


template<class MatrixType>
typename MatrixType::global_ordinal_type
LocalFilter<MatrixType>::getIndexBase () const
{
  return A_->getIndexBase ();
}


template<class MatrixType>
global_size_t LocalFilter<MatrixType>::getGlobalNumEntries () const
{
  return NumNonzeros_;
}


template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeNumEntries () const
{
  return NumNonzeros_;
}


template<class MatrixType>
size_t
LocalFilter<MatrixType>::
getNumEntriesInGlobalRow (global_ordinal_type globalRow) const
{
  const local_ordinal_type localRow = getRowMap ()->getLocalElement (globalRow);
  if (localRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid ()) {
    // NOTE (mfh 26 Mar 2014) We return zero if globalRow is not in
    // the row Map on this process, since "get the number of entries
    // in the global row" refers only to what the calling process owns
    // in that row.  In this case, it owns no entries in that row,
    // since it doesn't own the row.
    return 0;
  } else {
    return NumEntries_[localRow];
  }
}


template<class MatrixType>
size_t
LocalFilter<MatrixType>::
getNumEntriesInLocalRow (local_ordinal_type localRow) const
{
  // FIXME (mfh 07 Jul 2014) Shouldn't localRow be a local row index
  // in the matrix's row Map, not in the LocalFilter's row Map?  The
  // latter is different; it even has different global indices!
  // (Maybe _that_'s the bug.)

  if (getRowMap ()->isNodeLocalElement (localRow)) {
    return NumEntries_[localRow];
  } else {
    // NOTE (mfh 26 Mar 2014) We return zero if localRow is not in the
    // row Map on this process, since "get the number of entries in
    // the local row" refers only to what the calling process owns in
    // that row.  In this case, it owns no entries in that row, since
    // it doesn't own the row.
    return 0;
  }
}


template<class MatrixType>
global_size_t LocalFilter<MatrixType>::getGlobalNumDiags () const
{
  return A_->getGlobalNumDiags ();
}


template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeNumDiags () const
{
  return A_->getNodeNumDiags ();
}


template<class MatrixType>
size_t LocalFilter<MatrixType>::getGlobalMaxNumRowEntries () const
{
  return MaxNumEntries_;
}


template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeMaxNumRowEntries() const
{
  return MaxNumEntries_;
}


template<class MatrixType>
bool LocalFilter<MatrixType>::hasColMap () const
{
  return true;
}


template<class MatrixType>
bool LocalFilter<MatrixType>::isLowerTriangular () const
{
  return A_->isLowerTriangular();
}


template<class MatrixType>
bool LocalFilter<MatrixType>::isUpperTriangular () const
{
  return A_->isUpperTriangular();
}


template<class MatrixType>
bool LocalFilter<MatrixType>::isLocallyIndexed () const
{
  return A_->isLocallyIndexed ();
}


template<class MatrixType>
bool LocalFilter<MatrixType>::isGloballyIndexed () const
{
  return A_->isGloballyIndexed();
}


template<class MatrixType>
bool LocalFilter<MatrixType>::isFillComplete () const
{
  return A_->isFillComplete ();
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
getGlobalRowCopy (global_ordinal_type globalRow,
                  const Teuchos::ArrayView<global_ordinal_type>& globalIndices,
                  const Teuchos::ArrayView<scalar_type>& values,
                  size_t& numEntries) const
{
  typedef local_ordinal_type LO;
  typedef typename Teuchos::Array<LO>::size_type size_type;

  const LO localRow = this->getRowMap ()->getLocalElement (globalRow);
  if (localRow == Teuchos::OrdinalTraits<LO>::invalid ()) {
    // NOTE (mfh 26 Mar 2014) We return no entries if globalRow is not
    // in the row Map on this process, since "get a copy of the
    // entries in the global row" refers only to what the calling
    // process owns in that row.  In this case, it owns no entries in
    // that row, since it doesn't own the row.
    numEntries = 0;
  }
  else {
    // First get a copy of the current row using local indices.  Then,
    // convert to global indices using the input matrix's column Map.
    //
    numEntries = this->getNumEntriesInLocalRow (localRow);
    // FIXME (mfh 26 Mar 2014) If local_ordinal_type ==
    // global_ordinal_type, we could just alias the input array
    // instead of allocating a temporary array.
    Teuchos::Array<LO> localIndices (numEntries);
    this->getLocalRowCopy (localRow, localIndices (), values, numEntries);

    const map_type& colMap = * (this->getColMap ());

    // Don't fill the output array beyond its size.
    const size_type numEnt =
      std::min (static_cast<size_type> (numEntries),
                std::min (globalIndices.size (), values.size ()));
    for (size_type k = 0; k < numEnt; ++k) {
      globalIndices[k] = colMap.getGlobalElement (localIndices[k]);
    }
  }
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
getLocalRowCopy (local_ordinal_type LocalRow,
                 const Teuchos::ArrayView<local_ordinal_type> &Indices,
                 const Teuchos::ArrayView<scalar_type> &Values,
                 size_t &NumEntries) const
{
  typedef local_ordinal_type LO;
  typedef global_ordinal_type GO;

  if (! A_->getRowMap ()->isNodeLocalElement (LocalRow)) {
    // The calling process owns zero entries in the row.
    NumEntries = 0;
    return;
  }

  const size_t numEntInLclRow = NumEntries_[LocalRow];
  if (static_cast<size_t> (Indices.size ()) < numEntInLclRow ||
      static_cast<size_t> (Values.size ()) < numEntInLclRow) {
    // FIXME (mfh 07 Jul 2014) Return an error code instead of
    // throwing.  We should really attempt to fill as much space as
    // we're given, in this case.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error,
      "Ifpack2::LocalFilter::getLocalRowCopy: Invalid output array length.  "
      "The output arrays must each have length at least " << numEntInLclRow
      << " for local row " << LocalRow << " on Process "
      << localRowMap_->getComm ()->getRank () << ".");
  }
  else if (numEntInLclRow == static_cast<size_t> (0)) {
    // getNumEntriesInLocalRow() returns zero if LocalRow is not owned
    // by the calling process.  In that case, the calling process owns
    // zero entries in the row.
    NumEntries = 0;
    return;
  }

  // Always extract using the temporary arrays Values_ and
  // localIndices_.  This is because I may need more space than that
  // given by the user.  The users expects only the local (in the
  // domain Map) column indices, but I have to extract both local and
  // remote (not in the domain Map) column indices.
  //
  // FIXME (mfh 07 Jul 2014) Use of temporary local storage is not
  // conducive to thread parallelism.  A better way would be to change
  // the interface so that it only extracts values for the "local"
  // column indices.  CrsMatrix could take a set of column indices,
  // and return their corresponding values.
  size_t numEntInMat = 0;
  A_->getLocalRowCopy (LocalRow, localIndices_ (), Values_ (), numEntInMat);

  // Fill the user's arrays with the "local" indices and values in
  // that row.  Note that the matrix might have a different column Map
  // than the local filter.
  const map_type& matrixDomMap = * (A_->getDomainMap ());
  const map_type& matrixColMap = * (A_->getColMap ());

  const size_t capacity = static_cast<size_t> (std::min (Indices.size (),
                                                         Values.size ()));
  NumEntries = 0;
  const size_t numRows = localRowMap_->getNodeNumElements (); // superfluous
  const bool buggy = true; // mfh 07 Jul 2014: See FIXME below.
  for (size_t j = 0; j < numEntInMat; ++j) {
    // The LocalFilter only includes entries in the domain Map on
    // the calling process.  We figure out whether an entry is in
    // the domain Map by converting the (matrix column Map) index to
    // a global index, and then asking whether that global index is
    // in the domain Map.
    const LO matrixLclCol = localIndices_[j];
    const GO gblCol = matrixColMap.getGlobalElement (matrixLclCol);

    // FIXME (mfh 07 Jul 2014) This is the likely center of Bug 5992
    // and perhaps other bugs, like Bug 6117.  If 'buggy' is true,
    // Ifpack2 tests pass; if 'buggy' is false, the tests don't pass.
    // This suggests that Ifpack2 classes could be using LocalFilter
    // incorrectly, perhaps by giving it an incorrect domain Map.
    if (buggy) {
      // only local indices
      if ((size_t) localIndices_[j] < numRows) {
        Indices[NumEntries] = localIndices_[j];
        Values[NumEntries]  = Values_[j];
        NumEntries++;
      }
    } else {
      if (matrixDomMap.isNodeGlobalElement (gblCol)) {
        // Don't fill more space than the user gave us.  It's an error
        // for them not to give us enough space, but we still shouldn't
        // overwrite memory that doesn't belong to us.
        if (NumEntries < capacity) {
          Indices[NumEntries] = matrixLclCol;
          Values[NumEntries]  = Values_[j];
        }
        NumEntries++;
      }
    }
  }
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
getGlobalRowView (global_ordinal_type GlobalRow,
                  Teuchos::ArrayView<const global_ordinal_type> &indices,
                  Teuchos::ArrayView<const scalar_type> &values) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "Ifpack2::LocalFilter does not implement getGlobalRowView.");
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
getLocalRowView (local_ordinal_type LocalRow,
                 Teuchos::ArrayView<const local_ordinal_type> &indices,
                 Teuchos::ArrayView<const scalar_type> &values) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "Ifpack2::LocalFilter does not implement getLocalRowView.");
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
getLocalDiagCopy (Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& diag) const
{
  using Teuchos::RCP;
  typedef Tpetra::Vector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> vector_type;
  // This is always correct, and doesn't require a collective check
  // for sameness of Maps.
  RCP<vector_type> diagView = diag.offsetViewNonConst (A_->getRowMap (), 0);
  A_->getLocalDiagCopy (*diagView);
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
leftScale (const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "Ifpack2::LocalFilter does not implement leftScale.");
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
rightScale (const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "Ifpack2::LocalFilter does not implement rightScale.");
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Ifpack2::LocalFilter::apply: X and Y must have the same number of columns.  "
    "X has " << X.getNumVectors () << " columns, but Y has "
    << Y.getNumVectors () << " columns.");

  if (&X == &Y) {
    // FIXME (mfh 23 Apr 2014) We have to do more work to figure out
    // if X and Y alias one another.  For example, we should check
    // whether one is a noncontiguous view of the other.
    //
    // X and Y alias one another, so we have to copy X.
    MV X_copy (X, Teuchos::Copy);
    applyNonAliased (X_copy, Y, mode, alpha, beta);
  } else {
    applyNonAliased (X, Y, mode, alpha, beta);
  }
}

template<class MatrixType>
void
LocalFilter<MatrixType>::
applyNonAliased (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
                 Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
                 Teuchos::ETransp mode,
                 scalar_type alpha,
                 scalar_type beta) const
{
  using Teuchos::ArrayView;
  using Teuchos::ArrayRCP;
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  const scalar_type zero = STS::zero ();
  const scalar_type one = STS::one ();

  if (beta == zero) {
    Y.putScalar (zero);
  }
  else if (beta != one) {
    Y.scale (beta);
  }

  const size_t NumVectors = Y.getNumVectors ();
  const size_t numRows = localRowMap_->getNodeNumElements ();

  // FIXME (mfh 14 Apr 2014) At some point, we would like to
  // parallelize this using Kokkos.  This would require a
  // Kokkos-friendly version of getLocalRowCopy, perhaps with
  // thread-local storage.

  const bool constantStride = X.isConstantStride () && Y.isConstantStride ();
  if (constantStride) {
    // Since both X and Y have constant stride, we can work with "1-D"
    // views of their data.
    const size_t              x_stride = X.getStride();
    const size_t              y_stride = Y.getStride();
    ArrayRCP<scalar_type>        y_rcp = Y.get1dViewNonConst();
    ArrayRCP<const scalar_type>  x_rcp = X.get1dView();
    ArrayView<scalar_type>       y_ptr = y_rcp();
    ArrayView<const scalar_type> x_ptr = x_rcp();
    for (size_t i = 0; i < numRows; ++i) {
      size_t Nnz;
      // Use this class's getrow to make the below code simpler
      getLocalRowCopy (i, localIndices_ (), Values_ (), Nnz);
      if (mode == Teuchos::NO_TRANS) {
        for (size_t j = 0; j < Nnz; ++j) {
          const local_ordinal_type col = localIndices_[j];
          for (size_t k = 0; k < NumVectors; ++k) {
            y_ptr[i + y_stride*k] +=
              alpha * Values_[j] * x_ptr[col + x_stride*k];
          }
        }
      }
      else if (mode == Teuchos::TRANS) {
        for (size_t j = 0; j < Nnz; ++j) {
          const local_ordinal_type col = localIndices_[j];
          for (size_t k = 0; k < NumVectors; ++k) {
            y_ptr[col + y_stride*k] +=
              alpha * Values_[j] * x_ptr[i + x_stride*k];
          }
        }
      }
      else { //mode==Teuchos::CONJ_TRANS
        for (size_t j = 0; j < Nnz; ++j) {
          const local_ordinal_type col = localIndices_[j];
          for (size_t k = 0; k < NumVectors; ++k) {
            y_ptr[col + y_stride*k] +=
              alpha * STS::conjugate (Values_[j]) * x_ptr[i + x_stride*k];
          }
        }
      }
    }
  }
  else {
    // At least one of X or Y does not have constant stride.
    // This means we must work with "2-D" views of their data.
    ArrayRCP<ArrayRCP<const scalar_type> > x_ptr = X.get2dView();
    ArrayRCP<ArrayRCP<scalar_type> >       y_ptr = Y.get2dViewNonConst();

    for (size_t i = 0; i < numRows; ++i) {
      size_t Nnz;
      // Use this class's getrow to make the below code simpler
      getLocalRowCopy (i, localIndices_ (), Values_ (), Nnz);
      if (mode == Teuchos::NO_TRANS) {
        for (size_t k = 0; k < NumVectors; ++k) {
          ArrayView<const scalar_type> x_local = (x_ptr())[k]();
          ArrayView<scalar_type>       y_local = (y_ptr())[k]();
          for (size_t j = 0; j < Nnz; ++j) {
            y_local[i] += alpha * Values_[j] * x_local[localIndices_[j]];
          }
        }
      }
      else if (mode == Teuchos::TRANS) {
        for (size_t k = 0; k < NumVectors; ++k) {
          ArrayView<const scalar_type> x_local = (x_ptr())[k]();
          ArrayView<scalar_type>       y_local = (y_ptr())[k]();
          for (size_t j = 0; j < Nnz; ++j) {
            y_local[localIndices_[j]] += alpha * Values_[j] * x_local[i];
          }
        }
      }
      else { //mode==Teuchos::CONJ_TRANS
        for (size_t k = 0; k < NumVectors; ++k) {
          ArrayView<const scalar_type> x_local = (x_ptr())[k]();
          ArrayView<scalar_type>       y_local = (y_ptr())[k]();
          for (size_t j = 0; j < Nnz; ++j) {
            y_local[localIndices_[j]] +=
              alpha * STS::conjugate (Values_[j]) * x_local[i];
          }
        }
      }
    }
  }
}

template<class MatrixType>
bool LocalFilter<MatrixType>::hasTransposeApply () const
{
  return true;
}


template<class MatrixType>
bool LocalFilter<MatrixType>::supportsRowViews () const
{
  return false;
}


template<class MatrixType>
typename
LocalFilter<MatrixType>::mag_type
LocalFilter<MatrixType>::getFrobeniusNorm () const
{
#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
  typedef Kokkos::Details::ArithTraits<scalar_type> STS;
  typedef Kokkos::Details::ArithTraits<mag_type> STM;
#else
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
#endif
  typedef typename Teuchos::Array<scalar_type>::size_type size_type;

  const size_type maxNumRowEnt = getNodeMaxNumRowEntries ();
  Teuchos::Array<local_ordinal_type> ind (maxNumRowEnt);
  Teuchos::Array<scalar_type> val (maxNumRowEnt);
  const size_t numRows = static_cast<size_t> (localRowMap_->getNodeNumElements ());

  // FIXME (mfh 03 Apr 2013) Scale during sum to avoid overflow.
  mag_type sumSquared = STM::zero ();
  for (size_t i = 0; i < numRows; ++i) {
    size_t numEntries = 0;
    this->getLocalRowCopy (i, ind (), val (), numEntries);
    for (size_type k = 0; k < static_cast<size_type> (numEntries); ++k) {
      const mag_type v_k_abs = STS::magnitude (val[k]);
      sumSquared += v_k_abs * v_k_abs;
    }
  }
  return STM::squareroot (sumSquared); // Different for each process; that's OK.
}


template<class MatrixType>
TPETRA_DEPRECATED void
LocalFilter<MatrixType>::
getGlobalRowView (global_ordinal_type GlobalRow,
                  Teuchos::ArrayRCP<const global_ordinal_type> &indices,
                  Teuchos::ArrayRCP<const scalar_type> &values) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "Ifpack2::LocalFilter does not implement getGlobalRowView.");
}


template<class MatrixType>
TPETRA_DEPRECATED
void
LocalFilter<MatrixType>::
getLocalRowView (local_ordinal_type LocalRow,
                 Teuchos::ArrayRCP<const local_ordinal_type> &indices,
                 Teuchos::ArrayRCP<const scalar_type> &values) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "Ifpack2::LocalFilter does not implement getLocalRowView.");
}


template<class MatrixType>
std::string
LocalFilter<MatrixType>::description () const
{
  using Teuchos::TypeNameTraits;
  std::ostringstream os;

  os << "Ifpack2::LocalFilter: {";
  os << "MatrixType: " << TypeNameTraits<MatrixType>::name ();
  if (this->getObjectLabel () != "") {
    os << ", Label: \"" << this->getObjectLabel () << "\"";
  }
  os << ", Number of rows: " << getGlobalNumRows ()
     << ", Number of columns: " << getGlobalNumCols ()
     << "}";
  return os.str ();
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::OSTab;
  using Teuchos::TypeNameTraits;
  using std::endl;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

  if (vl > Teuchos::VERB_NONE) {
    // describe() starts with a tab, by convention.
    OSTab tab0 (out);

    out << "Ifpack2::LocalFilter:" << endl;
    OSTab tab1 (out);
    out << "MatrixType: " << TypeNameTraits<MatrixType>::name () << endl;
    if (this->getObjectLabel () != "") {
      out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
    }
    out << "Number of rows: " << getGlobalNumRows () << endl
        << "Number of columns: " << getGlobalNumCols () << endl
        << "Number of nonzeros: " << NumNonzeros_ << endl;

    if (vl > Teuchos::VERB_LOW) {
      out << "Row Map:" << endl;
      localRowMap_->describe (out, vl);
      out << "Domain Map:" << endl;
      localDomainMap_->describe (out, vl);
      out << "Range Map:" << endl;
      localRangeMap_->describe (out, vl);
    }
  }
}


} // namespace Ifpack2

#define IFPACK2_LOCALFILTER_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::LocalFilter< Tpetra::CrsMatrix<S, LO, GO, N> >;

#endif
