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

#ifndef IFPACK2_REORDERFILTER_DEF_HPP
#define IFPACK2_REORDERFILTER_DEF_HPP
#include "Ifpack2_ReorderFilter_decl.hpp"
#include <vector>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

namespace Ifpack2 {

template<class MatrixType>
ReorderFilter<MatrixType>::
ReorderFilter (const Teuchos::RCP<const row_matrix_type>& A,
               const Teuchos::ArrayRCP<local_ordinal_type>& perm,
               const Teuchos::ArrayRCP<local_ordinal_type>& reverseperm)
  : A_ (A),
    perm_ (perm),
    reverseperm_ (reverseperm)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::invalid_argument,
    "Ifpack2::ReorderFilter: The input matrix is null.");

  // use this filter only on serial matrices
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_->getComm()->getSize() != 1, std::invalid_argument,
    "Ifpack2::ReorderFilter: This class may only be used if the input matrix's "
    "communicator has one process.  This class is an implementation detail of "
    "Ifpack2::AdditiveSchwarz, and it is not meant to be used otherwise.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_->getNodeNumRows () != A_->getGlobalNumRows (),
    std::invalid_argument,
    "Ifpack2::ReorderFilter: The input matrix is not square.");

  // Temp arrays for apply
  Indices_.resize (A_->getNodeMaxNumRowEntries ());
  Values_.resize (A_->getNodeMaxNumRowEntries ());
}


template<class MatrixType>
ReorderFilter<MatrixType>::~ReorderFilter() {}


template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> > ReorderFilter<MatrixType>::getComm() const
{
  return A_->getComm();
}


template<class MatrixType>
Teuchos::RCP<typename ReorderFilter<MatrixType>::node_type>
ReorderFilter<MatrixType>::getNode () const
{
  return A_->getNode ();
}


template<class MatrixType>
Teuchos::RCP<const typename ReorderFilter<MatrixType>::map_type>
ReorderFilter<MatrixType>::getRowMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ReorderFilter::"
    "getRowMap: The matrix A is null, so there is no row Map.");

  return A_->getRowMap ();
}


template<class MatrixType>
Teuchos::RCP<const typename ReorderFilter<MatrixType>::map_type>
ReorderFilter<MatrixType>::getColMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ReorderFilter::"
    "getColMap: The matrix A is null, so there is no column Map.");

  return A_->getColMap();
}


template<class MatrixType>
Teuchos::RCP<const typename ReorderFilter<MatrixType>::map_type>
ReorderFilter<MatrixType>::getDomainMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ReorderFilter::"
    "getDomainMap: The matrix A is null, so there is no domain Map.");

  return A_->getDomainMap();
}


template<class MatrixType>
Teuchos::RCP<const typename ReorderFilter<MatrixType>::map_type>
ReorderFilter<MatrixType>::getRangeMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ReorderFilter::"
    "getRangeMap: The matrix A is null, so there is no range Map.");

  return A_->getRangeMap();
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type> >
ReorderFilter<MatrixType>::getGraph() const
{
  throw std::runtime_error("Ifpack2::ReorderFilter: does not support getGraph.");
}


template<class MatrixType>
global_size_t ReorderFilter<MatrixType>::getGlobalNumRows() const
{
  return A_->getGlobalNumRows();
}


template<class MatrixType>
global_size_t ReorderFilter<MatrixType>::getGlobalNumCols() const
{
  return A_->getGlobalNumCols();
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::getNodeNumRows() const
{
  return A_->getNodeNumRows();
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::getNodeNumCols() const
{
  return A_->getNodeNumCols();
}


template<class MatrixType>
typename MatrixType::global_ordinal_type ReorderFilter<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}


template<class MatrixType>
global_size_t ReorderFilter<MatrixType>::getGlobalNumEntries() const
{
  return A_->getGlobalNumEntries();
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::getNodeNumEntries() const
{
  return A_->getNodeNumEntries();
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::
getNumEntriesInGlobalRow (global_ordinal_type globalRow) const
{
  if (A_.is_null () || A_->getRowMap ().is_null ()) {
    return Teuchos::OrdinalTraits<size_t>::invalid ();
  }
  else {
    const local_ordinal_type lclRow =
      A_->getRowMap ()->getLocalElement (globalRow);
    if (lclRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid ()) {
      // The calling process doesn't own any entries in this row.
      return static_cast<size_t> (0);
    } else {
      const local_ordinal_type origLclRow = reverseperm_[lclRow];
      return A_->getNumEntriesInLocalRow (origLclRow);
    }
  }
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::
getNumEntriesInLocalRow (local_ordinal_type localRow) const
{
  // Make sure that localRow is in bounds before using it to index
  // into the permutation.
  if (A_->getRowMap ()->isNodeLocalElement (localRow)) {
    // localRow is a valid index into reverseperm_.
    const local_ordinal_type localReorderedRow = reverseperm_[localRow];
    return A_->getNumEntriesInLocalRow (localReorderedRow);
  } else {
    // The calling process doesn't own any entries in this row.
    return static_cast<size_t> (0);
  }
}


template<class MatrixType>
global_size_t ReorderFilter<MatrixType>::getGlobalNumDiags() const
{
  return A_->getGlobalNumDiags();
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::getNodeNumDiags() const
{
  return A_->getNodeNumDiags();
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return A_->getGlobalMaxNumRowEntries();
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::getNodeMaxNumRowEntries() const
{
  return A_->getNodeMaxNumRowEntries();
}


template<class MatrixType>
bool ReorderFilter<MatrixType>::hasColMap() const
{
  return true;
}


template<class MatrixType>
bool ReorderFilter<MatrixType>::isLowerTriangular() const
{
  return A_->isLowerTriangular();
}


template<class MatrixType>
bool ReorderFilter<MatrixType>::isUpperTriangular() const
{
  return A_->isUpperTriangular();
}


template<class MatrixType>
bool ReorderFilter<MatrixType>::isLocallyIndexed() const
{
  return A_->isLocallyIndexed();
}


template<class MatrixType>
bool ReorderFilter<MatrixType>::isGloballyIndexed() const
{
  return A_->isGloballyIndexed();
}


template<class MatrixType>
bool ReorderFilter<MatrixType>::isFillComplete() const
{
  return A_->isFillComplete();
}


template<class MatrixType>
void ReorderFilter<MatrixType>::
getGlobalRowCopy (global_ordinal_type globalRow,
                  const Teuchos::ArrayView<global_ordinal_type>& globalInd,
                  const Teuchos::ArrayView<scalar_type>& val,
                  size_t& numEntries) const
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::av_reinterpret_cast;
  typedef local_ordinal_type LO;
  typedef global_ordinal_type GO;
  typedef Teuchos::OrdinalTraits<LO> OTLO;

  const map_type& rowMap = * (A_->getRowMap ());
  const local_ordinal_type localRow = rowMap.getLocalElement (globalRow);
  TEUCHOS_TEST_FOR_EXCEPTION(
    localRow == OTLO::invalid (), std::invalid_argument, "Ifpack2::Reorder"
    "Filter::getGlobalRowCopy: The given global row index " << globalRow
    << " is not owned by the calling process with rank "
    << rowMap.getComm ()->getRank () << ".");

  if (sizeof (GO) == sizeof (LO)) {
    // This means we can convert local to global in place.
    ArrayView<LO> localInd = av_reinterpret_cast<LO> (globalInd);
    this->getLocalRowCopy (localRow, localInd, val, numEntries);

    // Convert local indices back to global indices.
    for (size_t k = 0; k < numEntries; ++k) {
      globalInd[k] = rowMap.getGlobalElement (localInd[k]);
    }
  }
  else {
    // LO and GO have different sizes, so we need a temp array
    // for converting local to global.
    numEntries = this->getNumEntriesInLocalRow (localRow);
    Array<LO> localInd (numEntries);
    this->getLocalRowCopy (localRow, localInd, val, numEntries);

    // Convert local indices back to global indices.
    for (size_t k = 0; k < numEntries; ++k) {
      globalInd[k] = rowMap.getGlobalElement (localInd[k]);
    }
  }
}


template<class MatrixType>
void ReorderFilter<MatrixType>::
getLocalRowCopy (local_ordinal_type LocalRow,
                 const Teuchos::ArrayView<local_ordinal_type> &Indices,
                 const Teuchos::ArrayView<scalar_type> &Values,
                 size_t &NumEntries) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A_->getRowMap ()->isNodeLocalElement (LocalRow),
    std::invalid_argument,
    "Ifpack2::ReorderFilter::getLocalRowCopy: The given local row index "
    << LocalRow << " is not a valid local row index on the calling process "
    "with rank " << A_->getRowMap ()->getComm ()->getRank () << ".");

  // This duplicates code in getNumEntriesInGlobalRow, but avoids an
  // extra array lookup and some extra tests.
  const local_ordinal_type origLclRow = reverseperm_[LocalRow];
  const size_t numEntries = A_->getNumEntriesInLocalRow (origLclRow);

  TEUCHOS_TEST_FOR_EXCEPTION(
    static_cast<size_t> (Indices.size ()) < numEntries ||
    static_cast<size_t> (Values.size ()) < numEntries,
    std::invalid_argument,
    "Ifpack2::ReorderFilter::getLocalRowCopy: The given array views are not "
    "long enough to store all the data in the given row " << LocalRow
    << ".  Indices.size() = " << Indices.size () << ", Values.size() = "
    << Values.size () << ", but the (original) row has " << numEntries
    << " entry/ies.");

  A_->getLocalRowCopy (origLclRow, Indices, Values, NumEntries);
  // Do a col reindex via perm
  //
  // FIXME (mfh 30 Jan 2014) This assumes that the row and column
  // indices are the same.
  for (size_t i = 0; i < NumEntries; ++i) {
    Indices[i] = perm_[Indices[i]];
  }
}


template<class MatrixType>
void ReorderFilter<MatrixType>::
getGlobalRowView (global_ordinal_type GlobalRow,
                  Teuchos::ArrayView<const global_ordinal_type> &indices,
                  Teuchos::ArrayView<const scalar_type> &values) const
{
  throw std::runtime_error("Ifpack2::ReorderFilter: does not support getGlobalRowView.");
}


template<class MatrixType>
void ReorderFilter<MatrixType>::
getLocalRowView (local_ordinal_type LocalRow,
                 Teuchos::ArrayView<const local_ordinal_type> &indices,
                 Teuchos::ArrayView<const scalar_type> &values) const
{
  throw std::runtime_error("Ifpack2::ReorderFilter: does not support getLocalRowView.");
}


template<class MatrixType>
void ReorderFilter<MatrixType>::
getLocalDiagCopy (Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &diag) const
{
  // This is somewhat dubious as to how the maps match.
  return A_->getLocalDiagCopy(diag);
}


template<class MatrixType>
void ReorderFilter<MatrixType>::leftScale(const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x)
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not support leftScale.");
}


template<class MatrixType>
void ReorderFilter<MatrixType>::rightScale(const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x)
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not support rightScale.");
}


template<class MatrixType>
void ReorderFilter<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  TEUCHOS_TEST_FOR_EXCEPTION(
    alpha != STS::one () || beta != STS::zero (), std::logic_error,
    "Ifpack2::ReorderFilter::apply is only implemented for alpha = 1 and "
    "beta = 0.  You set alpha = " << alpha << " and beta = " << beta << ".");

  // Note: This isn't AztecOO compliant.  But neither was Ifpack's version.
  // Note: The localized maps mean the matvec is trivial (and has no import)
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Ifpack2::ReorderFilter::apply: X.getNumVectors() != Y.getNumVectors().");

  const scalar_type zero = STS::zero ();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const scalar_type> > x_ptr = X.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > y_ptr = Y.get2dViewNonConst();

  Y.putScalar (zero);
  const size_t NumVectors = Y.getNumVectors ();

  for (size_t i = 0; i < A_->getNodeNumRows (); ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    getLocalRowCopy (i, Indices_ (), Values_ (), Nnz);
    if (mode == Teuchos::NO_TRANS) {
      for (size_t j = 0; j < Nnz; ++j) {
        for (size_t k = 0; k < NumVectors; ++k) {
          y_ptr[k][i] += Values_[j] * x_ptr[k][Indices_[j]];
        }
      }
    }
    else if (mode == Teuchos::TRANS) {
      for (size_t j = 0; j < Nnz; ++j) {
        for (size_t k = 0; k < NumVectors; ++k) {
          y_ptr[k][Indices_[j]] += Values_[j] * x_ptr[k][i];
        }
      }
    }
    else { //mode==Teuchos::CONJ_TRANS
      for (size_t j = 0; j < Nnz; ++j) {
        for (size_t k = 0; k < NumVectors; ++k) {
          y_ptr[k][Indices_[j]] += STS::conjugate(Values_[j]) * x_ptr[k][i];
        }
      }
    }
  }
}


template<class MatrixType>
bool ReorderFilter<MatrixType>::hasTransposeApply() const
{
  return true;
}


template<class MatrixType>
bool ReorderFilter<MatrixType>::supportsRowViews() const
{
  return false;
}


template<class MatrixType>
typename ReorderFilter<MatrixType>::mag_type ReorderFilter<MatrixType>::getFrobeniusNorm() const
{
  // Reordering doesn't change the Frobenius norm.
  return A_->getFrobeniusNorm ();
}


template<class MatrixType>
void ReorderFilter<MatrixType>::
permuteOriginalToReordered (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &originalX,
                            Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &reorderedY) const
{
  this->template permuteOriginalToReorderedTempl<scalar_type,scalar_type>(originalX, reorderedY);
}


template<class MatrixType>
template<class DomainScalar, class RangeScalar>
void ReorderFilter<MatrixType>::permuteOriginalToReorderedTempl(const Tpetra::MultiVector<DomainScalar,local_ordinal_type,global_ordinal_type,node_type> &originalX,
                                                                Tpetra::MultiVector<RangeScalar,local_ordinal_type,global_ordinal_type,node_type> &reorderedY) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(originalX.getNumVectors() != reorderedY.getNumVectors(), std::runtime_error,
                             "Ifpack2::ReorderFilter::permuteOriginalToReordered ERROR: X.getNumVectors() != Y.getNumVectors().");

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const DomainScalar> > x_ptr = originalX.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<RangeScalar> >        y_ptr = reorderedY.get2dViewNonConst();

  for(size_t k=0; k < originalX.getNumVectors(); k++)
    for(local_ordinal_type i=0; (size_t)i< originalX.getLocalLength(); i++)
      y_ptr[k][perm_[i]] = (RangeScalar)x_ptr[k][i];
}


template<class MatrixType>
void ReorderFilter<MatrixType>::permuteReorderedToOriginal(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &reorderedX,
                                                           Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &originalY) const
{
  this->template permuteReorderedToOriginalTempl<scalar_type,scalar_type>(reorderedX, originalY);
}


template<class MatrixType>
template<class DomainScalar, class RangeScalar>
void ReorderFilter<MatrixType>::
permuteReorderedToOriginalTempl (const Tpetra::MultiVector<DomainScalar,local_ordinal_type,global_ordinal_type,node_type> &reorderedX,
                                 Tpetra::MultiVector<RangeScalar,local_ordinal_type,global_ordinal_type,node_type> &originalY) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    reorderedX.getNumVectors() != originalY.getNumVectors(),
    std::runtime_error,
    "Ifpack2::ReorderFilter::permuteReorderedToOriginal: "
    "X.getNumVectors() != Y.getNumVectors().");

#ifdef HAVE_IFPACK2_DEBUG
  {
    typedef Teuchos::ScalarTraits<DomainScalar> STS;
    typedef typename STS::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    Teuchos::Array<magnitude_type> norms (reorderedX.getNumVectors ());
    reorderedX.norm2 (norms ());
    bool good = true;
    for (typename Teuchos::Array<magnitude_type>::size_type j = 0;
         j < reorderedX.getNumVectors (); ++j) {
      if (STM::isnaninf (norms[j])) {
        good = false;
        break;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! good, std::runtime_error, "Ifpack2::ReorderFilter::"
      "permuteReorderedToOriginalTempl: The 2-norm of the input reorderedX is "
      "NaN or Inf.");
  }
#endif // HAVE_IFPACK2_DEBUG

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const DomainScalar> > x_ptr = reorderedX.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<RangeScalar> >        y_ptr = originalY.get2dViewNonConst();

  for (size_t k = 0; k < reorderedX.getNumVectors (); ++k) {
    for (local_ordinal_type i = 0; (size_t)i < reorderedX.getLocalLength (); ++i) {
      y_ptr[k][reverseperm_[i]] = (RangeScalar) x_ptr[k][i];
    }
  }

#ifdef HAVE_IFPACK2_DEBUG
  {
    typedef Teuchos::ScalarTraits<RangeScalar> STS;
    typedef typename STS::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    Teuchos::Array<magnitude_type> norms (originalY.getNumVectors ());
    originalY.norm2 (norms ());
    bool good = true;
    for (typename Teuchos::Array<magnitude_type>::size_type j = 0;
         j < originalY.getNumVectors (); ++j) {
      if (STM::isnaninf (norms[j])) {
        good = false;
        break;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! good, std::runtime_error, "Ifpack2::ReorderFilter::"
      "permuteReorderedToOriginalTempl: The 2-norm of the output originalY is "
      "NaN or Inf.");
  }
#endif // HAVE_IFPACK2_DEBUG
}

} // namespace Ifpack2

#define IFPACK2_REORDERFILTER_INSTANT(S,LO,GO,N)                        \
  template class Ifpack2::ReorderFilter< Tpetra::CrsMatrix<S, LO, GO, N> >;

#endif
