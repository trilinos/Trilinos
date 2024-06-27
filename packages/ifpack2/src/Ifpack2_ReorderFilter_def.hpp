// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
    A_->getLocalNumRows () != A_->getGlobalNumRows (),
    std::invalid_argument,
    "Ifpack2::ReorderFilter: The input matrix is not square.");

  // Temp arrays for apply
  Kokkos::resize(Indices_,A_->getLocalMaxNumRowEntries ());
  Kokkos::resize(Values_,A_->getLocalMaxNumRowEntries ());
}


template<class MatrixType>
ReorderFilter<MatrixType>::~ReorderFilter() {}


template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> > ReorderFilter<MatrixType>::getComm() const
{
  return A_->getComm();
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
size_t ReorderFilter<MatrixType>::getLocalNumRows() const
{
  return A_->getLocalNumRows();
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::getLocalNumCols() const
{
  return A_->getLocalNumCols();
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
size_t ReorderFilter<MatrixType>::getLocalNumEntries() const
{
  return A_->getLocalNumEntries();
}

template<class MatrixType>
typename MatrixType::local_ordinal_type ReorderFilter<MatrixType>::getBlockSize() const
{
  return A_->getBlockSize();
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
size_t ReorderFilter<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return A_->getGlobalMaxNumRowEntries();
}


template<class MatrixType>
size_t ReorderFilter<MatrixType>::getLocalMaxNumRowEntries() const
{
  return A_->getLocalMaxNumRowEntries();
}


template<class MatrixType>
bool ReorderFilter<MatrixType>::hasColMap() const
{
  return true;
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
                   nonconst_global_inds_host_view_type &globalInd,
                   nonconst_values_host_view_type &val,
                   size_t& numEntries) const
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::av_reinterpret_cast;
  typedef local_ordinal_type LO;
  typedef Teuchos::OrdinalTraits<LO> OTLO;

  const map_type& rowMap = * (A_->getRowMap ());
  const local_ordinal_type localRow = rowMap.getLocalElement (globalRow);
  TEUCHOS_TEST_FOR_EXCEPTION(
    localRow == OTLO::invalid (), std::invalid_argument, "Ifpack2::Reorder"
    "Filter::getGlobalRowCopy: The given global row index " << globalRow
    << " is not owned by the calling process with rank "
    << rowMap.getComm ()->getRank () << ".");

  // The Indices_ temp array is only used in apply, not getLocalRowCopy, so this is safe
  numEntries = this->getNumEntriesInLocalRow (localRow);
  this->getLocalRowCopy (localRow, Indices_, val, numEntries);

  // Convert local indices back to global indices.
  for (size_t k = 0; k < numEntries; ++k) {
    globalInd[k] = rowMap.getGlobalElement (Indices_[k]);
  }
}


template<class MatrixType>
void ReorderFilter<MatrixType>::
getLocalRowCopy (local_ordinal_type LocalRow,
    nonconst_local_inds_host_view_type &Indices,
    nonconst_values_host_view_type &Values,
    size_t& NumEntries) const

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
void ReorderFilter<MatrixType>::getGlobalRowView(global_ordinal_type /* GlobalRow */,
                                                  global_inds_host_view_type &/*indices*/,
                                                  values_host_view_type &/*values*/) const
{
  throw std::runtime_error("Ifpack2::ReorderFilter: does not support getGlobalRowView.");
}



template<class MatrixType>
void ReorderFilter<MatrixType>::getLocalRowView(local_ordinal_type /* LocalRow */,
    local_inds_host_view_type & /*indices*/,
    values_host_view_type & /*values*/) const
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
void ReorderFilter<MatrixType>::leftScale(const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& /* x */)
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not support leftScale.");
}


template<class MatrixType>
void ReorderFilter<MatrixType>::rightScale(const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& /* x */)
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

  for (size_t i = 0; i < A_->getLocalNumRows (); ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    getLocalRowCopy (i, Indices_ , Values_ , Nnz);
    scalar_type* Values = reinterpret_cast<scalar_type*>(Values_.data());
    if (mode == Teuchos::NO_TRANS) {
      for (size_t j = 0; j < Nnz; ++j) {
        for (size_t k = 0; k < NumVectors; ++k) {
          y_ptr[k][i] += Values[j] * x_ptr[k][Indices_[j]];
        }
      }
    }
    else if (mode == Teuchos::TRANS) {
      for (size_t j = 0; j < Nnz; ++j) {
        for (size_t k = 0; k < NumVectors; ++k) {
          y_ptr[k][Indices_[j]] += Values[j] * x_ptr[k][i];
        }
      }
    }
    else { //mode==Teuchos::CONJ_TRANS
      for (size_t j = 0; j < Nnz; ++j) {
        for (size_t k = 0; k < NumVectors; ++k) {
          y_ptr[k][Indices_[j]] += STS::conjugate(Values[j]) * x_ptr[k][i];
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

  const local_ordinal_type blockSize = getBlockSize();
  const local_ordinal_type numRows = originalX.getLocalLength() / blockSize;
  for(size_t k=0; k < originalX.getNumVectors(); k++)
    for(local_ordinal_type i=0; i< numRows; i++)
      for(local_ordinal_type j=0; j< blockSize; ++j)
        y_ptr[k][perm_[i]*blockSize + j] = (RangeScalar)x_ptr[k][i*blockSize + j];
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
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    Teuchos::Array<magnitude_type> norms (reorderedX.getNumVectors ());
    reorderedX.norm2 (norms ());
    bool good = true;
    for (size_t j = 0;
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

  const local_ordinal_type blockSize = getBlockSize();
  const local_ordinal_type numRows = reorderedX.getLocalLength() / blockSize;
  for (size_t k = 0; k < reorderedX.getNumVectors (); ++k) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      for(local_ordinal_type j = 0; j < blockSize; ++j) {
        y_ptr[k][reverseperm_[i]*blockSize + j] = (RangeScalar) x_ptr[k][i*blockSize + j];
      }
    }
  }

#ifdef HAVE_IFPACK2_DEBUG
  {
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    Teuchos::Array<magnitude_type> norms (originalY.getNumVectors ());
    originalY.norm2 (norms ());
    bool good = true;
    for (size_t j = 0;
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
  template class Ifpack2::ReorderFilter< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif
