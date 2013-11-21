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
LocalFilter<MatrixType>::
LocalFilter (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
  NumRows_ (0),
  NumNonzeros_ (0),
  MaxNumEntries_ (0),
  MaxNumEntriesA_ (0)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Communicator containing this process only.
  RCP<const Teuchos::Comm<int> > localComm;
#ifdef HAVE_MPI
  localComm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));
#else
  localComm = rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_MPI

  // localized matrix has all the local rows of Matrix
  NumRows_ = A_->getNodeNumRows ();

  // build a linear map, based on the serial communicator
  localMap_ = rcp (new map_type (NumRows_, 0, localComm));

  // NodeNumEntries_ will contain the actual number of nonzeros
  // for each localized row (that is, without external nodes,
  // and always with the diagonal entry)
  NumEntries_.resize (NumRows_);

  // tentative value for MaxNumEntries. This is the number of
  // nonzeros in the local matrix
  MaxNumEntries_  = A_->getNodeMaxNumRowEntries ();
  MaxNumEntriesA_ = A_->getNodeMaxNumRowEntries ();

  // Allocate temporary arrays for getLocalRowCopy().
  Indices_.resize (MaxNumEntries_);
  Values_.resize (MaxNumEntries_);

  // now compute:
  // - the number of nonzero per row
  // - the total number of nonzeros
  // - the diagonal entries

  // compute nonzeros (total and per-row), and store the
  // diagonal entries (already modified)
  size_t ActualMaxNumEntries = 0;

  for (size_t i = 0; i < NumRows_; ++i) {
    NumEntries_[i] = 0;
    size_t Nnz, NewNnz = 0;
    A_->getLocalRowCopy (i, Indices_, Values_, Nnz);
    for (size_t j = 0 ; j < Nnz ; ++j) {
      // FIXME (mfh 03 Apr 2013) This assumes the following:
      //
      // 1. Row Map, range Map, and domain Map are all the same.
      //
      // 2. The column Map's list of GIDs on this process is the
      //    domain Map's list of GIDs, followed by remote GIDs.  Thus,
      //    for any GID in the domain Map on this process, its LID in
      //    the domain Map (and therefore in the row Map, by (1)) is
      //    the same as its LID in the column Map.  (Hence the
      //    less-than test, which if true, means that Indices_[j]
      //    belongs to the row Map.)
      if (Teuchos::as<size_t> (Indices_[j]) < NumRows_) {
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
  return localMap_->getComm ();
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
  return localMap_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
LocalFilter<MatrixType>::getColMap() const
{
  return localMap_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
LocalFilter<MatrixType>::getDomainMap() const
{
  return localMap_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
LocalFilter<MatrixType>::getRangeMap() const
{
  return localMap_;
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
  return NumRows_;
}


template<class MatrixType>
global_size_t LocalFilter<MatrixType>::getGlobalNumCols() const
{
  return NumRows_;
}


template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeNumRows() const
{
  return NumRows_;
}


template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeNumCols() const
{
  return NumRows_;
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
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "Ifpack2::LocalFilter does not implement getNumEntriesInGlobalRow.");
}


template<class MatrixType>
size_t
LocalFilter<MatrixType>::
getNumEntriesInLocalRow (local_ordinal_type localRow) const
{
  return NumEntries_[localRow];
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
getGlobalRowCopy (global_ordinal_type GlobalRow,
                  const Teuchos::ArrayView<global_ordinal_type> &Indices,
                  const Teuchos::ArrayView<scalar_type> &Values,
                  size_t &NumEntries) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "Ifpack2::LocalFilter does not implement getGlobalRowCopy.");
}


template<class MatrixType>
void
LocalFilter<MatrixType>::
getLocalRowCopy (local_ordinal_type LocalRow,
                 const Teuchos::ArrayView<local_ordinal_type> &Indices,
                 const Teuchos::ArrayView<scalar_type> &Values,
                 size_t &NumEntries) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    NumRows_ == 0, std::invalid_argument,
    "Ifpack2::LocalFilter::getLocalRowCopy: Invalid local row index "
    << LocalRow << ".  This process " << localMap_->getComm ()->getSize ()
    << " owns no rows of the matrix.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    LocalRow < 0 || (size_t) LocalRow >= NumRows_,
    std::invalid_argument,
    "Ifpack2::LocalFilter::getLocalRowCopy: Invalid local row index "
    << LocalRow << ".  The valid range of row indices on this process "
    << localMap_->getComm ()->getSize () << " is "
    << "[0, " << (NumRows_ - 1) << "].");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (size_t) Indices.size() <  NumEntries_[LocalRow], std::runtime_error,
    "Ifpack2::LocalFilter::getLocalRowCopy: Invalid output array length.  "
    "The output arrays must each have length at least " << NumEntries_[LocalRow]
    << " for local row " << LocalRow << " on process "
    << localMap_->getComm ()->getSize () << ".");

  size_t A_NumEntries=0;
  // Always extract using the object Values_ and Indices_.  This is
  // because I may need more space than that given by the user.  The
  // users expects only the local (in the domain Map) column indices,
  // but I have to extract both local and remote (not in the domain
  // Map) column indices.
  A_->getLocalRowCopy (LocalRow, Indices_ (), Values_ (), A_NumEntries);

  // populate the user's vectors
  NumEntries = 0;
  for (size_t j = 0 ; j < A_NumEntries; ++j) {
    // only local indices
    if ((size_t) Indices_[j] < NumRows_) {
      Indices[NumEntries] = Indices_[j];
      Values[NumEntries]  = Values_[j];
      NumEntries++;
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
  using Teuchos::ArrayRCP;
  typedef Tpetra::Vector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> vector_type;

  vector_type temp (A_->getRowMap ());
  A_->getLocalDiagCopy (temp);

  // FIXME (mfh 12 July 2013) WHY DO WE NEED ANYTHING MORE AFTER THE
  // ABOVE???  AND WHY CAN'T WE USE Vector::operator= INSTEAD OF
  // COPYING ALL THE DATA BY HAND???

  ArrayRCP<ArrayRCP<scalar_type> >       d_ptr = diag.get2dViewNonConst();
  ArrayRCP<ArrayRCP<const scalar_type> > t_ptr = temp.get2dView();

  for (size_t i = 0; i < NumRows_; ++i) {
    d_ptr[0][i] = t_ptr[0][i];
  }
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
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Ifpack2::LocalFilter::apply: X and Y must have the same number of columns.  "
    "X has " << X.getNumVectors () << " columns, but Y has "
    << Y.getNumVectors () << " columns.");

  const scalar_type zero = STS::zero();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const scalar_type> > x_ptr = X.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> >       y_ptr = Y.get2dViewNonConst();

  if (beta == zero) {
    Y.putScalar (zero);
  }
  else {
    Y.scale (beta);
  }

  const size_t NumVectors = Y.getNumVectors();

  // FIXME (mfh 12 July 2013) This would be a good candidate for
  // parallelization via Kokkos.

  for (size_t i = 0 ; i < NumRows_ ; ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    getLocalRowCopy (i, Indices_ (), Values_ (), Nnz);
    if (mode == Teuchos::NO_TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j)
        for (size_t k = 0 ; k < NumVectors ; ++k)
          y_ptr[k][i] += alpha * Values_[j] * x_ptr[k][Indices_[j]];
    }
    else if (mode == Teuchos::TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j)
        for (size_t k = 0 ; k < NumVectors ; ++k)
          y_ptr[k][Indices_[j]] += alpha * Values_[j] * x_ptr[k][i];
    }
    else { //mode==Teuchos::CONJ_TRANS
      for (size_t j = 0 ; j < Nnz ; ++j)
        for (size_t k = 0 ; k < NumVectors ; ++k)
          y_ptr[k][Indices_[j]] += alpha * STS::conjugate(Values_[j]) * x_ptr[k][i];
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
Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
LocalFilter<MatrixType>::getFrobeniusNorm () const
{
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  Teuchos::ArrayView<local_ordinal_type> ind;
  Teuchos::ArrayView<scalar_type> val;
  const size_t numRows = NumRows_; // compilers like constant loop bounds

  // FIXME (mfh 03 Apr 2013) Scale during sum to avoid overflow.
  magnitude_type sumSquared = STM::zero ();
  if (STS::isComplex) {
    for (size_t i = 0; i < numRows; ++i) {
      size_t numEntries = 0;
      this->getLocalRowCopy (i, ind, val, numEntries);
      for (size_t k = 0; k < numEntries; ++k) {
        sumSquared += STS::real (val[k]) * STS::real (val[k]) +
          STS::imag (val[k]) * STS::imag (val[k]);
      }
    }
  }
  else {
    for (size_t i = 0; i < numRows; ++i) {
      size_t numEntries = 0;
      this->getLocalRowCopy (i, ind, val, numEntries);
      for (size_t k = 0; k < numEntries; ++k) {
        sumSquared += STS::magnitude(val[k]) * STS::magnitude(val[k]);
      }
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

} // namespace Ifpack2

#endif
