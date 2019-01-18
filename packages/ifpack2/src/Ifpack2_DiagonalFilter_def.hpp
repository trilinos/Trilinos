/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_DIAGONALFILTER_DEF_HPP
#define IFPACK2_DIAGONALFILTER_DEF_HPP
#include "Ifpack2_DiagonalFilter_decl.hpp"
#include <vector>

#include "Teuchos_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"


namespace Ifpack2 {

template<class MatrixType>
DiagonalFilter<MatrixType>::
DiagonalFilter (const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix,
                typename Teuchos::ScalarTraits<Scalar>::magnitudeType AbsoluteThreshold,
                typename Teuchos::ScalarTraits<Scalar>::magnitudeType RelativeThreshold):
  A_(Matrix),
  AbsoluteThreshold_(AbsoluteThreshold),
  RelativeThreshold_(RelativeThreshold)
{
  pos_.resize(getNodeNumRows());
  val_=Teuchos::rcp(new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A_->getRowMap()));

  std::vector<LocalOrdinal> Indices(getNodeMaxNumRowEntries());
  std::vector<Scalar> Values(getNodeMaxNumRowEntries());
  size_t NumEntries;
  magnitudeType mysign;


  for (size_t MyRow = 0 ; MyRow < getNodeNumRows() ; ++MyRow) {
    pos_[MyRow] = -1;
    A_->getLocalRowCopy(MyRow,Indices,Values,NumEntries);

    for (size_t i = 0 ; i < NumEntries ; ++i) {
      if ((size_t)Indices[i] == MyRow) {
        pos_[MyRow] = i;
        if (Teuchos::ScalarTraits<Scalar>::real(Values[i]) < 0) {
          mysign=-Teuchos::ScalarTraits<magnitudeType>::one();
        }
        else {
          mysign=Teuchos::ScalarTraits<magnitudeType>::one();
        }
        val_->replaceLocalValue (MyRow, Values[i] * (RelativeThreshold_ - 1) +
                                 AbsoluteThreshold_ * mysign);
        break;
      }
    }
  }
}

template<class MatrixType>
DiagonalFilter<MatrixType>::~DiagonalFilter() { }

template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> > DiagonalFilter<MatrixType>::getComm() const
{
  return A_->getComm();
}

template<class MatrixType>
Teuchos::RCP <typename MatrixType::node_type> DiagonalFilter<MatrixType>::getNode() const
{
  return A_->getNode();
}

template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
DiagonalFilter<MatrixType>::getRowMap() const
{
  return A_->getRowMap();
}

template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
DiagonalFilter<MatrixType>::getColMap() const
{
  return A_->getColMap();
}

template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
DiagonalFilter<MatrixType>::getDomainMap() const
{
  return A_->getDomainMap();
}

template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
DiagonalFilter<MatrixType>::getRangeMap() const
{
  return A_->getRangeMap();
}

template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
DiagonalFilter<MatrixType>::getGraph() const
{
  return A_->getGraph();
}

template<class MatrixType>
global_size_t DiagonalFilter<MatrixType>::getGlobalNumRows() const
{
  return A_->getGlobalNumRows();
}

template<class MatrixType>
global_size_t DiagonalFilter<MatrixType>::getGlobalNumCols() const
{
  return A_->getGlobalNumCols();
}

template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getNodeNumRows() const
{
  return A_->getNodeNumRows();
}

template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getNodeNumCols() const
{
  return A_->getNodeNumCols();
}

template<class MatrixType>
typename MatrixType::global_ordinal_type DiagonalFilter<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}

template<class MatrixType>
global_size_t DiagonalFilter<MatrixType>::getGlobalNumEntries() const
{
  return A_->getGlobalNumEntries();
}

template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getNodeNumEntries() const
{
  return A_->getNodeNumEntries();
}

template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
{
  return A_->getNumEntriesInGlobalRow(globalRow);
}

template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
{
  return A_->getNumEntriesInLocalRow(localRow);
}

template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return A_->getGlobalMaxNumRowEntries();
}

template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getNodeMaxNumRowEntries() const
{
  return A_->getNodeMaxNumRowEntries();
}

template<class MatrixType>
bool DiagonalFilter<MatrixType>::hasColMap() const
{
  return A_->hasColMap();
}

template<class MatrixType>
bool DiagonalFilter<MatrixType>::isLocallyIndexed() const
{
  return A_->isLocallyIndexed();
}

template<class MatrixType>
bool DiagonalFilter<MatrixType>::isGloballyIndexed() const
{
  return A_->isGloballyIndexed();
}

template<class MatrixType>
bool DiagonalFilter<MatrixType>::isFillComplete() const
{
  return A_->isFillComplete();
}

template<class MatrixType>
void DiagonalFilter<MatrixType>::
getGlobalRowCopy (GlobalOrdinal GlobalRow,
                  const Teuchos::ArrayView<GlobalOrdinal> &Indices,
                  const Teuchos::ArrayView<Scalar> &Values,
                  size_t &NumEntries) const
{
  Teuchos::ArrayRCP< const Scalar > myvals=val_->get1dView();
  LocalOrdinal LocalRow=getRowMap()->getLocalElement(GlobalRow);

  A_->getGlobalRowCopy(GlobalRow, Indices,Values,NumEntries);

  if (pos_[LocalRow] != -1)
    Values[pos_[LocalRow]] += myvals[LocalRow];
}

template<class MatrixType>
void DiagonalFilter<MatrixType>::
getLocalRowCopy (LocalOrdinal LocalRow,
                 const Teuchos::ArrayView<LocalOrdinal> &Indices,
                 const Teuchos::ArrayView<Scalar> &Values,
                 size_t &NumEntries) const
{
  Teuchos::ArrayRCP< const Scalar > myvals=val_->get1dView();

  A_->getLocalRowCopy(LocalRow, Indices,Values,NumEntries);

  if (pos_[LocalRow] != -1)
    Values[pos_[LocalRow]] += myvals[LocalRow];
}

template<class MatrixType>
void DiagonalFilter<MatrixType>::
getGlobalRowView (GlobalOrdinal /* GlobalRow */,
                  Teuchos::ArrayView<const GlobalOrdinal> &/* indices */,
                  Teuchos::ArrayView<const Scalar> &/* values */) const
{
  throw std::runtime_error("Ifpack2::DiagonalFilter: does not support getGlobalRowView.");
}

template<class MatrixType>
void DiagonalFilter<MatrixType>::
getLocalRowView (LocalOrdinal /* LocalRow */,
                 Teuchos::ArrayView<const LocalOrdinal> &/* indices */,
                 Teuchos::ArrayView<const Scalar> &/* values */) const
{
  throw std::runtime_error("Ifpack2::DiagonalFilter: does not support getLocalRowView.");
}

template<class MatrixType>
void DiagonalFilter<MatrixType>::getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const
{
  // Returns the matrix's actual diagonal, rather than the "filtered" diagonal.
  // This is dubious, but it duplicates the functionality of Old Ifpack.
  return A_->getLocalDiagCopy(diag);
}

template<class MatrixType>
void DiagonalFilter<MatrixType>::leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* x */)
{
  throw std::runtime_error("Ifpack2::DiagonalFilter does not support leftScale.");
}

template<class MatrixType>
void DiagonalFilter<MatrixType>::rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* x */)
{
  throw std::runtime_error("Ifpack2::DiagonalFilter does not support rightScale.");
}

template<class MatrixType>
void DiagonalFilter<MatrixType>::
apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
       Teuchos::ETransp mode,
       Scalar alpha,
       Scalar beta) const
{
  // FIXME (mfh 27 Jul 2015) If this is a "diagonal filter," why do we
  // need to apply the whole matrix???

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  A_->apply(X,Y,mode,alpha,beta);
  Y.elementWiseMultiply(one,*val_,X,one);
}

template<class MatrixType>
bool DiagonalFilter<MatrixType>::hasTransposeApply() const
{
  return A_->hasTransposeApply();
}

template<class MatrixType>
bool DiagonalFilter<MatrixType>::supportsRowViews() const
{
  return false;
}

template<class MatrixType>
typename DiagonalFilter<MatrixType>::mag_type DiagonalFilter<MatrixType>::getFrobeniusNorm() const
{
  throw std::runtime_error("Ifpack2::DiagonalFilter does not implement getFrobeniusNorm.");
}

} // namespace Ifpack2

#define IFPACK2_DIAGONALFILTER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::DiagonalFilter< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif
