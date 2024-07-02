// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  pos_.resize(getLocalNumRows());
  val_=Teuchos::rcp(new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A_->getRowMap()));

  nonconst_local_inds_host_view_type Indices("Indices",getLocalMaxNumRowEntries());
  nonconst_values_host_view_type Values("Values",getLocalMaxNumRowEntries());
  size_t NumEntries;
  magnitudeType mysign;


  for (size_t MyRow = 0 ; MyRow < getLocalNumRows() ; ++MyRow) {
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
size_t DiagonalFilter<MatrixType>::getLocalNumRows() const
{
  return A_->getLocalNumRows();
}

template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getLocalNumCols() const
{
  return A_->getLocalNumCols();
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
size_t DiagonalFilter<MatrixType>::getLocalNumEntries() const
{
  return A_->getLocalNumEntries();
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
size_t DiagonalFilter<MatrixType>::getLocalMaxNumRowEntries() const
{
  return A_->getLocalMaxNumRowEntries();
}

template<class MatrixType>
typename MatrixType::local_ordinal_type DiagonalFilter<MatrixType>::getBlockSize() const
{
  return A_->getBlockSize();
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
                   nonconst_global_inds_host_view_type &Indices,
                   nonconst_values_host_view_type &Values,
                   size_t& NumEntries) const

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
                   nonconst_local_inds_host_view_type &Indices,
                   nonconst_values_host_view_type &Values,
                   size_t& NumEntries) const
{
  Teuchos::ArrayRCP< const Scalar > myvals=val_->get1dView();

  A_->getLocalRowCopy(LocalRow, Indices,Values,NumEntries);

  if (pos_[LocalRow] != -1)
    Values[pos_[LocalRow]] += myvals[LocalRow];
}


template<class MatrixType>
void DiagonalFilter<MatrixType>::getGlobalRowView(GlobalOrdinal /* GlobalRow */,
                                                  global_inds_host_view_type &/*indices*/,
                                                  values_host_view_type &/*values*/) const
{
  throw std::runtime_error("Ifpack2::DiagonalFilter: does not support getGlobalRowView.");
}


template<class MatrixType>
void DiagonalFilter<MatrixType>::getLocalRowView(LocalOrdinal /* LocalRow */,
    local_inds_host_view_type & /*indices*/,
    values_host_view_type & /*values*/) const
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
