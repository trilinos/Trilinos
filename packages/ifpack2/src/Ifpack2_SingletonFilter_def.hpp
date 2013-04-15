/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_SINGLETONFILTER_DEF_HPP
#define IFPACK2_SINGLETONFILTER_DEF_HPP
#include "Ifpack2_SingletonFilter_decl.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

namespace Ifpack2 {

//==========================================================================
template<class MatrixType>
SingletonFilter<MatrixType>::SingletonFilter(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix):
  A_(Matrix),
  NumSingletons_(0),
  NumRows_(0),
  NumNonzeros_(0),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0)
{

  // use this filter only on serial matrices
  if (A_->getComm()->getSize() != 1 || A_->getNodeNumRows() != A_->getGlobalNumRows()) {
    throw std::runtime_error("Ifpack2::SingeltonFilter can be used with Comm().getSize() == 1 only. This class is a tool for Ifpack2_AdditiveSchwarz, and it is not meant to be used otherwise.");
  }

  // Number of rows in A
  size_t NumRowsA_ = A_->getNodeNumRows();

  // tentative value for MaxNumEntries. This is the number of
  // nonzeros in the local matrix
  MaxNumEntriesA_ = A_->getNodeMaxNumRowEntries();

  // ExtractMyRowCopy() will use these vectors
  Indices_.resize(MaxNumEntriesA_);
  Values_.resize(MaxNumEntriesA_);

  // Initialize reordering vector to -1
  Reorder_.resize(NumRowsA_);
  Reorder_.assign(Reorder_.size(),-1);

  // first check how may singletons I do have
  NumRows_=0;
  for (size_t i = 0 ; i < NumRowsA_ ; ++i) {
    size_t Nnz;
    A_->getLocalRowCopy(i,Indices_,Values_,Nnz);
    if (Nnz != 1) {
      Reorder_[i] = NumRows_++;
    }
    else {
      NumSingletons_++;
    }
  }

  // build the inverse reordering
  InvReorder_.resize(NumRows_);
  for (size_t i = 0 ; i < NumRowsA_ ; ++i) {
    if (Reorder_[i] < 0)
      continue;
    InvReorder_[Reorder_[i]] = i;
  }
  NumEntries_.resize(NumRows_);
  SingletonIndex_.resize(NumSingletons_);

  
  // now compute the nonzeros per row
  size_t count = 0;
  for (size_t i = 0 ; i < NumRowsA_ ; ++i) {
    size_t Nnz;
    A_->getLocalRowCopy(i,Indices_,Values_,Nnz);
    LocalOrdinal ii = Reorder_[i];
    if (ii >= 0) {
      NumEntries_[ii] = Nnz;
      NumNonzeros_ += Nnz;
      if (Nnz > MaxNumEntries_)
	MaxNumEntries_ = Nnz;
    }
    else {
      SingletonIndex_[count] = i;
      count++;
    }
  }

  // Build the reduced map.  This map should be serial
  ReducedMap_ = Teuchos::rcp( new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(NumRows_,0,A_->getComm()) );

  // and finish up with the diagonal entry
  Diagonal_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(ReducedMap_) );

  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> DiagonalA(A_->getRowMap());
  A_->getLocalDiagCopy(DiagonalA);
  const Teuchos::ArrayRCP<const Scalar> & DiagonalAview = DiagonalA.get1dView();
  for (size_t i = 0 ; i < NumRows_ ; ++i) {
    LocalOrdinal ii = InvReorder_[i];
    Diagonal_->replaceLocalValue((LocalOrdinal)i,DiagonalAview[ii]);
  }
}

//=========================================================================
template<class MatrixType>
SingletonFilter<MatrixType>::~SingletonFilter() { }

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & SingletonFilter<MatrixType>::getComm() const
{
  return A_->getComm();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP <typename MatrixType::node_type> SingletonFilter<MatrixType>::getNode() const
{
  return A_->getNode();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
SingletonFilter<MatrixType>::getRowMap() const
{
  return ReducedMap_;
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
SingletonFilter<MatrixType>::getColMap() const
{
  return ReducedMap_;
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
SingletonFilter<MatrixType>::getDomainMap() const
{
  return ReducedMap_;
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
SingletonFilter<MatrixType>::getRangeMap() const
{
  return ReducedMap_;
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
SingletonFilter<MatrixType>::getGraph() const
{
  throw std::runtime_error("Ifpack2::SingletonFilter: does not support getGraph.");
}

//==========================================================================
template<class MatrixType>
global_size_t SingletonFilter<MatrixType>::getGlobalNumRows() const
{
  return NumRows_;
}

//==========================================================================
template<class MatrixType>
global_size_t SingletonFilter<MatrixType>::getGlobalNumCols() const
{
  return NumRows_;
}

//==========================================================================
template<class MatrixType>
size_t SingletonFilter<MatrixType>::getNodeNumRows() const
{
  return NumRows_;
}

//==========================================================================
 
template<class MatrixType>
size_t SingletonFilter<MatrixType>::getNodeNumCols() const
{
  return NumRows_;
}

//==========================================================================  
template<class MatrixType>
typename MatrixType::global_ordinal_type SingletonFilter<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}

//========================================================================== 
template<class MatrixType>
global_size_t SingletonFilter<MatrixType>::getGlobalNumEntries() const
{
  return NumNonzeros_;
}

//========================================================================== 
template<class MatrixType>
size_t SingletonFilter<MatrixType>::getNodeNumEntries() const
{
  return NumNonzeros_;
}

//==========================================================================
template<class MatrixType> 
size_t SingletonFilter<MatrixType>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
{
  throw std::runtime_error("Ifpack2::SingletonFilter does not implement getNumEntriesInGlobalRow.");
}

//==========================================================================  
template<class MatrixType> 
size_t SingletonFilter<MatrixType>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
{
  return NumEntries_[localRow];
}

//==========================================================================  
template<class MatrixType>   
global_size_t SingletonFilter<MatrixType>::getGlobalNumDiags() const
{
  return A_->getGlobalNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t SingletonFilter<MatrixType>::getNodeNumDiags() const
{
  return A_->getNodeNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t SingletonFilter<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return MaxNumEntries_;
}

//==========================================================================  
template<class MatrixType>   
size_t SingletonFilter<MatrixType>::getNodeMaxNumRowEntries() const
{
  return MaxNumEntries_;
}

//==========================================================================  
template<class MatrixType>   
bool SingletonFilter<MatrixType>::hasColMap() const
{
  return true;
}

//==========================================================================  
template<class MatrixType>   
bool SingletonFilter<MatrixType>::isLowerTriangular() const
{
  return A_->isLowerTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool SingletonFilter<MatrixType>::isUpperTriangular() const
{
  return A_->isUpperTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool SingletonFilter<MatrixType>::isLocallyIndexed() const
{
  return A_->isLocallyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool SingletonFilter<MatrixType>::isGloballyIndexed() const
{
  return A_->isGloballyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool SingletonFilter<MatrixType>::isFillComplete() const
{
  return A_->isFillComplete();
}
  
//==========================================================================
template<class MatrixType> 
void SingletonFilter<MatrixType>::getGlobalRowCopy(GlobalOrdinal GlobalRow,
						  const Teuchos::ArrayView<GlobalOrdinal> &Indices,
						  const Teuchos::ArrayView<Scalar> &Values,
						  size_t &NumEntries) const
{
  throw std::runtime_error("Ifpack2::SingletonFilter does not implement getGlobalRowCopy.");
}

//==========================================================================  
template<class MatrixType> 
void SingletonFilter<MatrixType>::getLocalRowCopy(LocalOrdinal LocalRow, 
					      const Teuchos::ArrayView<LocalOrdinal> &Indices, 
					      const Teuchos::ArrayView<Scalar> &Values,
					      size_t &NumEntries) const 
{ 
  TEUCHOS_TEST_FOR_EXCEPTION((LocalRow < 0 || (size_t) LocalRow >=  NumRows_ || (size_t) Indices.size() <  NumEntries_[LocalRow]), std::runtime_error, "Ifpack2::SingletonFilter::getLocalRowCopy invalid row or array size.");

  size_t Nnz;
  LocalOrdinal ARow = InvReorder_[LocalRow];
  A_->getLocalRowCopy(ARow,Indices_(),Values_(),Nnz);

  // populate the user's vectors
  NumEntries = 0;
  for (size_t i = 0 ; i < Nnz ; ++i) {
    LocalOrdinal ii = Reorder_[Indices_[i]];
    if ( ii >= 0) {
      Indices[NumEntries] = ii;
      Values[NumEntries] = Values_[i];
      NumEntries++;
    }
  }

}

//==========================================================================  
template<class MatrixType> 
void SingletonFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
						  Teuchos::ArrayView<const GlobalOrdinal> &indices, 
						  Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::SingletonFilter: does not support getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void SingletonFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow, 
						 Teuchos::ArrayView<const LocalOrdinal> &indices, 
						 Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::SingletonFilter: does not support getLocalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void SingletonFilter<MatrixType>::getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const
{
  // This is somewhat dubious as to how the maps match.
  return A_->getLocalDiagCopy(diag);
}

//========================================================================== 
template<class MatrixType> 
void SingletonFilter<MatrixType>::leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::SingletonFilter does not support leftScale.");
}

//==========================================================================  
template<class MatrixType> 
void SingletonFilter<MatrixType>::rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::SingletonFilter does not support rightScale.");
}

//==========================================================================  
template<class MatrixType> 
void SingletonFilter<MatrixType>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
				       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
				       Teuchos::ETransp mode, 
				       Scalar alpha,
				       Scalar beta) const
{  
  this->template applyTempl<Scalar,Scalar>(X, Y, mode, alpha, beta);
}


//==========================================================================  
template<class MatrixType> 
template<class DomainScalar, class RangeScalar>
void SingletonFilter<MatrixType>::applyTempl(const Tpetra::MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
					     Tpetra::MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
					     Teuchos::ETransp mode, 
					     RangeScalar alpha,
					     RangeScalar beta) const
{  
  // Note: This isn't AztecOO compliant.  But neither was Ifpack's version.

  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
			     "Ifpack2::SingletonFilter::apply ERROR: X.getNumVectors() != Y.getNumVectors().");
  
  RangeScalar zero = Teuchos::ScalarTraits<RangeScalar>::zero();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const DomainScalar> > x_ptr = X.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<RangeScalar> >        y_ptr = Y.get2dViewNonConst();

  Y.putScalar(zero);
  size_t NumVectors = Y.getNumVectors();


  for (size_t i = 0 ; i < NumRows_ ; ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    getLocalRowCopy(i,Indices_(),Values_(),Nnz);
    if (mode==Teuchos::NO_TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][i] += (RangeScalar)Values_[j] * (RangeScalar)x_ptr[k][Indices_[j]];      
    }
    else if (mode==Teuchos::TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][Indices_[j]] += (RangeScalar)Values_[j] * (RangeScalar)x_ptr[k][i];
    }
    else { //mode==Teuchos::CONJ_TRANS
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][Indices_[j]] += Teuchos::ScalarTraits<RangeScalar>::conjugate((RangeScalar)Values_[j]) * (RangeScalar)x_ptr[k][i];
    }
  }
}
  

//==========================================================================  
template<class MatrixType> 
bool SingletonFilter<MatrixType>::hasTransposeApply() const
{
  return true;
}

//==========================================================================  
template<class MatrixType> 
bool SingletonFilter<MatrixType>::supportsRowViews() const
{
  return false;
}

//==============================================================================
template<class MatrixType> 
void SingletonFilter<MatrixType>::SolveSingletons(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, 
						  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& LHS)
{
  this->template SolveSingletonsTempl<Scalar,Scalar>(RHS, LHS);
}

//==============================================================================
template<class MatrixType> 
template<class DomainScalar, class RangeScalar>
void SingletonFilter<MatrixType>::SolveSingletonsTempl(const Tpetra::MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, 
						       Tpetra::MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>& LHS)
{
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const DomainScalar> > RHS_ptr = RHS.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<RangeScalar> >        LHS_ptr = LHS.get2dViewNonConst();

  for (size_t i = 0 ; i < NumSingletons_ ; ++i) {
    LocalOrdinal ii = SingletonIndex_[i];
    // get the diagonal value for the singleton
    size_t Nnz;
    A_->getLocalRowCopy(ii,Indices_(),Values_(),Nnz);
    for (size_t j = 0 ; j < Nnz ; ++j) {
      if (Indices_[j] == ii) {
	for (size_t k = 0 ; k < LHS.getNumVectors() ; ++k)
	  LHS_ptr[k][ii] = (RangeScalar)RHS_ptr[k][ii] / (RangeScalar)Values_[j];
      }
    }
  }
}

//==============================================================================
template<class MatrixType> 
void SingletonFilter<MatrixType>::CreateReducedRHS(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& LHS,
						   const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, 
						   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& ReducedRHS)
{
  this->template CreateReducedRHSTempl<Scalar,Scalar>(LHS, RHS, ReducedRHS);
}

//==============================================================================
template<class MatrixType> 
template<class DomainScalar, class RangeScalar>
void SingletonFilter<MatrixType>::CreateReducedRHSTempl(const Tpetra::MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node>& LHS,
							const Tpetra::MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, 
							Tpetra::MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>& ReducedRHS)
{
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const RangeScalar > > RHS_ptr = RHS.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const DomainScalar> > LHS_ptr = LHS.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<RangeScalar> >        ReducedRHS_ptr = ReducedRHS.get2dViewNonConst();

  size_t NumVectors = LHS.getNumVectors();

  for (size_t i = 0 ; i < NumRows_ ; ++i)
    for (size_t k = 0 ; k < NumVectors ; ++k)
      ReducedRHS_ptr[k][i] = RHS_ptr[k][InvReorder_[i]];

  for (size_t i = 0 ; i < NumRows_ ; ++i) {
    LocalOrdinal ii = InvReorder_[i];
    size_t Nnz; 
    A_->getLocalRowCopy(ii,Indices_(),Values_(),Nnz);

    for (size_t j = 0 ; j < Nnz ; ++j) {
      if (Reorder_[Indices_[j]] == -1) {
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  ReducedRHS_ptr[k][i] -= (RangeScalar)Values_[j] * (RangeScalar)LHS_ptr[k][Indices_[j]];
      }
    }
  }
}

//==============================================================================
template<class MatrixType> 
void SingletonFilter<MatrixType>::UpdateLHS(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& ReducedLHS,
					    Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& LHS)
{
  this->template UpdateLHSTempl<Scalar,Scalar>(ReducedLHS, LHS);
}

//==============================================================================
template<class MatrixType> 
template<class DomainScalar, class RangeScalar>
void SingletonFilter<MatrixType>::UpdateLHSTempl(const Tpetra::MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node>& ReducedLHS,
						 Tpetra::MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>& LHS)
{
  
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<RangeScalar> >        LHS_ptr = LHS.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const DomainScalar> >  ReducedLHS_ptr = ReducedLHS.get2dView();

  for (size_t i = 0 ; i < NumRows_ ; ++i)
    for (size_t k = 0 ; k < LHS.getNumVectors() ; ++k)
      LHS_ptr[k][InvReorder_[i]] = (RangeScalar)ReducedLHS_ptr[k][i];
}

//==========================================================================  
template<class MatrixType> 
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType SingletonFilter<MatrixType>::getFrobeniusNorm() const
{
  throw std::runtime_error("Ifpack2::SingletonFilter does not implement getFrobeniusNorm.");
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void SingletonFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
								     Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
								     Teuchos::ArrayRCP<const Scalar>        &values) const
{
  throw std::runtime_error("Ifpack2::SingletonFilter does not implement getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void SingletonFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow,
								    Teuchos::ArrayRCP<const LocalOrdinal> &indices,
								    Teuchos::ArrayRCP<const Scalar>       &values) const
{
  throw std::runtime_error("Ifpack2::SingletonFilter does not implement getLocalRowView.");
}




}// namespace Ifpack2

#endif
