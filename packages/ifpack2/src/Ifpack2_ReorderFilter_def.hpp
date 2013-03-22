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

//==========================================================================

#ifdef HAVE_IFPACK2_ZOLTAN2
template<class MatrixType>
ReorderFilter<MatrixType>::ReorderFilter(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & Matrix,
					 const Teuchos::RCP<const Zoltan2::OrderingSolution<GlobalOrdinal,LocalOrdinal> > & Reordering):
  A_(Matrix)
{

 // use this filter only on serial matrices
  if (A_->getComm()->getSize() != 1 || A_->getNodeNumRows() != A_->getGlobalNumRows()) {
    throw std::runtime_error("Ifpack2::ReorderFilter can be used with Comm().getSize() == 1 only. This class is a tool for Ifpack2_AdditiveSchwarz, and it is not meant to be used otherwise.");
  }

  // perm_[i]         gives the where OLD index i shows up in the NEW ordering
  // reverseperm_[i]  gives the where NEW index i shows up in the OLD ordering
  perm_=Reordering->getPermutationRCPConst();

  // Generate the reverse permutation
  size_t N=A_->getNodeNumRows();
  reverseperm_.resize(N);

  for(size_t i=0; i<N; i++) {
    reverseperm_[perm_[i]] = i;
  }

  // Temp arrays for apply
  Indices_.resize(A_->getNodeMaxNumRowEntries());
  Values_.resize(A_->getNodeMaxNumRowEntries());
}
#endif

//=========================================================================
template<class MatrixType>
ReorderFilter<MatrixType>::~ReorderFilter() { }

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & ReorderFilter<MatrixType>::getComm() const
{
  return A_->getComm();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP <typename MatrixType::node_type> ReorderFilter<MatrixType>::getNode() const
{
  return A_->getNode();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
ReorderFilter<MatrixType>::getRowMap() const
{
  return A_->getRowMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
ReorderFilter<MatrixType>::getColMap() const
{
  return A_->getColMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
ReorderFilter<MatrixType>::getDomainMap() const
{
  return A_->getDomainMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
ReorderFilter<MatrixType>::getRangeMap() const
{
  return A_->getRangeMap();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
ReorderFilter<MatrixType>::getGraph() const
{
  throw std::runtime_error("Ifpack2::ReorderFilter: does not support getGraph.");
}

//==========================================================================
template<class MatrixType>
global_size_t ReorderFilter<MatrixType>::getGlobalNumRows() const
{
  return A_->getGlobalNumRows();
}

//==========================================================================
template<class MatrixType>
global_size_t ReorderFilter<MatrixType>::getGlobalNumCols() const
{
  return A_->getGlobalNumCols();
}

//==========================================================================
template<class MatrixType>
size_t ReorderFilter<MatrixType>::getNodeNumRows() const
{
  return A_->getNodeNumRows();
}

//==========================================================================
 
template<class MatrixType>
size_t ReorderFilter<MatrixType>::getNodeNumCols() const
{
  return A_->getNodeNumCols();
}

//==========================================================================  
template<class MatrixType>
typename MatrixType::global_ordinal_type ReorderFilter<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}

//========================================================================== 
template<class MatrixType>
global_size_t ReorderFilter<MatrixType>::getGlobalNumEntries() const
{
  return A_->getGlobalNumEntries();
}

//========================================================================== 
template<class MatrixType>
size_t ReorderFilter<MatrixType>::getNodeNumEntries() const
{
  return A_->getNodeNumEntries();
}

//==========================================================================
template<class MatrixType> 
size_t ReorderFilter<MatrixType>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not implement getNumEntriesInGlobalRow");
}

//==========================================================================  
template<class MatrixType> 
size_t ReorderFilter<MatrixType>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
{
  LocalOrdinal MyReorderedRow = reverseperm_[localRow];
  return A_->getNumEntriesInLocalRow(MyReorderedRow);
}

//==========================================================================  
template<class MatrixType>   
global_size_t ReorderFilter<MatrixType>::getGlobalNumDiags() const
{
  return A_->getGlobalNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t ReorderFilter<MatrixType>::getNodeNumDiags() const
{
  return A_->getNodeNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t ReorderFilter<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return A_->getGlobalMaxNumRowEntries();
}

//==========================================================================  
template<class MatrixType>   
size_t ReorderFilter<MatrixType>::getNodeMaxNumRowEntries() const
{
  return A_->getNodeMaxNumRowEntries();
}

//==========================================================================  
template<class MatrixType>   
bool ReorderFilter<MatrixType>::hasColMap() const
{
  return true;
}

//==========================================================================  
template<class MatrixType>   
bool ReorderFilter<MatrixType>::isLowerTriangular() const
{
  return A_->isLowerTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool ReorderFilter<MatrixType>::isUpperTriangular() const
{
  return A_->isUpperTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool ReorderFilter<MatrixType>::isLocallyIndexed() const
{
  return A_->isLocallyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool ReorderFilter<MatrixType>::isGloballyIndexed() const
{
  return A_->isGloballyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool ReorderFilter<MatrixType>::isFillComplete() const
{
  return A_->isFillComplete();
}
  
//==========================================================================
template<class MatrixType> 
void ReorderFilter<MatrixType>::getGlobalRowCopy(GlobalOrdinal GlobalRow,
						  const Teuchos::ArrayView<GlobalOrdinal> &Indices,
						  const Teuchos::ArrayView<Scalar> &Values,
						  size_t &NumEntries) const
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not implement getGlobalRowCopy.");
}

//==========================================================================  
template<class MatrixType> 
void ReorderFilter<MatrixType>::getLocalRowCopy(LocalOrdinal LocalRow, 
					      const Teuchos::ArrayView<LocalOrdinal> &Indices, 
					      const Teuchos::ArrayView<Scalar> &Values,
					      size_t &NumEntries) const 
{ 
  TEUCHOS_TEST_FOR_EXCEPTION((LocalRow < 0 || (size_t) LocalRow >=  A_->getNodeNumRows() || (size_t) Indices.size() <  A_->getNumEntriesInLocalRow(LocalRow)), std::runtime_error, "Ifpack2::ReorderFilter::getLocalRowCopy invalid row or array size.");


  LocalOrdinal MyOriginalRow = reverseperm_[LocalRow];
  A_->getLocalRowCopy(MyOriginalRow,Indices,Values,NumEntries);
  // Do a col reindex via perm
  for (size_t i = 0 ; i < NumEntries ; ++i) {
    Indices[i] = perm_[Indices[i]];
  }
}

//==========================================================================  
template<class MatrixType> 
void ReorderFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
						  Teuchos::ArrayView<const GlobalOrdinal> &indices, 
						  Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::ReorderFilter: does not support getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void ReorderFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow, 
						 Teuchos::ArrayView<const LocalOrdinal> &indices, 
						 Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::ReorderFilter: does not support getLocalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void ReorderFilter<MatrixType>::getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const
{
  // This is somewhat dubious as to how the maps match.
  return A_->getLocalDiagCopy(diag);
}

//========================================================================== 
template<class MatrixType> 
void ReorderFilter<MatrixType>::leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not support leftScale.");
}

//==========================================================================  
template<class MatrixType> 
void ReorderFilter<MatrixType>::rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not support rightScale.");
}

//==========================================================================  
template<class MatrixType> 
void ReorderFilter<MatrixType>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
				       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
				       Teuchos::ETransp mode, 
				       Scalar alpha,
				       Scalar beta) const
{  
  // Note: This isn't AztecOO compliant.  But neither was Ifpack's version.
  // Note: The localized maps mean the matvec is trivial (and has no import)
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Ifpack2::ReorderFilter::apply ERROR: X.getNumVectors() != Y.getNumVectors().");
 
  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr = X.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       y_ptr = Y.get2dViewNonConst();

  Y.putScalar(zero);
  size_t NumVectors = Y.getNumVectors();

  for (size_t i = 0 ; i < A_->getNodeNumRows() ; ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    getLocalRowCopy(i,Indices_(),Values_(),Nnz);
    if (mode==Teuchos::NO_TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][i] += Values_[j] * x_ptr[k][Indices_[j]];      
    }
    else if (mode==Teuchos::TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][Indices_[j]] += Values_[j] * x_ptr[k][i];
    }
    else { //mode==Teuchos::CONJ_TRANS
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][Indices_[j]] += Teuchos::ScalarTraits<Scalar>::conjugate(Values_[j]) * x_ptr[k][i];
    }
  }
}
  
//==========================================================================  
template<class MatrixType> 
bool ReorderFilter<MatrixType>::hasTransposeApply() const
{
  return true;
}

//==========================================================================  
template<class MatrixType> 
bool ReorderFilter<MatrixType>::supportsRowViews() const
{
  return false;
}

//==========================================================================  
template<class MatrixType> 
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType ReorderFilter<MatrixType>::getFrobeniusNorm() const
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not implement getFrobeniusNorm.");
}

//==========================================================================  
//! Permute multivector: original-to-reordered
template<class MatrixType> 
void ReorderFilter<MatrixType>::permuteOriginalToReordered(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &originalX, 
							   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &reorderedY) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(originalX.getNumVectors() != reorderedY.getNumVectors(), std::runtime_error,
			     "Ifpack2::ReorderFilter::permuteOriginalToReordered ERROR: X.getNumVectors() != Y.getNumVectors().");

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr = originalX.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       y_ptr = reorderedY.get2dViewNonConst();

  for(size_t k=0; k < originalX.getNumVectors(); k++)
    for(LocalOrdinal i=0; (size_t)i< originalX.getLocalLength(); i++)
      y_ptr[k][perm_[i]] = x_ptr[k][i];
}							  
  
//==========================================================================  
//! Permute multivector: reordered-to-original
template<class MatrixType> 
void ReorderFilter<MatrixType>::permuteReorderedToOriginal(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &reorderedX, 
							   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &originalY) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(reorderedX.getNumVectors() != originalY.getNumVectors(), std::runtime_error,
			     "Ifpack2::ReorderFilter::permuteReorderedToOriginal ERROR: X.getNumVectors() != Y.getNumVectors().");

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr = reorderedX.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       y_ptr = originalY.get2dViewNonConst();

  for(size_t k=0; k < reorderedX.getNumVectors(); k++)
    for(LocalOrdinal i=0; (size_t)i< reorderedX.getLocalLength(); i++)
      y_ptr[k][reverseperm_[i]] = x_ptr[k][i];
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void ReorderFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
								     Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
								     Teuchos::ArrayRCP<const Scalar>        &values) const
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not implement getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void ReorderFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow,
								    Teuchos::ArrayRCP<const LocalOrdinal> &indices,
								    Teuchos::ArrayRCP<const Scalar>       &values) const
{
  throw std::runtime_error("Ifpack2::ReorderFilter does not implement getLocalRowView.");
}

}// namespace Ifpack2

#endif
