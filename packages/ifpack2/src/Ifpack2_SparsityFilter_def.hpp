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

#ifndef IFPACK2_SPARSITYFILTER_DEF_HPP
#define IFPACK2_SPARSITYFILTER_DEF_HPP
#include "Ifpack2_SparsityFilter_decl.hpp"
#include <vector>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

#include <vector>
#include <algorithm>

namespace Ifpack2 {

//==========================================================================
template<class MatrixType>
SparsityFilter<MatrixType>::SparsityFilter(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix, 
					   size_t AllowedNumEntries,
					   LocalOrdinal AllowedBandwidth):
  A_(Matrix),
  AllowedNumEntries_(AllowedNumEntries),
  AllowedBandwidth_(AllowedBandwidth),
  NumRows_(0),
  NumNonzeros_(0),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0)
{

 // use this filter only on serial matrices
  if (A_->getComm()->getSize() != 1 || A_->getNodeNumRows() != A_->getGlobalNumRows()) {
    throw std::runtime_error("Ifpack2::SparsityFilter can be used with Comm().getSize() == 1 only. This class is a tool for Ifpack2_AdditiveSchwarz, and it is not meant to be used otherwise.");
  }


  // localized matrix has all the local rows of Matrix
  NumRows_ = A_->getNodeNumRows();

  // NodeNumEntries_ will contain the actual number of nonzeros
  // for each localized row (that is, without external nodes,
  // and always with the diagonal entry)
  NumEntries_.resize(NumRows_);

  // tentative value for MaxNumEntries. This is the number of
  // nonzeros in the local matrix
  MaxNumEntries_  = A_->getNodeMaxNumRowEntries();
  MaxNumEntriesA_ = A_->getNodeMaxNumRowEntries();

  // ExtractMyRowCopy() will use these vectors
  Indices_.resize(MaxNumEntries_);
  Values_.resize(MaxNumEntries_);

  size_t ActualMaxNumEntries = 0;
  for (size_t i = 0 ; i < NumRows_ ; ++i) {
    NumEntries_[i] = MaxNumEntriesA_;
    size_t Nnz;
    A_->getLocalRowCopy(i,Indices_,Values_,Nnz);

    NumNonzeros_ += Nnz;
    NumEntries_[i] = Nnz;
    if (Nnz > ActualMaxNumEntries)
      ActualMaxNumEntries = Nnz;
  }
  
  MaxNumEntries_ = ActualMaxNumEntries;
}

//=========================================================================
template<class MatrixType>
SparsityFilter<MatrixType>::~SparsityFilter() { }

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & SparsityFilter<MatrixType>::getComm() const
{
  return A_->getComm();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP <typename MatrixType::node_type> SparsityFilter<MatrixType>::getNode() const
{
  return A_->getNode();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
SparsityFilter<MatrixType>::getRowMap() const
{
  return A_->getRowMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
SparsityFilter<MatrixType>::getColMap() const
{
  return A_->getColMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
SparsityFilter<MatrixType>::getDomainMap() const
{
  return A_->getDomainMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
SparsityFilter<MatrixType>::getRangeMap() const
{
  return A_->getRangeMap();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
SparsityFilter<MatrixType>::getGraph() const
{
  throw std::runtime_error("Ifpack2::SparsityFilter: does not support getGraph.");
}

//==========================================================================
template<class MatrixType>
global_size_t SparsityFilter<MatrixType>::getGlobalNumRows() const
{
  return NumRows_;
}

//==========================================================================
template<class MatrixType>
global_size_t SparsityFilter<MatrixType>::getGlobalNumCols() const
{
  return NumRows_;
}

//==========================================================================
template<class MatrixType>
size_t SparsityFilter<MatrixType>::getNodeNumRows() const
{
  return NumRows_;
}

//==========================================================================
 
template<class MatrixType>
size_t SparsityFilter<MatrixType>::getNodeNumCols() const
{
  return NumRows_;
}

//==========================================================================  
template<class MatrixType>
typename MatrixType::global_ordinal_type SparsityFilter<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}

//========================================================================== 
template<class MatrixType>
global_size_t SparsityFilter<MatrixType>::getGlobalNumEntries() const
{
  return NumNonzeros_;
}

//========================================================================== 
template<class MatrixType>
size_t SparsityFilter<MatrixType>::getNodeNumEntries() const
{
  return NumNonzeros_;
}

//==========================================================================
template<class MatrixType> 
size_t SparsityFilter<MatrixType>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
{
  throw std::runtime_error("Ifpack2::SparsityFilter does not implement getNumEntriesInGlobalRow.");
}

//==========================================================================  
template<class MatrixType> 
size_t SparsityFilter<MatrixType>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
{
  return NumEntries_[localRow];
}

//==========================================================================  
template<class MatrixType>   
global_size_t SparsityFilter<MatrixType>::getGlobalNumDiags() const
{
  return A_->getGlobalNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t SparsityFilter<MatrixType>::getNodeNumDiags() const
{
  return A_->getNodeNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t SparsityFilter<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return MaxNumEntries_;
}

//==========================================================================  
template<class MatrixType>   
size_t SparsityFilter<MatrixType>::getNodeMaxNumRowEntries() const
{
  return MaxNumEntries_;
}

//==========================================================================  
template<class MatrixType>   
bool SparsityFilter<MatrixType>::hasColMap() const
{
  return true;
}

//==========================================================================  
template<class MatrixType>   
bool SparsityFilter<MatrixType>::isLowerTriangular() const
{
  return A_->isLowerTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool SparsityFilter<MatrixType>::isUpperTriangular() const
{
  return A_->isUpperTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool SparsityFilter<MatrixType>::isLocallyIndexed() const
{
  return A_->isLocallyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool SparsityFilter<MatrixType>::isGloballyIndexed() const
{
  return A_->isGloballyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool SparsityFilter<MatrixType>::isFillComplete() const
{
  return A_->isFillComplete();
}
  
//==========================================================================
template<class MatrixType> 
void SparsityFilter<MatrixType>::getGlobalRowCopy(GlobalOrdinal GlobalRow,
						  const Teuchos::ArrayView<GlobalOrdinal> &Indices,
						  const Teuchos::ArrayView<Scalar> &Values,
						  size_t &NumEntries) const
{
  throw std::runtime_error("Ifpack2::SparsityFilter does not implement getGlobalRowCopy.");
}

//==========================================================================  
template<class MatrixType> 
void SparsityFilter<MatrixType>::getLocalRowCopy(LocalOrdinal LocalRow, 
					      const Teuchos::ArrayView<LocalOrdinal> &Indices, 
					      const Teuchos::ArrayView<Scalar> &Values,
					      size_t &NumEntries) const 
{ 
  TEUCHOS_TEST_FOR_EXCEPTION((LocalRow < 0 || (size_t) LocalRow >=  NumRows_ || (size_t) Indices.size() <  NumEntries_[LocalRow]), std::runtime_error, "Ifpack2::SparsityFilter::getLocalRowCopy invalid row or array size.");

  // Note: This function will work correctly if called by apply, say, with Indices_ and Values_ as
  // parameters.  The structure of the loop below should make that obvious.

  // always extract using the object Values_ and Indices_.
  // This is because I need more space than that given by
  // the user (for the external nodes)
  size_t A_NumEntries=0;
  A_->getLocalRowCopy(LocalRow,Indices_(),Values_(),A_NumEntries);
  magnitudeType Threshold = Teuchos::ScalarTraits<magnitudeType>::zero();
  std::vector<magnitudeType> Values2(A_NumEntries,Teuchos::ScalarTraits<magnitudeType>::zero());

  if (A_NumEntries > AllowedNumEntries_) {
    size_t count = 0;
    for (size_t i = 0 ; i < A_NumEntries ; ++i) {
      // skip diagonal entry (which is always inserted)
      if (Indices_[i] == LocalRow)
	continue;
      // put magnitude
      Values2[count] = Teuchos::ScalarTraits<Scalar>::magnitude(Values_[i]);
      count++;
    }

    // sort in descending order
    std::sort(Values2.rbegin(),Values2.rend());
    // get the cut-off value
    Threshold = Values2[AllowedNumEntries_ - 1];
  }
  

  // loop over all nonzero elements of row MyRow,
  // and drop elements below specified threshold.
  // Also, add value to zero diagonal
  NumEntries = 0;
  for (size_t i = 0 ; i < A_NumEntries ; ++i) {
    if (std::abs(Indices_[i] - LocalRow) > AllowedBandwidth_)
      continue;
    
    if ((Indices_[i] != LocalRow) && (Teuchos::ScalarTraits<Scalar>::magnitude(Values_[i]) < Threshold))
      continue;

    Values[NumEntries] = Values_[i];
    Indices[NumEntries] = Indices_[i];

    NumEntries++;
    if (NumEntries > AllowedNumEntries_)
      break;
  }

}

//==========================================================================  
template<class MatrixType> 
void SparsityFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
						  Teuchos::ArrayView<const GlobalOrdinal> &indices, 
						  Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::SparsityFilter: does not support getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void SparsityFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow, 
						 Teuchos::ArrayView<const LocalOrdinal> &indices, 
						 Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::SparsityFilter: does not support getLocalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void SparsityFilter<MatrixType>::getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const
{
  // This is somewhat dubious as to how the maps match.
  return A_->getLocalDiagCopy(diag);
}

//========================================================================== 
template<class MatrixType> 
void SparsityFilter<MatrixType>::leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::SparsityFilter does not support leftScale.");
}

//==========================================================================  
template<class MatrixType> 
void SparsityFilter<MatrixType>::rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::SparsityFilter does not support rightScale.");
}

//==========================================================================  
template<class MatrixType> 
void SparsityFilter<MatrixType>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
				       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
				       Teuchos::ETransp mode, 
				       Scalar alpha,
				       Scalar beta) const
{  
  // Note: This isn't AztecOO compliant.  But neither was Ifpack's version.
  // Note: The localized maps mean the matvec is trivial (and has no import)
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Ifpack2::SparsityFilter::apply ERROR: X.getNumVectors() != Y.getNumVectors().");
 
  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr = X.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       y_ptr = Y.get2dViewNonConst();

  Y.putScalar(zero);
  size_t NumVectors = Y.getNumVectors();

  for (size_t i = 0 ; i < NumRows_ ; ++i) {
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
bool SparsityFilter<MatrixType>::hasTransposeApply() const
{
  return true;
}

//==========================================================================
template<class MatrixType> 
bool SparsityFilter<MatrixType>::supportsRowViews() const
{
  return false;
}

//==========================================================================  
template<class MatrixType> 
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType SparsityFilter<MatrixType>::getFrobeniusNorm() const
{
  throw std::runtime_error("Ifpack2::SparsityFilter does not implement getFrobeniusNorm.");
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void SparsityFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
								     Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
								     Teuchos::ArrayRCP<const Scalar>        &values) const
{
  throw std::runtime_error("Ifpack2::SparsityFilter does not implement getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void SparsityFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow,
								    Teuchos::ArrayRCP<const LocalOrdinal> &indices,
								    Teuchos::ArrayRCP<const Scalar>       &values) const
{
  throw std::runtime_error("Ifpack2::SparsityFilter does not implement getLocalRowView.");
}

}// namespace Ifpack2

#endif
