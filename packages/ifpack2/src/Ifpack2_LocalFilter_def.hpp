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

#ifndef IFPACK2_LOCALFILTER_DEF_HPP
#define IFPACK2_LOCALFILTER_DEF_HPP
#include "Ifpack2_LocalFilter_decl.hpp"
#include <vector>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif
namespace Ifpack2 {

//==========================================================================
template<class MatrixType>
LocalFilter<MatrixType>::LocalFilter(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix):
  A_(Matrix),
  NumRows_(0),
  NumNonzeros_(0),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0)
{

#ifdef HAVE_MPI
  LocalComm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper((MPI_Comm)MPI_COMM_SELF)));
#else
  LocalComm_ = Teuchos::rcp( new Teuchos::SerialComm<int>() );
#endif

  // localized matrix has all the local rows of Matrix
  NumRows_ = A_->getNodeNumRows();

  // build a linear map, based on the serial communicator
  LocalMap_ = Teuchos::rcp( new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(NumRows_,0,LocalComm_) );

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

  // now compute:
  // - the number of nonzero per row
  // - the total number of nonzeros
  // - the diagonal entries

  // compute nonzeros (total and per-row), and store the
  // diagonal entries (already modified)
  size_t ActualMaxNumEntries = 0;

  for (size_t i = 0 ; i < NumRows_ ; ++i) {
    
    NumEntries_[i] = 0;
    size_t Nnz, NewNnz = 0;
    A_->getLocalRowCopy(i,Indices_,Values_,Nnz);
    for (size_t j = 0 ; j < Nnz ; ++j) {
      if ((size_t) Indices_[j] < NumRows_ ) ++NewNnz;
    }

    if (NewNnz > ActualMaxNumEntries)
      ActualMaxNumEntries = NewNnz;

    NumNonzeros_ += NewNnz;
    NumEntries_[i] = NewNnz;

  }
 
  MaxNumEntries_ = ActualMaxNumEntries;
}

//=========================================================================
template<class MatrixType>
LocalFilter<MatrixType>::~LocalFilter() { }

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & LocalFilter<MatrixType>::getComm() const
{
  return LocalComm_;
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP <typename MatrixType::node_type> LocalFilter<MatrixType>::getNode() const
{
  return A_->getNode();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
LocalFilter<MatrixType>::getRowMap() const
{
  return LocalMap_;
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
LocalFilter<MatrixType>::getColMap() const
{
  return LocalMap_;
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
LocalFilter<MatrixType>::getDomainMap() const
{
  return LocalMap_;
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
LocalFilter<MatrixType>::getRangeMap() const
{
  return LocalMap_;
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
LocalFilter<MatrixType>::getGraph() const
{
  throw std::runtime_error("Ifpack2::LocalFilter: does not support getGraph.");
}

//==========================================================================
template<class MatrixType>
global_size_t LocalFilter<MatrixType>::getGlobalNumRows() const
{
  return NumRows_;
}

//==========================================================================
template<class MatrixType>
global_size_t LocalFilter<MatrixType>::getGlobalNumCols() const
{
  return NumRows_;
}

//==========================================================================
template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeNumRows() const
{
  return NumRows_;
}

//==========================================================================
 
template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeNumCols() const
{
  return NumRows_;
}

//==========================================================================  
template<class MatrixType>
typename MatrixType::global_ordinal_type LocalFilter<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}

//========================================================================== 
template<class MatrixType>
global_size_t LocalFilter<MatrixType>::getGlobalNumEntries() const
{
  return NumNonzeros_;
}

//========================================================================== 
template<class MatrixType>
size_t LocalFilter<MatrixType>::getNodeNumEntries() const
{
  return NumNonzeros_;
}

//==========================================================================
template<class MatrixType> 
size_t LocalFilter<MatrixType>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
{
  throw std::runtime_error("Ifpack2::LocalFilter does not implement getNumEntriesInGlobalRow.");
}

//==========================================================================  
template<class MatrixType> 
size_t LocalFilter<MatrixType>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
{
  return NumEntries_[localRow];
}

//==========================================================================  
template<class MatrixType>   
global_size_t LocalFilter<MatrixType>::getGlobalNumDiags() const
{
  return A_->getGlobalNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t LocalFilter<MatrixType>::getNodeNumDiags() const
{
  return A_->getNodeNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t LocalFilter<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return MaxNumEntries_;
}

//==========================================================================  
template<class MatrixType>   
size_t LocalFilter<MatrixType>::getNodeMaxNumRowEntries() const
{
  return MaxNumEntries_;
}

//==========================================================================  
template<class MatrixType>   
bool LocalFilter<MatrixType>::hasColMap() const
{
  return true;
}

//==========================================================================  
template<class MatrixType>   
bool LocalFilter<MatrixType>::isLowerTriangular() const
{
  return A_->isLowerTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool LocalFilter<MatrixType>::isUpperTriangular() const
{
  return A_->isUpperTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool LocalFilter<MatrixType>::isLocallyIndexed() const
{
  return A_->isLocallyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool LocalFilter<MatrixType>::isGloballyIndexed() const
{
  return A_->isGloballyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool LocalFilter<MatrixType>::isFillComplete() const
{
  return A_->isFillComplete();
}
  
//==========================================================================
template<class MatrixType> 
void LocalFilter<MatrixType>::getGlobalRowCopy(GlobalOrdinal GlobalRow,
						  const Teuchos::ArrayView<GlobalOrdinal> &Indices,
						  const Teuchos::ArrayView<Scalar> &Values,
						  size_t &NumEntries) const
{
  throw std::runtime_error("Ifpack2::LocalFilter does not implement getGlobalRowCopy.");
}

//==========================================================================  
template<class MatrixType> 
void LocalFilter<MatrixType>::getLocalRowCopy(LocalOrdinal LocalRow, 
					      const Teuchos::ArrayView<LocalOrdinal> &Indices, 
					      const Teuchos::ArrayView<Scalar> &Values,
					      size_t &NumEntries) const 
{ 
  TEUCHOS_TEST_FOR_EXCEPTION((LocalRow < 0 || (size_t) LocalRow >=  NumRows_ || (size_t) Indices.size() <  NumEntries_[LocalRow]), std::runtime_error, "Ifpack2::LocalFilter::getLocalRowCopy invalid row or array size.");

  size_t A_NumEntries=0;
  // always extract using the object Values_ and Indices_.
  // This is because I need more space than that given by
  // the user (for the external nodes)
  A_->getLocalRowCopy(LocalRow,Indices_(),Values_(),A_NumEntries);

  // populate the user's vectors
  NumEntries=0;
  for (size_t j = 0 ; j < A_NumEntries; ++j) {
    // only local indices
    if ((size_t)Indices_[j] < NumRows_ ) {
      Indices[NumEntries] = Indices_[j];
      Values[NumEntries]  = Values_[j];
      NumEntries++;
    }
  }

}

//==========================================================================  
template<class MatrixType> 
void LocalFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
						  Teuchos::ArrayView<const GlobalOrdinal> &indices, 
						  Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::LocalFilter: does not support getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void LocalFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow, 
						 Teuchos::ArrayView<const LocalOrdinal> &indices, 
						 Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::LocalFilter: does not support getLocalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void LocalFilter<MatrixType>::getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const
{
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> temp(A_->getRowMap());
  A_->getLocalDiagCopy(temp);
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       d_ptr = diag.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > t_ptr = temp.get2dView();

  for(size_t i=0; i<NumRows_; i++)
    d_ptr[0][i] = t_ptr[0][i];
}

//========================================================================== 
template<class MatrixType> 
void LocalFilter<MatrixType>::leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::LocalFilter does not support leftScale.");
}

//==========================================================================  
template<class MatrixType> 
void LocalFilter<MatrixType>::rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::LocalFilter does not support rightScale.");
}

//==========================================================================  
template<class MatrixType> 
void LocalFilter<MatrixType>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
				       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
				       Teuchos::ETransp mode, 
				       Scalar alpha,
				       Scalar beta) const
{  
  // Note: This isn't AztecOO compliant.  But neither was Ifpack's version.

  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
			     "Ifpack2::LocalFilter::apply ERROR: X.getNumVectors() != Y.getNumVectors().");

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
bool LocalFilter<MatrixType>::hasTransposeApply() const
{
  return true;
}

//==========================================================================  
template<class MatrixType> 
bool LocalFilter<MatrixType>::supportsRowViews() const
{
  return false;
}

//==========================================================================  
template<class MatrixType> 
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType LocalFilter<MatrixType>::getFrobeniusNorm() const
{
  throw std::runtime_error("Ifpack2::LocalFilter does not implement getFrobeniusNorm.");
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void LocalFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
								     Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
								     Teuchos::ArrayRCP<const Scalar>        &values) const
{
  throw std::runtime_error("Ifpack2::LocalFilter does not implement getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void LocalFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow,
								    Teuchos::ArrayRCP<const LocalOrdinal> &indices,
								    Teuchos::ArrayRCP<const Scalar>       &values) const
{
  throw std::runtime_error("Ifpack2::LocalFilter does not implement getLocalRowView.");
}

}// namespace Ifpack2

#endif
