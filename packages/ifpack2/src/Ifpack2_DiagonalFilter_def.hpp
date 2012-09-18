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

#ifndef IFPACK2_DIAGONALFILTER_DEF_HPP
#define IFPACK2_DIAGONALFILTER_DEF_HPP
#include "Ifpack2_DiagonalFilter_decl.hpp"
#include <vector>

#include "Teuchos_Comm.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"


namespace Ifpack2 {

//==========================================================================
template<class MatrixType>
DiagonalFilter<MatrixType>::DiagonalFilter(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix,
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
	if(Teuchos::ScalarTraits<Scalar>::real(Values[i]) < 0)
	  mysign=-Teuchos::ScalarTraits<magnitudeType>::one();
	else
	  mysign=Teuchos::ScalarTraits<magnitudeType>::one();


	val_->replaceLocalValue(MyRow, Values[i] * (RelativeThreshold_ - 1) +
				AbsoluteThreshold_ * mysign);
	break;
      }
    }
  }
}

//=========================================================================
template<class MatrixType>
DiagonalFilter<MatrixType>::~DiagonalFilter() { }

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & DiagonalFilter<MatrixType>::getComm() const
{
  return A_->getComm();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP <typename MatrixType::node_type> DiagonalFilter<MatrixType>::getNode() const
{
  return A_->getNode();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
DiagonalFilter<MatrixType>::getRowMap() const
{
  return A_->getRowMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
DiagonalFilter<MatrixType>::getColMap() const
{
  return A_->getColMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
DiagonalFilter<MatrixType>::getDomainMap() const
{
  return A_->getDomainMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
DiagonalFilter<MatrixType>::getRangeMap() const
{
  return A_->getRangeMap();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >

DiagonalFilter<MatrixType>::getGraph() const
{
  return A_->getGraph();
}

//==========================================================================
template<class MatrixType>
global_size_t DiagonalFilter<MatrixType>::getGlobalNumRows() const
{
  return A_->getGlobalNumRows();
}

//==========================================================================
template<class MatrixType>
global_size_t DiagonalFilter<MatrixType>::getGlobalNumCols() const
{
  return A_->getGlobalNumCols();
}

//==========================================================================
template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getNodeNumRows() const
{
  return A_->getNodeNumRows();
}

//==========================================================================
 
template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getNodeNumCols() const
{
  return A_->getNodeNumCols();
}

//==========================================================================  
template<class MatrixType>
typename MatrixType::global_ordinal_type DiagonalFilter<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}

//========================================================================== 
template<class MatrixType>
global_size_t DiagonalFilter<MatrixType>::getGlobalNumEntries() const
{
  return A_->getGlobalNumEntries();
}

//========================================================================== 
template<class MatrixType>
size_t DiagonalFilter<MatrixType>::getNodeNumEntries() const
{
  return A_->getNodeNumEntries();
}

//==========================================================================
template<class MatrixType> 
size_t DiagonalFilter<MatrixType>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
{
  return A_->getNumEntriesInGlobalRow(globalRow);
}

//==========================================================================  
template<class MatrixType> 
size_t DiagonalFilter<MatrixType>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
{
  return A_->getNumEntriesInLocalRow(localRow);
}

//==========================================================================  
template<class MatrixType>   
global_size_t DiagonalFilter<MatrixType>::getGlobalNumDiags() const
{
  return A_->getGlobalNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t DiagonalFilter<MatrixType>::getNodeNumDiags() const
{
  return A_->getNodeNumDiags();
}

//==========================================================================  
template<class MatrixType>   
size_t DiagonalFilter<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return A_->getGlobalMaxNumRowEntries();
}

//==========================================================================  
template<class MatrixType>   
size_t DiagonalFilter<MatrixType>::getNodeMaxNumRowEntries() const
{
  return A_->getNodeMaxNumRowEntries();
}

//==========================================================================  
template<class MatrixType>   
bool DiagonalFilter<MatrixType>::hasColMap() const
{
  return A_->hasColMap();
}

//==========================================================================  
template<class MatrixType>   
bool DiagonalFilter<MatrixType>::isLowerTriangular() const
{
  return A_->isLowerTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool DiagonalFilter<MatrixType>::isUpperTriangular() const
{
  return A_->isUpperTriangular();
}

//==========================================================================  
template<class MatrixType>   
bool DiagonalFilter<MatrixType>::isLocallyIndexed() const
{
  return A_->isLocallyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool DiagonalFilter<MatrixType>::isGloballyIndexed() const
{
  return A_->isGloballyIndexed();
}

//==========================================================================  
template<class MatrixType>   
bool DiagonalFilter<MatrixType>::isFillComplete() const
{
  return A_->isFillComplete();
}
  
//==========================================================================
template<class MatrixType> 
void DiagonalFilter<MatrixType>::getGlobalRowCopy(GlobalOrdinal GlobalRow,
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

//==========================================================================  
template<class MatrixType> 
void DiagonalFilter<MatrixType>::getLocalRowCopy(LocalOrdinal LocalRow, 
						 const Teuchos::ArrayView<LocalOrdinal> &Indices, 
						 const Teuchos::ArrayView<Scalar> &Values,
						 size_t &NumEntries) const 
{ 
  Teuchos::ArrayRCP< const Scalar > myvals=val_->get1dView();

  A_->getLocalRowCopy(LocalRow, Indices,Values,NumEntries);

  if (pos_[LocalRow] != -1)
    Values[pos_[LocalRow]] += myvals[LocalRow];
}

//==========================================================================  
template<class MatrixType> 
void DiagonalFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
						  Teuchos::ArrayView<const GlobalOrdinal> &indices, 
						  Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::DiagonalFilter: does not support getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void DiagonalFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow, 
						 Teuchos::ArrayView<const LocalOrdinal> &indices, 
						 Teuchos::ArrayView<const Scalar> &values) const
{
  throw std::runtime_error("Ifpack2::DiagonalFilter: does not support getLocalRowView.");
}

//==========================================================================  
template<class MatrixType> 
void DiagonalFilter<MatrixType>::getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const
{
  // Returns the matrix's actual diagonal, rather than the "filtered" diagonal.
  // This is dubious, but it duplicates the functionality of Old Ifpack.
  return A_->getLocalDiagCopy(diag);
}

//========================================================================== 
template<class MatrixType> 
void DiagonalFilter<MatrixType>::leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::DiagonalFilter does not support leftScale.");
}

//==========================================================================  
template<class MatrixType> 
void DiagonalFilter<MatrixType>::rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) 
{
  throw std::runtime_error("Ifpack2::DiagonalFilter does not support rightScale.");
}

//==========================================================================  
template<class MatrixType> 
void DiagonalFilter<MatrixType>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
				       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
				       Teuchos::ETransp mode, 
				       Scalar alpha,
				       Scalar beta) const
{
  Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  A_->apply(X,Y,mode,alpha,beta);
  Y.elementWiseMultiply(one,*val_,X,one);
}
  

//==========================================================================  
template<class MatrixType> 
bool DiagonalFilter<MatrixType>::hasTransposeApply() const
{
  return A_->hasTransposeApply();
}

//==========================================================================
template<class MatrixType> 
bool DiagonalFilter<MatrixType>::supportsRowViews() const
{
  return false;
}

//==========================================================================  
template<class MatrixType> 
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType DiagonalFilter<MatrixType>::getFrobeniusNorm() const
{
  throw std::runtime_error("Ifpack2::DiagonalFilter does not implement getFrobeniusNorm.");
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void DiagonalFilter<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
								     Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
								     Teuchos::ArrayRCP<const Scalar>        &values) const
{
  throw std::runtime_error("Ifpack2::DiagonalFilter does not implement getGlobalRowView.");
}

//==========================================================================  
template<class MatrixType> 
TPETRA_DEPRECATED  void DiagonalFilter<MatrixType>::getLocalRowView(LocalOrdinal LocalRow,
								    Teuchos::ArrayRCP<const LocalOrdinal> &indices,
								    Teuchos::ArrayRCP<const Scalar>       &values) const
{
  throw std::runtime_error("Ifpack2::DiagonalFilter does not implement getLocalRowView.");
}

}// namespace Ifpack2

#endif
