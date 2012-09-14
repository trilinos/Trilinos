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

#ifndef IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP
#define IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP
#include "Ifpack2_OverlappingRowMatrix_decl.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Teuchos_CommHelpers.hpp"

namespace Ifpack2 {
//==========================================================================
// Standard constructor
template<class MatrixType>
OverlappingRowMatrix<MatrixType>::OverlappingRowMatrix(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix_in,
						       int OverlapLevel_in) :
  A_(Matrix_in),
  OverlapLevel_(OverlapLevel_in),
  UseSubComm_(false)
{  
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;

  // Short form names
  typedef typename Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>        MapType;
  typedef typename Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>           ImportType;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsMatrixType;
  typedef typename Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> RowMatrixType;

  TEUCHOS_TEST_FOR_EXCEPTION(OverlapLevel_ <= 0, std::runtime_error, "Ifpack2::OverlappingRowMatrix: OverlapLevel must be > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(A_->getComm()->getSize() == 1, std::runtime_error, "Ifpack2::OverlappingRowMatrix: Matrix must be parallel.");

  RCP<const CrsMatrixType> ACRS = Teuchos::rcp_dynamic_cast<const CrsMatrixType,const RowMatrixType>(A_);
  TEUCHOS_TEST_FOR_EXCEPTION(ACRS == Teuchos::null, std::runtime_error, "Ifpack2::OverlappingRowMatrix: Matrix must be Tpetra::CrsMatrix.");
  
  NumMyRowsA_ = A_->getNodeNumRows();

  GlobalOrdinal global_invalid = Teuchos::OrdinalTraits<global_size_t>::invalid();

  // Temp arrays
  Array<GlobalOrdinal> ExtElements;  
  RCP<MapType>         TmpMap;
  RCP<CrsMatrixType>   TmpMatrix; 
  RCP<ImportType>      TmpImporter;
  RCP<const MapType>   RowMap, ColMap;

  // The big import loop
  for (int overlap = 0 ; overlap < OverlapLevel_ ; ++overlap) {
    // Get the current maps
    if(overlap==0){
      RowMap = A_->getRowMap();
      ColMap = A_->getColMap(); 
    }
    else {
      RowMap = TmpMatrix->getRowMap();
      ColMap = TmpMatrix->getColMap(); 
    }

    size_t size = ColMap->getNodeNumElements() - RowMap->getNodeNumElements(); 
    Array<GlobalOrdinal> mylist(size); 
    size_t count = 0; 

    // define the set of rows that are in ColMap but not in RowMap 
    for (LocalOrdinal i = 0 ; (size_t) i < ColMap->getNodeNumElements() ; ++i) { 
      GlobalOrdinal GID = ColMap->getGlobalElement(i); 
      if (A_->getRowMap()->getLocalElement(GID) ==  global_invalid) { 
	typename Array<GlobalOrdinal>::iterator pos 
	  = find(ExtElements.begin(),ExtElements.end(),GID); 
        if (pos == ExtElements.end()) { 
          ExtElements.push_back(GID);
          mylist[count] = GID; 
          ++count; 
        } 
      } 
    }
    
    // Allocate & import new matrices, maps, etc.
    TmpMap      = rcp(new MapType(global_invalid,mylist(0,count-1),Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),A_->getComm(),A_->getNode()));
    TmpMatrix   = rcp(new CrsMatrixType(TmpMap,0));
    TmpImporter = rcp(new ImportType(A_->getRowMap(),TmpMap));

    TmpMatrix->doImport(*ACRS,*TmpImporter,Tpetra::INSERT);
    TmpMatrix->fillComplete(A_->getDomainMap(),TmpMap);
  }

  // build the map containing all the nodes (original
  // matrix + extended matrix)
  Array<GlobalOrdinal> mylist(NumMyRowsA_ + ExtElements.size());
  for (LocalOrdinal i = 0 ; (size_t)i < NumMyRowsA_ ; ++i)
    mylist[i] = A_->getRowMap()->getGlobalElement(i);
  for (LocalOrdinal i = 0 ; i < ExtElements.size() ; ++i)
    mylist[i + NumMyRowsA_] = ExtElements[i];

  RowMap_= rcp(new MapType(global_invalid,mylist(),Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),A_->getComm(),A_->getNode()));
  ColMap_= RowMap_;

  // now build the map corresponding to all the external nodes
  // (with respect to A().RowMatrixRowMap().
  ExtMap_      = rcp(new MapType(global_invalid,ExtElements(),Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),A_->getComm(),A_->getNode()));
  ExtMatrix_   = rcp(new CrsMatrixType(ExtMap_,ColMap_,0));
  ExtImporter_ = rcp(new ImportType(A_->getRowMap(),ExtMap_));

  RCP<CrsMatrixType> ExtMatrixCRS = Teuchos::rcp_dynamic_cast<CrsMatrixType,RowMatrixType>(ExtMatrix_);
  ExtMatrixCRS->doImport(*ACRS,*ExtImporter_,Tpetra::INSERT);
  ExtMatrixCRS->fillComplete(A_->getDomainMap(),RowMap_);

  Importer_ = rcp(new ImportType(A_->getRowMap(),RowMap_));

  // fix indices for overlapping matrix
  NumMyRowsB_ = ExtMatrix_->getNodeNumRows();
  NumMyRows_ = NumMyRowsA_ + NumMyRowsB_;
  NumMyCols_ = NumMyRows_;

  NumMyDiagonals_ = A_->getNodeNumDiags() + ExtMatrix_->getNodeNumDiags(); 
  NumMyNonzeros_  = A_->getNodeNumEntries() + ExtMatrix_->getNodeNumEntries();

  // FIXME: Fix this later when Teuchos::Comm gets redone
  Tpetra::global_size_t NumMyNonzeros_tmp = NumMyNonzeros_;
  Teuchos::reduceAll<int,Tpetra::global_size_t>(*A_->getComm(),Teuchos::REDUCE_SUM,NumMyNonzeros_tmp,Teuchos::outArg(NumGlobalNonzeros_));  
  Tpetra::global_size_t NumMyRows_tmp = NumMyRows_;
  Teuchos::reduceAll<int,Tpetra::global_size_t>(*A_->getComm(),Teuchos::REDUCE_SUM,NumMyRows_tmp,Teuchos::outArg(NumGlobalRows_));  

  MaxNumEntries_ = A_->getNodeMaxNumRowEntries();  
  if (MaxNumEntries_ < ExtMatrix_->getNodeMaxNumRowEntries())
    MaxNumEntries_ = ExtMatrix_->getNodeMaxNumRowEntries();

  // Resize temp arrays
  Indices_.resize(MaxNumEntries_);
  Values_.resize(MaxNumEntries_);
}



//==========================================================================
// Sub-communicator constructor
template<class MatrixType>
OverlappingRowMatrix<MatrixType>::OverlappingRowMatrix(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix_in,
							      int OverlapLevel_in, int subdomainID)
{
  //FIXME
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix - subdomain code not implemented yet.");
}

//==========================================================================
// Destructor
template<class MatrixType>
OverlappingRowMatrix<MatrixType>::~OverlappingRowMatrix()
{

}
//==========================================================================
// Returns the communicator.
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & OverlappingRowMatrix<MatrixType>::getComm() const
{
  return A_->getComm();
}
  
//==========================================================================
// Returns the underlying node.
template<class MatrixType>
Teuchos::RCP<typename MatrixType::node_type> OverlappingRowMatrix<MatrixType>::getNode() const
{
  return A_->getNode();
}
  
//==========================================================================
// Returns the Map that describes the row distribution in this matrix.
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > & OverlappingRowMatrix<MatrixType>::getRowMap() const
{
  return RowMap_;
}
  
//==========================================================================
//  Returns the Map that describes the column distribution in this matrix.
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > & OverlappingRowMatrix<MatrixType>::getColMap() const
{
  return ColMap_;
}

//==========================================================================
//  Returns the Map that describes the column distribution in this matrix.
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > & OverlappingRowMatrix<MatrixType>::getDomainMap() const
{
  return A_->getDomainMap();
}

//==========================================================================
//  Returns the Map that describes the column distribution in this matrix.
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > & OverlappingRowMatrix<MatrixType>::getRangeMap() const
{
  return A_->getRangeMap();
}
  
//==========================================================================
  // Returns the RowGraph associated with this matrix. 
template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > OverlappingRowMatrix<MatrixType>::getGraph() const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix - ERROR getGraph() is not implemented.");
}
  
//==========================================================================
// Returns the number of global rows in this matrix.
template<class MatrixType>
global_size_t OverlappingRowMatrix<MatrixType>::getGlobalNumRows() const
{
  return NumGlobalRows_;
}
  
//==========================================================================
//  Returns the number of global columns in this matrix.
template<class MatrixType>
global_size_t OverlappingRowMatrix<MatrixType>::getGlobalNumCols() const
{
  return NumGlobalRows_;
}
  
//==========================================================================
// Returns the number of rows owned on the calling node.
template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNodeNumRows() const
{
  return NumMyRows_;
}
  
//==========================================================================
// Returns the number of columns needed to apply the forward operator on this node, i.e., the number of elements listed in the column map.
template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNodeNumCols() const
{
  return NumMyCols_;
}
  
//==========================================================================
// Returns the index base for global indices for this matrix. 
template<class MatrixType>
typename MatrixType::global_ordinal_type OverlappingRowMatrix<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}
  
//==========================================================================
// Returns the global number of entries in this matrix.
template<class MatrixType>
Tpetra::global_size_t OverlappingRowMatrix<MatrixType>::getGlobalNumEntries() const
{
  return NumGlobalNonzeros_;
}
  
//==========================================================================
// Returns the local number of entries in this matrix.
template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNodeNumEntries() const
{
  return NumMyNonzeros_;
}
  
//==========================================================================
//  Returns the current number of entries on this node in the specified global row.
/* Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row is not valid for this graph. */
template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
{
  LocalOrdinal localRow = RowMap_->getLocalElement(globalRow);
  if(localRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) return Teuchos::OrdinalTraits<size_t>::invalid();
  else return getNumEntriesInLocalRow(localRow);
}

  
//==========================================================================
// Returns the current number of entries on this node in the specified local row.
/* Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
{
  if ((size_t)localRow < NumMyRowsA_)
    return A_->getNumEntriesInLocalRow(localRow);
  else
    return ExtMatrix_->getNumEntriesInLocalRow((LocalOrdinal)(localRow-NumMyRowsA_));
}
  
//==========================================================================
//  Returns the number of global diagonal entries, based on global row/column index comparisons. 
template<class MatrixType>
global_size_t OverlappingRowMatrix<MatrixType>::getGlobalNumDiags() const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix::getGlobalNumDiags() not supported.");
}
  
//==========================================================================
//  Returns the number of local diagonal entries, based on global row/column index comparisons. 
template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNodeNumDiags() const
{
  return A_->getNodeNumDiags();
}
  
//==========================================================================
//  Returns the maximum number of entries across all rows/columns on all nodes.
template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getGlobalMaxNumRowEntries() const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix::getGlobalMaxNumRowEntries() not supported.");
}
  
//==========================================================================
template<class MatrixType>
//  Returns the maximum number of entries across all rows/columns on this node.
size_t OverlappingRowMatrix<MatrixType>::getNodeMaxNumRowEntries() const
{
  return MaxNumEntries_;
}
  
//==========================================================================
//  Indicates whether this matrix has a well-defined column map
template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::hasColMap() const
{
  return true;
}
  
//==========================================================================
//  Indicates whether this matrix is lower triangular.
template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isLowerTriangular() const
{
  return A_->isLowerTriangular();
}
  
//==========================================================================
//  Indicates whether this matrix is upper triangular.
template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isUpperTriangular() const
{
  return A_->isUpperTriangular();
} 
//==========================================================================
//  If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isLocallyIndexed() const
{
  return true;
}
   
//==========================================================================
//  If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isGloballyIndexed() const
{
  return false;
}
  
//==========================================================================
// Returns \c true if fillComplete() has been called.
template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isFillComplete() const
{
  return true;
}
  
//==========================================================================
// Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
  /*
    \param LocalRow - (In) Global row number for which indices are desired.
    \param Indices - (Out) Global column indices corresponding to values.
    \param Values - (Out) Matrix values.
    \param NumEntries - (Out) Number of indices.
    
    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
    returned as Teuchos::OrdinalTraits<size_t>::invalid().
  */
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::getGlobalRowCopy(GlobalOrdinal GlobalRow,
				const Teuchos::ArrayView<GlobalOrdinal> &Indices,
				const Teuchos::ArrayView<Scalar> &Values,
				size_t &NumEntries) const
{
  LocalOrdinal LocalRow=RowMap_->getLocalElement(GlobalRow);
  if(LocalRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) NumEntries = Teuchos::OrdinalTraits<size_t>::invalid();
  else {
    if ((size_t)LocalRow < NumMyRowsA_) 
      A_->getGlobalRowCopy(GlobalRow,Indices,Values,NumEntries);
    else 
      ExtMatrix_->getGlobalRowCopy(GlobalRow,Indices,Values,NumEntries);
  }

}
  
//==========================================================================
// Extract a list of entries in a specified local row of the graph. Put into storage allocated by calling routine.
  /*
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices - (Out) Local column indices corresponding to values.
    \param Values - (Out) Matrix values.
    \param NumIndices - (Out) Number of indices.
    
    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
    returned as Teuchos::OrdinalTraits<size_t>::invalid().
  */
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::getLocalRowCopy(LocalOrdinal LocalRow, 
						       const Teuchos::ArrayView<LocalOrdinal> &Indices, 
						       const Teuchos::ArrayView<Scalar> &Values,
						       size_t &NumEntries) const
{
  if ((size_t)LocalRow < NumMyRowsA_) 
    A_->getLocalRowCopy(LocalRow,Indices,Values,NumEntries);
  else 
    ExtMatrix_->getLocalRowCopy(LocalRow-(LocalOrdinal)(NumMyRowsA_),Indices,Values,NumEntries);
}
  
//==========================================================================
// Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*
    \param GlobalRow - (In) Global row number for which indices are desired.
    \param Indices   - (Out) Global column indices corresponding to values.
    \param Values    - (Out) Row values
    \pre <tt>isLocallyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>
    
    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
  */
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::getGlobalRowView(GlobalOrdinal GlobalRow, 
							Teuchos::ArrayView<const GlobalOrdinal> &indices, 
							Teuchos::ArrayView<const Scalar> &values) const
{
  LocalOrdinal LocalRow=RowMap_->getLocalElement(GlobalRow);
  if(LocalRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid())  {
    indices=Teuchos::null; 
    values =Teuchos::null;
  }
  else {
    if ((size_t)LocalRow < NumMyRowsA_) 
      A_->getGlobalRowView(GlobalRow,indices,values);
    else 
      ExtMatrix_->getGlobalRowView(GlobalRow,indices,values);
  }
}
  
//==========================================================================
// Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices  - (Out) Global column indices corresponding to values.
    \param Values   - (Out) Row values
    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>
    
    Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
  */
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::getLocalRowView(LocalOrdinal LocalRow, 
						       Teuchos::ArrayView<const LocalOrdinal> &indices, 
						       Teuchos::ArrayView<const Scalar> &values) const
{
  if ((size_t)LocalRow < NumMyRowsA_) 
    A_->getLocalRowView(LocalRow,indices,values);
  else 
    ExtMatrix_->getLocalRowView(LocalRow-(LocalOrdinal)(NumMyRowsA_),indices,values);
}
  
//==========================================================================
//  Get a copy of the diagonal entries owned by this node, with local row indices.
  /* Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
    the zero and non-zero diagonals owned by this node. */
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix::getLocalDiagCopy not supported.");
}
  
//==========================================================================
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x)
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix does not support leftScale.");
}
  
//==========================================================================
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x)
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix does not support leftScale.");
}
  
//==========================================================================
// Returns the Frobenius norm of the matrix. 
template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType OverlappingRowMatrix<MatrixType>::getFrobeniusNorm() const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix does not support getFrobeniusNorm.");
}



//==========================================================================
  // \brief Computes the operator-multivector application.
  /* Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
    vary according to the values of \c alpha and \c beta. Specifically
    - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
    - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.

    This is analagous to the *Multiply* function in Ifpack, not the *Apply*
  */
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
					     Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
					     Teuchos::ETransp mode, 
					     Scalar alpha,
					     Scalar beta) const
{
  // Note: This isn't AztecOO compliant.  But neither was Ifpack's version.
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
			     "Ifpack2::OverlappingRowMatrix::apply ERROR: X.getNumVectors() != Y.getNumVectors().");

  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr = X.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       y_ptr = Y.get2dViewNonConst();
  Y.putScalar(zero);
  size_t NumVectors = Y.getNumVectors();

  for (size_t i = 0 ; i < NumMyRowsA_ ; ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    A_->getLocalRowCopy(i,Indices_(),Values_(),Nnz);
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

 for (size_t i = 0 ; i < NumMyRowsB_ ; ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    ExtMatrix_->getLocalRowCopy(i,Indices_(),Values_(),Nnz);
    if (mode==Teuchos::NO_TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][NumMyRowsA_+i] += Values_[j] * x_ptr[k][Indices_[j]];      
    }
    else if (mode==Teuchos::TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][NumMyRowsA_+Indices_[j]] += Values_[j] * x_ptr[k][i];
    }
    else { //mode==Teuchos::CONJ_TRANS
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][NumMyRowsA_+Indices_[j]] += Teuchos::ScalarTraits<Scalar>::conjugate(Values_[j]) * x_ptr[k][i];
    }
  }
}
  


// ======================================================================
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::importMultiVector(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
							 Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &OvX, 
							 Tpetra::CombineMode CM)
{
  OvX.doImport(X,*Importer_,CM);
}


// ======================================================================
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::exportMultiVector(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &OvX, 
							 Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
							 Tpetra::CombineMode CM)
{
  X.doExport(OvX,*Importer_,CM);
}

//==========================================================================
// Indicates whether this operator supports applying the adjoint operator.
template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::hasTransposeApply() const
{
  return true;
}


} // namespace Ifpack2

#endif // IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP





















#ifdef OLD_AND_BUSTED
// ====================================================================== 

#ifdef IFPACK_SUBCOMM_CODE

// ====================================================================== 
// Constructor for the case of two or more cores per subdomain
Ifpack_OverlappingRowMatrix::
Ifpack_OverlappingRowMatrix(const RCP<const Epetra_RowMatrix>& Matrix_in,
                            int OverlapLevel_in, int subdomainID)  :
  Matrix_(Matrix_in),
  OverlapLevel_(OverlapLevel_in)
{

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  double t0 = MPI_Wtime();
  double t1, true_t0=t0;
#endif
  //FIXME timing

  const int votesMax = INT_MAX;

  // should not be here if no overlap
  if (OverlapLevel_in == 0)
    IFPACK_CHK_ERRV(-1);

  // nothing to do as well with one process
  if (Comm().NumProc() == 1)
    IFPACK_CHK_ERRV(-1);

  // subdomainID is the node (aka socket) number, and so is system dependent
  // nodeComm is the communicator for all the processes on a particular node
  // these processes will have the same subdomainID.
# ifdef HAVE_MPI
  const Epetra_MpiComm *pComm = dynamic_cast<const Epetra_MpiComm*>( &Comm() );
  assert(pComm != NULL);
  MPI_Comm nodeMPIComm;
  MPI_Comm_split(pComm->Comm(),subdomainID,pComm->MyPID(),&nodeMPIComm);
  Epetra_MpiComm *nodeComm = new Epetra_MpiComm(nodeMPIComm);
# else
  Epetra_SerialComm *nodeComm =  dynamic_cast<Epetra_MpiComm*>( &(Matrix_in->RowMatrixRowMap().Comm()) );
# endif
  
  NumMyRowsA_ = A().NumMyRows();

  RCP<Epetra_Map> TmpMap;
  RCP<Epetra_CrsMatrix> TmpMatrix; 
  RCP<Epetra_Import> TmpImporter;

  // importing rows corresponding to elements that are 
  // in ColMap, but not in RowMap 
  const Epetra_Map *RowMap; 
  const Epetra_Map *ColMap; 
  const Epetra_Map *DomainMap;

  // TODO Count #connections from nodes I own to each ghost node


  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime(); 
  fprintf(stderr,"[node %d]: time for initialization %2.3e\n", subdomainID, t1-t0);
  t0=t1;
#endif
  //FIXME timing

#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: overlap hash table allocs %2.3e\n", subdomainID, t1-t0);
  t0=t1;
#endif

  Teuchos::Hashtable<int,int> colMapTable(3 * A().RowMatrixColMap().NumMyElements() );

  // ghostTable holds off-node GIDs that are connected to on-node rows and can potentially be this PID's overlap
  // TODO hopefully 3 times the # column entries is a reasonable table size
  Teuchos::Hashtable<int,int> ghostTable(3 * A().RowMatrixColMap().NumMyElements() );

  /* ** ************************************************************************** ** */
  /* ** ********************** start of main overlap loop ************************ ** */
  /* ** ************************************************************************** ** */
  for (int overlap = 0 ; overlap < OverlapLevel_in ; ++overlap)
  {
    // nbTable holds node buddy GIDs that are connected to current PID's rows, i.e., GIDs that should be in the overlapped
    // matrix's column map

    if (TmpMatrix != Teuchos::null) {
      RowMap = &(TmpMatrix->RowMatrixRowMap()); 
      ColMap = &(TmpMatrix->RowMatrixColMap()); 
      DomainMap = &(TmpMatrix->OperatorDomainMap());
    }
    else {
      RowMap = &(A().RowMatrixRowMap()); 
      ColMap = &(A().RowMatrixColMap()); 
      DomainMap = &(A().OperatorDomainMap());
    }

#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap pointer pulls %2.3e\n", subdomainID, t1-t0);
    t0=t1;
#endif
    
    // For each column ID, determine the owning node (as opposed to core)
    // ID of the corresponding row.
    Epetra_IntVector colIdList( *ColMap );
    Epetra_IntVector rowIdList(*DomainMap);
    rowIdList.PutValue(subdomainID);  
    Teuchos::RCP<Epetra_Import> nodeIdImporter = rcp(new Epetra_Import( *ColMap, *DomainMap ));

#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap intvector/importer allocs %2.3e\n", subdomainID, t1-t0);
    t0=t1;
#endif
    
    colIdList.Import(rowIdList,*nodeIdImporter,Insert);
    
    int size = ColMap->NumMyElements() - RowMap->NumMyElements(); 
    int count = 0; 

#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap col/row ID imports %2.3e\n", subdomainID, t1-t0);
    t0=t1;
#endif

    
    // define the set of off-node rows that are in ColMap but not in RowMap
    // This naturally does not include off-core rows that are on the same node as me, i.e., node buddy rows.
    for (int i = 0 ; i < ColMap->NumMyElements() ; ++i) {
      int GID = ColMap->GID(i); 
      int myvotes=0;
      if (ghostTable.containsKey(GID)) myvotes = ghostTable.get(GID);
      if ( colIdList[i] != subdomainID && myvotes < votesMax) // don't include anybody found in a previous round
      {
        int votes;
        if (ghostTable.containsKey(GID)) {
          votes = ghostTable.get(GID);
          votes++;
          ghostTable.put(GID,votes);
        } else {
          ghostTable.put(GID,1);
        }
      }
    } //for (int i = 0 ; i < ColMap->NumMyElements() ; ++i)

    Teuchos::Array<int> gidsarray,votesarray;
    ghostTable.arrayify(gidsarray,votesarray);
    int *gids = gidsarray.getRawPtr();
    int *votes = votesarray.getRawPtr();

    /*
       This next bit of code decides which node buddy (NB) gets which ghost points.  Everyone sends their
       list of ghost points to pid 0 of the local subcommunicator.  Pid 0 decides who gets what:

          - if a ghost point is touched by only one NB, that NB gets the ghost point
          - if two or more NBs share a ghost point, whichever NB has the most connections to the ghost
            point gets it.
    */

#   ifdef HAVE_MPI  //FIXME What if we build in serial?!  This file will likely not compile.
    int *cullList;
    int ncull;
    //int mypid = nodeComm->MyPID();

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap before culling %2.3e\n", subdomainID, t1-t0);
    t0=t1;
#endif
    //FIXME timing


    if (nodeComm->MyPID() == 0)
    {
      // Figure out how much pid 0 is to receive
      MPI_Status status;
      int *lengths = new int[nodeComm->NumProc()+1];
      lengths[0] = 0;
      lengths[1] = ghostTable.size();
      for (int i=1; i<nodeComm->NumProc(); i++) {
        int leng;
        MPI_Recv( &leng, 1, MPI_INT, i, MPI_ANY_TAG, nodeComm->Comm(), &status);
        lengths[i+1] = lengths[i] + leng;
      }
      int total = lengths[nodeComm->NumProc()];

      int* ghosts = new int[total];
      for (int i=0; i<total; i++) ghosts[i] = -9;
      int *round  = new int[total];
      int *owningpid  = new int[total];

      for (int i=1; i<nodeComm->NumProc(); i++) {
        int count = lengths[i+1] - lengths[i];
        MPI_Recv( ghosts+lengths[i], count, MPI_INT, i, MPI_ANY_TAG, nodeComm->Comm(), &status);
        MPI_Recv( round+lengths[i], count, MPI_INT, i, MPI_ANY_TAG, nodeComm->Comm(), &status);
      }

      // put in pid 0's info
      for (int i=0; i<lengths[1]; i++) {
        ghosts[i] = gids[i];
        round[i] = votes[i];
        owningpid[i] = 0;
      }

      // put in the pid associated with each ghost
      for (int j=1; j<nodeComm->NumProc(); j++) {
        for (int i=lengths[j]; i<lengths[j+1]; i++) {
          owningpid[i] = j;
        }
      }

      // sort everything based on the ghost gids
      int* roundpid[2];
      roundpid[0] = round; roundpid[1] = owningpid;
      Epetra_Util epetraUtil;
      epetraUtil.Sort(true,total,ghosts,0,0,2,roundpid);

      //set up arrays that get sent back to node buddies and that tell them what ghosts to cull
      int *nlosers = new int[nodeComm->NumProc()];
      int** losers = new int*[nodeComm->NumProc()];
      for (int i=0; i<nodeComm->NumProc(); i++) {
        nlosers[i] = 0;
        losers[i] = new int[ lengths[i+1]-lengths[i] ];
      }

      // Walk through ghosts array and and for each sequence of duplicate ghost GIDs, choose just one NB to keep it.
      // The logic is pretty simple.  The ghost list is sorted, so all duplicate PIDs are together.
      // The list is traversed.  As duplicates are found, node pid 0 keeps track of the current "winning"
      // pid.  When a pid is determined to have "lost" (less votes/connections to the current GID), the
      // GID is added to that pid's list of GIDs to be culled.  At the end of the repeated sequence, we have
      // a winner, and other NBs know whether they need to delete it from their import list.
      int max=0;   //for duplicated ghosts, index of pid with most votes and who hence keeps the ghost.
                   // TODO to break ties randomly

      for (int i=1; i<total; i++) {
        if (ghosts[i] == ghosts[i-1]) {
          int idx = i; // pid associated with idx is current "loser"
          if (round[i] > round[max]) {
            idx = max;
            max=i;
          }
          int j = owningpid[idx];
          losers[j][nlosers[j]++] = ghosts[idx];
        } else {
          max=i;
        }
      } //for (int i=1; i<total; i++)

      delete [] round;
      delete [] ghosts;
      delete [] owningpid;

      // send the arrays of ghost GIDs to be culled back to the respective node buddies
      for (int i=1; i<nodeComm->NumProc(); i++) {
        MPI_Send( nlosers+i, 1, MPI_INT, i, 8675, nodeComm->Comm());
        MPI_Send( losers[i], nlosers[i], MPI_INT, i, 8675, nodeComm->Comm());
      }

      //FIXME Unnecessary memory allocation and copying, but makes culling code cleaner
      //Could we stick this info into gids and votes, since neither is being used anymore?
      //TODO Instead of using "losers" arrays, just use "cullList" as in the else clause
      ncull = nlosers[0];
      cullList = new int[ncull+1];
      for (int i=0; i<nlosers[0]; i++)
        cullList[i] = losers[0][i];

      for (int i=0; i<nodeComm->NumProc(); i++)
        delete [] losers[i];

      delete [] lengths;
      delete [] losers;
      delete [] nlosers;

    } else { //everyone but pid 0

      // send to node pid 0 all ghosts that this pid could potentially import
      int hashsize = ghostTable.size();
      MPI_Send( &hashsize, 1, MPI_INT, 0, 8675, nodeComm->Comm());
      MPI_Send( gids, hashsize, MPI_INT, 0, 8675, nodeComm->Comm());
      MPI_Send( votes, hashsize, MPI_INT, 0, 8675, nodeComm->Comm());

      // receive the ghost GIDs that should not be imported (subset of the list sent off just above)
      MPI_Status status;
      MPI_Recv( &ncull, 1, MPI_INT, 0, 8675, nodeComm->Comm(), &status);
      cullList = new int[ncull+1];
      MPI_Recv( cullList, ncull, MPI_INT, 0, 8675, nodeComm->Comm(), &status);
    }


    //TODO clean out hash table after each time through overlap loop   4/1/07 JJH done moved both hash tables to local scope

    // Remove from my row map all off-node ghosts that will be imported by a node buddy.
    for (int i=0; i<ncull; i++) {
      try{ghostTable.remove(cullList[i]);}

      catch(...) {
        printf("pid %d: In OverlappingRowMatrix ctr, problem removing ghost elt %d from ghostTable\n",
               Comm().MyPID(),cullList[i]);
        fflush(stdout);
      }
    }//for

    delete [] cullList;

    // Save off the remaining ghost GIDs from the current overlap round.
    // These are off-node GIDs (rows) that I will import.
    gidsarray.clear(); votesarray.clear();
    ghostTable.arrayify(gidsarray,votesarray);

    vector<int> list(size); 
    count=0;
    for (int i=0; i<ghostTable.size(); i++) {
      // if votesarray[i] == votesmax, then this GID was found during a previous overlap round
      if (votesarray[i] < votesMax) {
        list[count] = gidsarray[i];
        ghostTable.put(gidsarray[i],votesMax); //set this GID's vote to the maximum so that this PID alway wins in any subsequent overlap tiebreaking.
        count++;
      }
    }

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap duplicate removal %2.3e\n", subdomainID, t1-t0);
    t0=t1;
#endif
    //FIXME timing

#   endif //ifdef HAVE_MPI

    TmpMap = rcp( new Epetra_Map(-1,count, &list[0],0,Comm()) );

    TmpMatrix = rcp( new Epetra_CrsMatrix(Copy,*TmpMap,0) ); 

    TmpImporter = rcp( new Epetra_Import(*TmpMap,A().RowMatrixRowMap()) ); 

    TmpMatrix->Import(A(),*TmpImporter,Insert); 
    TmpMatrix->FillComplete(A().OperatorDomainMap(),*TmpMap); 

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap TmpMatrix fillcomplete %2.3e\n", subdomainID, t1-t0);
    t0=t1;
#endif
    //FIXME timing


    // These next two imports get the GIDs that need to go into the column map of the overlapped matrix.

    // For each column ID in the overlap, determine the owning node (as opposed to core)
    // ID of the corresponding row.  Save those column IDs whose owning node is the current one.
    // This should get all the imported ghost GIDs.
    Epetra_IntVector ov_colIdList( TmpMatrix->ColMap() );
    ov_colIdList.PutValue(-1);
    Epetra_IntVector ov_rowIdList( TmpMatrix->RowMap() );
    ov_rowIdList.PutValue(subdomainID);  
    Teuchos::RCP<Epetra_Import> ov_nodeIdImporter = rcp(new Epetra_Import( TmpMatrix->ColMap(), TmpMatrix->RowMap()));
    ov_colIdList.Import(ov_rowIdList,*ov_nodeIdImporter,Insert);

    for (int i=0 ; i<ov_colIdList.MyLength(); i++) {
       if (ov_colIdList[i] == subdomainID) {
         int GID = (ov_colIdList.Map()).GID(i);
         colMapTable.put(GID,1);
       }
    }

    // Do a second import of the owning node ID from A's rowmap to TmpMat's column map.  This ensures that
    // all GIDs that belong to a node buddy and are in a ghost row's sparsity pattern will be in the final
    // overlapped matrix's column map.
    ov_colIdList.PutValue(-1);
    Epetra_IntVector ArowIdList( A().RowMatrixRowMap() );
    ArowIdList.PutValue(subdomainID);
    nodeIdImporter = rcp(new Epetra_Import( TmpMatrix->ColMap(), A().RowMatrixRowMap() ));
    ov_colIdList.Import(ArowIdList,*nodeIdImporter,Insert);

    for (int i=0 ; i<ov_colIdList.MyLength(); i++) {
       if (ov_colIdList[i] == subdomainID) {
         int GID = (ov_colIdList.Map()).GID(i);
         colMapTable.put(GID,1);
       }
    }

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap 2 imports to fix up colmap %2.3e\n", subdomainID, t1-t0);
    t0=t1;
#endif
    //FIXME timing

  } //for (int overlap = 0 ; overlap < OverlapLevel_in ; ++overlap)

  /* ** ************************************************************************ ** */
  /* ** ********************** end of main overlap loop ************************ ** */
  /* ** ************************************************************************ ** */

  // off-node GIDs that will be in the overlap
  vector<int> ghostElements; 

  Teuchos::Array<int> gidsarray,votesarray;
  ghostTable.arrayify(gidsarray,votesarray);
  for (int i=0; i<ghostTable.size(); i++) {
    ghostElements.push_back(gidsarray[i]);
  }

    for (int i = 0 ; i < A().RowMatrixColMap().NumMyElements() ; ++i) {
      int GID = A().RowMatrixColMap().GID(i);
      // Remove any entries that are in A's original column map
      if (colMapTable.containsKey(GID)) {
        try{colMapTable.remove(GID);}
        catch(...) {
          printf("pid %d: In OverlappingRowMatrix ctr, problem removing entry %d from colMapTable\n", Comm().MyPID(),GID);
          fflush(stdout);
        }
      }
    }

    // GIDs that will be in the overlapped matrix's column map
    vector<int> colMapElements; 

  gidsarray.clear(); votesarray.clear();
  colMapTable.arrayify(gidsarray,votesarray);
  for (int i=0; i<colMapTable.size(); i++)
    colMapElements.push_back(gidsarray[i]);

/*
   We need two maps here.  The first is the row map, which we've got by using my original rows
   plus whatever I've picked up in ghostElements.

   The second is a column map.  This map should include my rows, plus any columns that could come from node buddies.
   These GIDs come from the std:array colMapElements, which in turn comes from colMapTable.
   This map should *not* omit ghosts that have been imported by my node buddies, i.e., for any row that I own,
   the stencil should include all column GIDs (including imported ghosts) that are on the node.
*/

  // build the row map containing all the nodes (original matrix + off-node matrix)
  vector<int> rowList(NumMyRowsA_ + ghostElements.size());
  for (int i = 0 ; i < NumMyRowsA_ ; ++i)
    rowList[i] = A().RowMatrixRowMap().GID(i);
  for (int i = 0 ; i < (int)ghostElements.size() ; ++i)
    rowList[i + NumMyRowsA_] = ghostElements[i];

  // row map for the overlapped matrix (local + overlap)
  Map_ = rcp( new Epetra_Map(-1, NumMyRowsA_ + ghostElements.size(), &rowList[0], 0, Comm()) );

  // build the column map for the overlapping matrix
  //vector<int> colList(colMapElements.size());
  // column map for the overlapped matrix (local + overlap)
  //colMap_ = rcp( new Epetra_Map(-1, colMapElements.size(), &colList[0], 0, Comm()) );
  //for (int i = 0 ; i < (int)colMapElements.size() ; i++)
  //  colList[i] = colMapElements[i];
  vector<int> colList(A().RowMatrixColMap().NumMyElements() + colMapElements.size());
  int nc = A().RowMatrixColMap().NumMyElements();
  for (int i = 0 ; i < nc; i++)
    colList[i] = A().RowMatrixColMap().GID(i);
  for (int i = 0 ; i < (int)colMapElements.size() ; i++)
    colList[nc+i] = colMapElements[i];

  // column map for the overlapped matrix (local + overlap)
  //colMap_ = rcp( new Epetra_Map(-1, A().RowMatrixColMap().NumMyElements() + colMapElements.size(), &colList[0], 0, Comm()) );
  colMap_ = new Epetra_Map(-1, A().RowMatrixColMap().NumMyElements() + colMapElements.size(), &colList[0], 0, Comm());
    

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time to init B() col/row maps %2.3e\n", subdomainID, t1-t0);
  t0=t1;
#endif
  //FIXME timing

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++ start of sort
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // build the column map, but don't use a copy constructor, b/c local communicator SubComm_ is
  // different from that of Matrix.
  try {
    // build row map based on the local communicator.  We need this temporarily to build the column map.
    Epetra_Map* nodeMap = new Epetra_Map(-1,NumMyRowsA_ + ghostElements.size(),&rowList[0],0,*nodeComm);
    int numMyElts = colMap_->NumMyElements();
    assert(numMyElts!=0);

    // The column map *must* be sorted: first locals, then ghosts.  The ghosts
    // must be further sorted so that they are contiguous by owning processor.

    // For each GID in column map, determine owning PID in local communicator.
    int* myGlobalElts = new int[numMyElts];
    colMap_->MyGlobalElements(myGlobalElts);
    int* pidList = new int[numMyElts];
    nodeMap->RemoteIDList(numMyElts, myGlobalElts, pidList, 0);

    // First sort column map GIDs based on the owning PID in local communicator.
    Epetra_Util Util;
    int *tt[1];
    tt[0] = myGlobalElts;
    Util.Sort(true, numMyElts, pidList, 0, (double**)0, 1, tt);

    // For each remote PID, sort the column map GIDs in ascending order.
    // Don't sort the local PID's entries.
    int localStart=0;
    while (localStart<numMyElts) {
      int currPID = (pidList)[localStart];
      int i=localStart;
      while (i<numMyElts && (pidList)[i] == currPID) i++;
      if (currPID != nodeComm->MyPID())
        Util.Sort(true, i-localStart, (myGlobalElts)+localStart, 0, 0, 0, 0);
      localStart = i;
    }

    // now move the local entries to the front of the list
    localStart=0;
    while (localStart<numMyElts && (pidList)[localStart] != nodeComm->MyPID()) localStart++;
    assert(localStart != numMyElts);
    int localEnd=localStart;
    while (localEnd<numMyElts && (pidList)[localEnd] == nodeComm->MyPID()) localEnd++;
    int* mySortedGlobalElts = new int[numMyElts];
    int leng = localEnd - localStart;
    /* This is a little gotcha.  It appears that the ordering of the column map's local entries
       must be the same as that of the domain map's local entries.  See the comment in method
       MakeColMap() in Epetra_CrsGraph.cpp, line 1072. */
    int *rowGlobalElts =  nodeMap->MyGlobalElements();

    /*TODO TODO TODO NTS rows 68 and 83 show up as local GIDs in rowGlobalElts for both pids 0 & 1 on node 0.
                    This seems to imply that there is something wrong with rowList!!
                    It is ok for the overlapped matrix's row map to have repeated entries (it's overlapped, after all).
                    But ... on a node, there should be no repeats.
                    Ok, here's the underlying cause.  On node 0, gpid 1 finds 83 on overlap round 0.  gpid 0 finds 83
                    on overlap round 1.  This shouldn't happen .. once a node buddy "owns" a row, no other node buddy
                    should ever import ("own") that row.

                    Here's a possible fix.  Extend tie breaking across overlap rounds.  This means sending to lpid 0
                    a list of *all* rows you've imported (this round plus previous rounds) for tie-breaking
                    If a nb won a ghost row in a previous round, his votes for that ghost row should be really high, i.e.,
                    he should probably win that ghost row forever.
                    7/17/09

      7/28/09 JJH I believe I've fixed this, but I thought it might be helpful to have this comment in here, in case I missed something.
    */

    //move locals to front of list
    for (int i=0; i<leng; i++) {
      /* //original code */
      (mySortedGlobalElts)[i] = rowGlobalElts[i];
      //(&*mySortedGlobalElts)[i] = (&*myGlobalElts)[localStart+i];
      //(&*mySortedPidList)[i] = (&*pidList)[localStart+i];
    }
    for (int i=leng; i< localEnd; i++) {
      (mySortedGlobalElts)[i] = (myGlobalElts)[i-leng];
    }
    for (int i=localEnd; i<numMyElts; i++) {
      (mySortedGlobalElts)[i] = (myGlobalElts)[i];
    }

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: time to sort col map arrays (cp=1) %2.3e\n", subdomainID, t1-t0);
    t0=t1;
#endif
    //FIXME timing

    int indexBase = colMap_->IndexBase();
    if (colMap_) delete colMap_;
    colMap_ = new Epetra_Map(-1,numMyElts,mySortedGlobalElts,indexBase,Comm());

    delete nodeMap;
    delete [] myGlobalElts;
    delete [] pidList;
    delete [] mySortedGlobalElts;

  } //try
  catch(...) {
    printf("** * Ifpack_OverlappingRowmatrix ctor: problem creating column map * **\n\n");
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++ end of sort
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time to sort col map arrays (cp=2) %2.3e\n", subdomainID, t1-t0);
  t0=t1;
#endif
  //FIXME timing


/*

   FIXME
   Does the column map need to be sorted for the overlapping matrix?

   The column map *must* be sorted:

        first locals
        then ghosts

   The ghosts must be further sorted so that they are contiguous by owning processor

  int* RemoteSizeList = 0
  int* RemoteColIndices = ColIndices.Values() + NumLocalColGIDs; // Points to back end of ColIndices

  EPETRA_CHK_ERR(DomainMap.RemoteIDList(NumRemoteColGIDs, RemoteColIndices, PIDList.Values(), 0, RemoteSizeList));
  Epetra_Util epetraUtil;
  SortLists[0] = RemoteColIndices;
  SortLists[1] = RemoteSizeList;
  epetraUtil.Sort(true, NumRemoteColGIDs, PIDList.Values(), 0, 0, NLists, SortLists);
*/

  // now build the map corresponding to all the external nodes
  // (with respect to A().RowMatrixRowMap().
  ExtMap_ = rcp( new Epetra_Map(-1,ghostElements.size(), &ghostElements[0],0,Comm()) );
  ExtMatrix_ = rcp( new Epetra_CrsMatrix(Copy,*ExtMap_,*colMap_,0) ); 

  ExtImporter_ = rcp( new Epetra_Import(*ExtMap_,A().RowMatrixRowMap()) ); 
  ExtMatrix_->Import(A(),*ExtImporter_,Insert); 

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time to import and setup ExtMap_ %2.3e\n", subdomainID, t1-t0);
  t0=t1;
#endif
  //FIXME timing


/*
  Notes to self:    In FillComplete, the range map does not have to be 1-1 as long as
                    (row map == range map).  Ditto for the domain map not being 1-1
                    if (col map == domain map).
                    
*/

  ExtMatrix_->FillComplete( *colMap_ , *Map_ ); //FIXME wrong

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time to FillComplete B() %2.3e\n", subdomainID, t1-t0);
  t0=t1;
#endif
  //FIXME timing


  // Note: B() = *ExtMatrix_ .

  Importer_ = rcp( new Epetra_Import(*Map_,A().RowMatrixRowMap()) ); //FIXME is this right?!

  // fix indices for overlapping matrix
  NumMyRowsB_ = B().NumMyRows();
  NumMyRows_ = NumMyRowsA_ + NumMyRowsB_;  //TODO is this wrong for a subdomain on >1 processor? // should be ok
  //NumMyCols_ = NumMyRows_;  //TODO is this wrong for a subdomain on >1 processor?  // YES!!!
  //NumMyCols_ = A().NumMyCols() + B().NumMyCols();
  NumMyCols_ = B().NumMyCols();

  /*FIXME*/ //somehow B's NumMyCols is the entire subdomain (local + overlap)

  NumMyDiagonals_ = A().NumMyDiagonals() + B().NumMyDiagonals();
  
  NumMyNonzeros_ = A().NumMyNonzeros() + B().NumMyNonzeros();
  Comm().SumAll(&NumMyNonzeros_,&NumGlobalNonzeros_,1);
  MaxNumEntries_ = A().MaxNumEntries();
  
  if (MaxNumEntries_ < B().MaxNumEntries())
    MaxNumEntries_ = B().MaxNumEntries();

# ifdef HAVE_MPI
  delete nodeComm;
# endif

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time for final calcs %2.3e\n", subdomainID, t1-t0);
  fprintf(stderr,"[node %d]: total IORM ctor time %2.3e\n", subdomainID, t1-true_t0);
#endif
  //FIXME timing


} //Ifpack_OverlappingRowMatrix() ctor for more than one core

#else
# ifdef IFPACK_NODE_AWARE_CODE

// ====================================================================== 
// Constructor for the case of two or more cores per subdomain
Ifpack_OverlappingRowMatrix::
Ifpack_OverlappingRowMatrix(const RCP<const Epetra_RowMatrix>& Matrix_in,
                            int OverlapLevel_in, int nodeID)  :
  Matrix_(Matrix_in),
  OverlapLevel_(OverlapLevel_in)
{

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  double t0 = MPI_Wtime();
  double t1, true_t0=t0;
#endif
  //FIXME timing

  const int votesMax = INT_MAX;

  // should not be here if no overlap
  if (OverlapLevel_in == 0)
    IFPACK_CHK_ERRV(-1);

  // nothing to do as well with one process
  if (Comm().NumProc() == 1)
    IFPACK_CHK_ERRV(-1);

  // nodeID is the node (aka socket) number, and so is system dependent
  // nodeComm is the communicator for all the processes on a particular node
  // these processes will have the same nodeID.
# ifdef HAVE_MPI
  const Epetra_MpiComm *pComm = dynamic_cast<const Epetra_MpiComm*>( &Comm() );
  assert(pComm != NULL);
  MPI_Comm nodeMPIComm;
  MPI_Comm_split(pComm->Comm(),nodeID,pComm->MyPID(),&nodeMPIComm);
  Epetra_MpiComm *nodeComm = new Epetra_MpiComm(nodeMPIComm);
# else
  Epetra_SerialComm *nodeComm =  dynamic_cast<Epetra_MpiComm*>( &(Matrix_in->RowMatrixRowMap().Comm()) );
# endif
  
  NumMyRowsA_ = A().NumMyRows();

  RCP<Epetra_Map> TmpMap;
  RCP<Epetra_CrsMatrix> TmpMatrix; 
  RCP<Epetra_Import> TmpImporter;

  // importing rows corresponding to elements that are 
  // in ColMap, but not in RowMap 
  const Epetra_Map *RowMap; 
  const Epetra_Map *ColMap; 
  const Epetra_Map *DomainMap;

  // TODO Count #connections from nodes I own to each ghost node


  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime(); 
  fprintf(stderr,"[node %d]: time for initialization %2.3e\n", nodeID, t1-t0);
  t0=t1;
#endif
  //FIXME timing

#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: overlap hash table allocs %2.3e\n", nodeID, t1-t0);
  t0=t1;
#endif

  Teuchos::Hashtable<int,int> colMapTable(3 * A().RowMatrixColMap().NumMyElements() );

  // ghostTable holds off-node GIDs that are connected to on-node rows and can potentially be this PID's overlap
  // TODO hopefully 3 times the # column entries is a reasonable table size
  Teuchos::Hashtable<int,int> ghostTable(3 * A().RowMatrixColMap().NumMyElements() );

  /* ** ************************************************************************** ** */
  /* ** ********************** start of main overlap loop ************************ ** */
  /* ** ************************************************************************** ** */
  for (int overlap = 0 ; overlap < OverlapLevel_in ; ++overlap)
  {
    // nbTable holds node buddy GIDs that are connected to current PID's rows, i.e., GIDs that should be in the overlapped
    // matrix's column map

    if (TmpMatrix != Teuchos::null) {
      RowMap = &(TmpMatrix->RowMatrixRowMap()); 
      ColMap = &(TmpMatrix->RowMatrixColMap()); 
      DomainMap = &(TmpMatrix->OperatorDomainMap());
    }
    else {
      RowMap = &(A().RowMatrixRowMap()); 
      ColMap = &(A().RowMatrixColMap()); 
      DomainMap = &(A().OperatorDomainMap());
    }

#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap pointer pulls %2.3e\n", nodeID, t1-t0);
    t0=t1;
#endif
    
    // For each column ID, determine the owning node (as opposed to core)
    // ID of the corresponding row.
    Epetra_IntVector colIdList( *ColMap );
    Epetra_IntVector rowIdList(*DomainMap);
    rowIdList.PutValue(nodeID);  
    Teuchos::RCP<Epetra_Import> nodeIdImporter = rcp(new Epetra_Import( *ColMap, *DomainMap ));

#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap intvector/importer allocs %2.3e\n", nodeID, t1-t0);
    t0=t1;
#endif
    
    colIdList.Import(rowIdList,*nodeIdImporter,Insert);
    
    int size = ColMap->NumMyElements() - RowMap->NumMyElements(); 
    int count = 0; 

#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap col/row ID imports %2.3e\n", nodeID, t1-t0);
    t0=t1;
#endif

    
    // define the set of off-node rows that are in ColMap but not in RowMap
    // This naturally does not include off-core rows that are on the same node as me, i.e., node buddy rows.
    for (int i = 0 ; i < ColMap->NumMyElements() ; ++i) {
      int GID = ColMap->GID(i); 
      int myvotes=0;
      if (ghostTable.containsKey(GID)) myvotes = ghostTable.get(GID);
      if ( colIdList[i] != nodeID && myvotes < votesMax) // don't include anybody found in a previous round
      {
        int votes;
        if (ghostTable.containsKey(GID)) {
          votes = ghostTable.get(GID);
          votes++;
          ghostTable.put(GID,votes);
        } else {
          ghostTable.put(GID,1);
        }
      }
    } //for (int i = 0 ; i < ColMap->NumMyElements() ; ++i)

    Teuchos::Array<int> gidsarray,votesarray;
    ghostTable.arrayify(gidsarray,votesarray);
    int *gids = gidsarray.getRawPtr();
    int *votes = votesarray.getRawPtr();

    /*
       This next bit of code decides which node buddy (NB) gets which ghost points.  Everyone sends their
       list of ghost points to pid 0 of the local subcommunicator.  Pid 0 decides who gets what:

          - if a ghost point is touched by only one NB, that NB gets the ghost point
          - if two or more NBs share a ghost point, whichever NB has the most connections to the ghost
            point gets it.
    */

#   ifdef HAVE_MPI  //FIXME What if we build in serial?!  This file will likely not compile.
    int *cullList;
    int ncull;
    //int mypid = nodeComm->MyPID();

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap before culling %2.3e\n", nodeID, t1-t0);
    t0=t1;
#endif
    //FIXME timing


    if (nodeComm->MyPID() == 0)
    {
      // Figure out how much pid 0 is to receive
      MPI_Status status;
      int *lengths = new int[nodeComm->NumProc()+1];
      lengths[0] = 0;
      lengths[1] = ghostTable.size();
      for (int i=1; i<nodeComm->NumProc(); i++) {
        int leng;
        MPI_Recv( &leng, 1, MPI_INT, i, MPI_ANY_TAG, nodeComm->Comm(), &status);
        lengths[i+1] = lengths[i] + leng;
      }
      int total = lengths[nodeComm->NumProc()];

      int* ghosts = new int[total];
      for (int i=0; i<total; i++) ghosts[i] = -9;
      int *round  = new int[total];
      int *owningpid  = new int[total];

      for (int i=1; i<nodeComm->NumProc(); i++) {
        int count = lengths[i+1] - lengths[i];
        MPI_Recv( ghosts+lengths[i], count, MPI_INT, i, MPI_ANY_TAG, nodeComm->Comm(), &status);
        MPI_Recv( round+lengths[i], count, MPI_INT, i, MPI_ANY_TAG, nodeComm->Comm(), &status);
      }

      // put in pid 0's info
      for (int i=0; i<lengths[1]; i++) {
        ghosts[i] = gids[i];
        round[i] = votes[i];
        owningpid[i] = 0;
      }

      // put in the pid associated with each ghost
      for (int j=1; j<nodeComm->NumProc(); j++) {
        for (int i=lengths[j]; i<lengths[j+1]; i++) {
          owningpid[i] = j;
        }
      }

      // sort everything based on the ghost gids
      int* roundpid[2];
      roundpid[0] = round; roundpid[1] = owningpid;
      Epetra_Util epetraUtil;
      epetraUtil.Sort(true,total,ghosts,0,0,2,roundpid);

      //set up arrays that get sent back to node buddies and that tell them what ghosts to cull
      int *nlosers = new int[nodeComm->NumProc()];
      int** losers = new int*[nodeComm->NumProc()];
      for (int i=0; i<nodeComm->NumProc(); i++) {
        nlosers[i] = 0;
        losers[i] = new int[ lengths[i+1]-lengths[i] ];
      }

      // Walk through ghosts array and and for each sequence of duplicate ghost GIDs, choose just one NB to keep it.
      // The logic is pretty simple.  The ghost list is sorted, so all duplicate PIDs are together.
      // The list is traversed.  As duplicates are found, node pid 0 keeps track of the current "winning"
      // pid.  When a pid is determined to have "lost" (less votes/connections to the current GID), the
      // GID is added to that pid's list of GIDs to be culled.  At the end of the repeated sequence, we have
      // a winner, and other NBs know whether they need to delete it from their import list.
      int max=0;   //for duplicated ghosts, index of pid with most votes and who hence keeps the ghost.
                   // TODO to break ties randomly

      for (int i=1; i<total; i++) {
        if (ghosts[i] == ghosts[i-1]) {
          int idx = i; // pid associated with idx is current "loser"
          if (round[i] > round[max]) {
            idx = max;
            max=i;
          }
          int j = owningpid[idx];
          losers[j][nlosers[j]++] = ghosts[idx];
        } else {
          max=i;
        }
      } //for (int i=1; i<total; i++)

      delete [] round;
      delete [] ghosts;
      delete [] owningpid;

      // send the arrays of ghost GIDs to be culled back to the respective node buddies
      for (int i=1; i<nodeComm->NumProc(); i++) {
        MPI_Send( nlosers+i, 1, MPI_INT, i, 8675, nodeComm->Comm());
        MPI_Send( losers[i], nlosers[i], MPI_INT, i, 8675, nodeComm->Comm());
      }

      //FIXME Unnecessary memory allocation and copying, but makes culling code cleaner
      //Could we stick this info into gids and votes, since neither is being used anymore?
      //TODO Instead of using "losers" arrays, just use "cullList" as in the else clause
      ncull = nlosers[0];
      cullList = new int[ncull+1];
      for (int i=0; i<nlosers[0]; i++)
        cullList[i] = losers[0][i];

      for (int i=0; i<nodeComm->NumProc(); i++)
        delete [] losers[i];

      delete [] lengths;
      delete [] losers;
      delete [] nlosers;

    } else { //everyone but pid 0

      // send to node pid 0 all ghosts that this pid could potentially import
      int hashsize = ghostTable.size();
      MPI_Send( &hashsize, 1, MPI_INT, 0, 8675, nodeComm->Comm());
      MPI_Send( gids, hashsize, MPI_INT, 0, 8675, nodeComm->Comm());
      MPI_Send( votes, hashsize, MPI_INT, 0, 8675, nodeComm->Comm());

      // receive the ghost GIDs that should not be imported (subset of the list sent off just above)
      MPI_Status status;
      MPI_Recv( &ncull, 1, MPI_INT, 0, 8675, nodeComm->Comm(), &status);
      cullList = new int[ncull+1];
      MPI_Recv( cullList, ncull, MPI_INT, 0, 8675, nodeComm->Comm(), &status);
    }


    //TODO clean out hash table after each time through overlap loop   4/1/07 JJH done moved both hash tables to local scope

    // Remove from my row map all off-node ghosts that will be imported by a node buddy.
    for (int i=0; i<ncull; i++) {
      try{ghostTable.remove(cullList[i]);}

      catch(...) {
        printf("pid %d: In OverlappingRowMatrix ctr, problem removing ghost elt %d from ghostTable\n",
               Comm().MyPID(),cullList[i]);
        fflush(stdout);
      }
    }//for

    delete [] cullList;

    // Save off the remaining ghost GIDs from the current overlap round.
    // These are off-node GIDs (rows) that I will import.
    gidsarray.clear(); votesarray.clear();
    ghostTable.arrayify(gidsarray,votesarray);

    vector<int> list(size); 
    count=0;
    for (int i=0; i<ghostTable.size(); i++) {
      // if votesarray[i] == votesmax, then this GID was found during a previous overlap round
      if (votesarray[i] < votesMax) {
        list[count] = gidsarray[i];
        ghostTable.put(gidsarray[i],votesMax); //set this GID's vote to the maximum so that this PID alway wins in any subsequent overlap tiebreaking.
        count++;
      }
    }

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap duplicate removal %2.3e\n", nodeID, t1-t0);
    t0=t1;
#endif
    //FIXME timing

#   endif //ifdef HAVE_MPI

    TmpMap = rcp( new Epetra_Map(-1,count, &list[0],0,Comm()) );

    TmpMatrix = rcp( new Epetra_CrsMatrix(Copy,*TmpMap,0) ); 

    TmpImporter = rcp( new Epetra_Import(*TmpMap,A().RowMatrixRowMap()) ); 

    TmpMatrix->Import(A(),*TmpImporter,Insert); 
    TmpMatrix->FillComplete(A().OperatorDomainMap(),*TmpMap); 

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap TmpMatrix fillcomplete %2.3e\n", nodeID, t1-t0);
    t0=t1;
#endif
    //FIXME timing


    // These next two imports get the GIDs that need to go into the column map of the overlapped matrix.

    // For each column ID in the overlap, determine the owning node (as opposed to core)
    // ID of the corresponding row.  Save those column IDs whose owning node is the current one.
    // This should get all the imported ghost GIDs.
    Epetra_IntVector ov_colIdList( TmpMatrix->ColMap() );
    ov_colIdList.PutValue(-1);
    Epetra_IntVector ov_rowIdList( TmpMatrix->RowMap() );
    ov_rowIdList.PutValue(nodeID);  
    Teuchos::RCP<Epetra_Import> ov_nodeIdImporter = rcp(new Epetra_Import( TmpMatrix->ColMap(), TmpMatrix->RowMap()));
    ov_colIdList.Import(ov_rowIdList,*ov_nodeIdImporter,Insert);

    for (int i=0 ; i<ov_colIdList.MyLength(); i++) {
       if (ov_colIdList[i] == nodeID) {
         int GID = (ov_colIdList.Map()).GID(i);
         colMapTable.put(GID,1);
       }
    }

    // Do a second import of the owning node ID from A's rowmap to TmpMat's column map.  This ensures that
    // all GIDs that belong to a node buddy and are in a ghost row's sparsity pattern will be in the final
    // overlapped matrix's column map.
    ov_colIdList.PutValue(-1);
    Epetra_IntVector ArowIdList( A().RowMatrixRowMap() );
    ArowIdList.PutValue(nodeID);
    nodeIdImporter = rcp(new Epetra_Import( TmpMatrix->ColMap(), A().RowMatrixRowMap() ));
    ov_colIdList.Import(ArowIdList,*nodeIdImporter,Insert);

    for (int i=0 ; i<ov_colIdList.MyLength(); i++) {
       if (ov_colIdList[i] == nodeID) {
         int GID = (ov_colIdList.Map()).GID(i);
         colMapTable.put(GID,1);
       }
    }

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: overlap 2 imports to fix up colmap %2.3e\n", nodeID, t1-t0);
    t0=t1;
#endif
    //FIXME timing

  } //for (int overlap = 0 ; overlap < OverlapLevel_in ; ++overlap)

  /* ** ************************************************************************ ** */
  /* ** ********************** end of main overlap loop ************************ ** */
  /* ** ************************************************************************ ** */

  // off-node GIDs that will be in the overlap
  vector<int> ghostElements; 

  Teuchos::Array<int> gidsarray,votesarray;
  ghostTable.arrayify(gidsarray,votesarray);
  for (int i=0; i<ghostTable.size(); i++) {
    ghostElements.push_back(gidsarray[i]);
  }

    for (int i = 0 ; i < A().RowMatrixColMap().NumMyElements() ; ++i) {
      int GID = A().RowMatrixColMap().GID(i);
      // Remove any entries that are in A's original column map
      if (colMapTable.containsKey(GID)) {
        try{colMapTable.remove(GID);}
        catch(...) {
          printf("pid %d: In OverlappingRowMatrix ctr, problem removing entry %d from colMapTable\n", Comm().MyPID(),GID);
          fflush(stdout);
        }
      }
    }

    // GIDs that will be in the overlapped matrix's column map
    vector<int> colMapElements; 

  gidsarray.clear(); votesarray.clear();
  colMapTable.arrayify(gidsarray,votesarray);
  for (int i=0; i<colMapTable.size(); i++)
    colMapElements.push_back(gidsarray[i]);

/*
   We need two maps here.  The first is the row map, which we've got by using my original rows
   plus whatever I've picked up in ghostElements.

   The second is a column map.  This map should include my rows, plus any columns that could come from node buddies.
   These GIDs come from the std:array colMapElements, which in turn comes from colMapTable.
   This map should *not* omit ghosts that have been imported by my node buddies, i.e., for any row that I own,
   the stencil should include all column GIDs (including imported ghosts) that are on the node.
*/

  // build the row map containing all the nodes (original matrix + off-node matrix)
  vector<int> rowList(NumMyRowsA_ + ghostElements.size());
  for (int i = 0 ; i < NumMyRowsA_ ; ++i)
    rowList[i] = A().RowMatrixRowMap().GID(i);
  for (int i = 0 ; i < (int)ghostElements.size() ; ++i)
    rowList[i + NumMyRowsA_] = ghostElements[i];

  // row map for the overlapped matrix (local + overlap)
  Map_ = rcp( new Epetra_Map(-1, NumMyRowsA_ + ghostElements.size(), &rowList[0], 0, Comm()) );

  // build the column map for the overlapping matrix
  //vector<int> colList(colMapElements.size());
  // column map for the overlapped matrix (local + overlap)
  //colMap_ = rcp( new Epetra_Map(-1, colMapElements.size(), &colList[0], 0, Comm()) );
  //for (int i = 0 ; i < (int)colMapElements.size() ; i++)
  //  colList[i] = colMapElements[i];
  vector<int> colList(A().RowMatrixColMap().NumMyElements() + colMapElements.size());
  int nc = A().RowMatrixColMap().NumMyElements();
  for (int i = 0 ; i < nc; i++)
    colList[i] = A().RowMatrixColMap().GID(i);
  for (int i = 0 ; i < (int)colMapElements.size() ; i++)
    colList[nc+i] = colMapElements[i];

  // column map for the overlapped matrix (local + overlap)
  //colMap_ = rcp( new Epetra_Map(-1, A().RowMatrixColMap().NumMyElements() + colMapElements.size(), &colList[0], 0, Comm()) );
  colMap_ = new Epetra_Map(-1, A().RowMatrixColMap().NumMyElements() + colMapElements.size(), &colList[0], 0, Comm());
    

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time to init B() col/row maps %2.3e\n", nodeID, t1-t0);
  t0=t1;
#endif
  //FIXME timing

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++ start of sort
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // build the column map, but don't use a copy constructor, b/c local communicator SubComm_ is
  // different from that of Matrix.
  try {
    // build row map based on the local communicator.  We need this temporarily to build the column map.
    Epetra_Map* nodeMap = new Epetra_Map(-1,NumMyRowsA_ + ghostElements.size(),&rowList[0],0,*nodeComm);
    int numMyElts = colMap_->NumMyElements();
    assert(numMyElts!=0);

    // The column map *must* be sorted: first locals, then ghosts.  The ghosts
    // must be further sorted so that they are contiguous by owning processor.

    // For each GID in column map, determine owning PID in local communicator.
    int* myGlobalElts = new int[numMyElts];
    colMap_->MyGlobalElements(myGlobalElts);
    int* pidList = new int[numMyElts];
    nodeMap->RemoteIDList(numMyElts, myGlobalElts, pidList, 0);

    // First sort column map GIDs based on the owning PID in local communicator.
    Epetra_Util Util;
    int *tt[1];
    tt[0] = myGlobalElts;
    Util.Sort(true, numMyElts, pidList, 0, (double**)0, 1, tt);

    // For each remote PID, sort the column map GIDs in ascending order.
    // Don't sort the local PID's entries.
    int localStart=0;
    while (localStart<numMyElts) {
      int currPID = (pidList)[localStart];
      int i=localStart;
      while (i<numMyElts && (pidList)[i] == currPID) i++;
      if (currPID != nodeComm->MyPID())
        Util.Sort(true, i-localStart, (myGlobalElts)+localStart, 0, 0, 0, 0);
      localStart = i;
    }

    // now move the local entries to the front of the list
    localStart=0;
    while (localStart<numMyElts && (pidList)[localStart] != nodeComm->MyPID()) localStart++;
    assert(localStart != numMyElts);
    int localEnd=localStart;
    while (localEnd<numMyElts && (pidList)[localEnd] == nodeComm->MyPID()) localEnd++;
    int* mySortedGlobalElts = new int[numMyElts];
    int leng = localEnd - localStart;
    /* This is a little gotcha.  It appears that the ordering of the column map's local entries
       must be the same as that of the domain map's local entries.  See the comment in method
       MakeColMap() in Epetra_CrsGraph.cpp, line 1072. */
    int *rowGlobalElts =  nodeMap->MyGlobalElements();

    /*TODO TODO TODO NTS rows 68 and 83 show up as local GIDs in rowGlobalElts for both pids 0 & 1 on node 0.
                    This seems to imply that there is something wrong with rowList!!
                    It is ok for the overlapped matrix's row map to have repeated entries (it's overlapped, after all).
                    But ... on a node, there should be no repeats.
                    Ok, here's the underlying cause.  On node 0, gpid 1 finds 83 on overlap round 0.  gpid 0 finds 83
                    on overlap round 1.  This shouldn't happen .. once a node buddy "owns" a row, no other node buddy
                    should ever import ("own") that row.

                    Here's a possible fix.  Extend tie breaking across overlap rounds.  This means sending to lpid 0
                    a list of *all* rows you've imported (this round plus previous rounds) for tie-breaking
                    If a nb won a ghost row in a previous round, his votes for that ghost row should be really high, i.e.,
                    he should probably win that ghost row forever.
                    7/17/09

      7/28/09 JJH I believe I've fixed this, but I thought it might be helpful to have this comment in here, in case I missed something.
    */

    //move locals to front of list
    for (int i=0; i<leng; i++) {
      /* //original code */
      (mySortedGlobalElts)[i] = rowGlobalElts[i];
      //(&*mySortedGlobalElts)[i] = (&*myGlobalElts)[localStart+i];
      //(&*mySortedPidList)[i] = (&*pidList)[localStart+i];
    }
    for (int i=leng; i< localEnd; i++) {
      (mySortedGlobalElts)[i] = (myGlobalElts)[i-leng];
    }
    for (int i=localEnd; i<numMyElts; i++) {
      (mySortedGlobalElts)[i] = (myGlobalElts)[i];
    }

    //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
    t1 = MPI_Wtime();
    fprintf(stderr,"[node %d]: time to sort col map arrays (cp=1) %2.3e\n", nodeID, t1-t0);
    t0=t1;
#endif
    //FIXME timing

    int indexBase = colMap_->IndexBase();
    if (colMap_) delete colMap_;
    colMap_ = new Epetra_Map(-1,numMyElts,mySortedGlobalElts,indexBase,Comm());

    delete nodeMap;
    delete [] myGlobalElts;
    delete [] pidList;
    delete [] mySortedGlobalElts;

  } //try
  catch(...) {
    printf("** * Ifpack_OverlappingRowmatrix ctor: problem creating column map * **\n\n");
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++ end of sort
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time to sort col map arrays (cp=2) %2.3e\n", nodeID, t1-t0);
  t0=t1;
#endif
  //FIXME timing


/*

   FIXME
   Does the column map need to be sorted for the overlapping matrix?

   The column map *must* be sorted:

        first locals
        then ghosts

   The ghosts must be further sorted so that they are contiguous by owning processor

  int* RemoteSizeList = 0
  int* RemoteColIndices = ColIndices.Values() + NumLocalColGIDs; // Points to back end of ColIndices

  EPETRA_CHK_ERR(DomainMap.RemoteIDList(NumRemoteColGIDs, RemoteColIndices, PIDList.Values(), 0, RemoteSizeList));
  Epetra_Util epetraUtil;
  SortLists[0] = RemoteColIndices;
  SortLists[1] = RemoteSizeList;
  epetraUtil.Sort(true, NumRemoteColGIDs, PIDList.Values(), 0, 0, NLists, SortLists);
*/

  // now build the map corresponding to all the external nodes
  // (with respect to A().RowMatrixRowMap().
  ExtMap_ = rcp( new Epetra_Map(-1,ghostElements.size(), &ghostElements[0],0,Comm()) );
  ExtMatrix_ = rcp( new Epetra_CrsMatrix(Copy,*ExtMap_,*colMap_,0) ); 

  ExtImporter_ = rcp( new Epetra_Import(*ExtMap_,A().RowMatrixRowMap()) ); 
  ExtMatrix_->Import(A(),*ExtImporter_,Insert); 

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time to import and setup ExtMap_ %2.3e\n", nodeID, t1-t0);
  t0=t1;
#endif
  //FIXME timing


/*
  Notes to self:    In FillComplete, the range map does not have to be 1-1 as long as
                    (row map == range map).  Ditto for the domain map not being 1-1
                    if (col map == domain map).
                    
*/

  ExtMatrix_->FillComplete( *colMap_ , *Map_ ); //FIXME wrong

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time to FillComplete B() %2.3e\n", nodeID, t1-t0);
  t0=t1;
#endif
  //FIXME timing


  // Note: B() = *ExtMatrix_ .

  Importer_ = rcp( new Epetra_Import(*Map_,A().RowMatrixRowMap()) ); //FIXME is this right?!

  // fix indices for overlapping matrix
  NumMyRowsB_ = B().NumMyRows();
  NumMyRows_ = NumMyRowsA_ + NumMyRowsB_;  //TODO is this wrong for a subdomain on >1 processor? // should be ok
  //NumMyCols_ = NumMyRows_;  //TODO is this wrong for a subdomain on >1 processor?  // YES!!!
  //NumMyCols_ = A().NumMyCols() + B().NumMyCols();
  NumMyCols_ = B().NumMyCols();

  /*FIXME*/ //somehow B's NumMyCols is the entire subdomain (local + overlap)

  NumMyDiagonals_ = A().NumMyDiagonals() + B().NumMyDiagonals();
  
  NumMyNonzeros_ = A().NumMyNonzeros() + B().NumMyNonzeros();
  Comm().SumAll(&NumMyNonzeros_,&NumGlobalNonzeros_,1);
  MaxNumEntries_ = A().MaxNumEntries();
  
  if (MaxNumEntries_ < B().MaxNumEntries())
    MaxNumEntries_ = B().MaxNumEntries();

# ifdef HAVE_MPI
  delete nodeComm;
# endif

  //FIXME timing
#ifdef IFPACK_OVA_TIME_BUILD
  t1 = MPI_Wtime();
  fprintf(stderr,"[node %d]: time for final calcs %2.3e\n", nodeID, t1-t0);
  fprintf(stderr,"[node %d]: total IORM ctor time %2.3e\n", nodeID, t1-true_t0);
#endif
  //FIXME timing


} //Ifpack_OverlappingRowMatrix() ctor for more than one core



# endif //ifdef IFPACK_NODE_AWARE_CODE
#endif // ifdef IFPACK_SUBCOMM_CODE



// ======================================================================



// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  IFPACK_CHK_ERR(-1);
}


// ======================================================================
int Ifpack_OverlappingRowMatrix::
Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();
  vector<int> Ind(MaxNumEntries_);
  vector<double> Val(MaxNumEntries_);

  Y.PutScalar(0.0);

  // matvec with A (local rows)
  for (int i = 0 ; i < NumMyRowsA_ ; ++i) {
    for (int k = 0 ; k < NumVectors ; ++k) {
      int Nnz;
      IFPACK_CHK_ERR(A().ExtractMyRowCopy(i,MaxNumEntries_,Nnz, 
                                          &Val[0], &Ind[0]));
      for (int j = 0 ; j < Nnz ; ++j) {
        Y[k][i] += Val[j] * X[k][Ind[j]];
      }
    }
  }

  // matvec with B (overlapping rows)
  for (int i = 0 ; i < NumMyRowsB_ ; ++i) {
    for (int k = 0 ; k < NumVectors ; ++k) {
      int Nnz;
      IFPACK_CHK_ERR(B().ExtractMyRowCopy(i,MaxNumEntries_,Nnz, 
                                          &Val[0], &Ind[0]));
      for (int j = 0 ; j < Nnz ; ++j) {
        Y[k][i + NumMyRowsA_] += Val[j] * X[k][Ind[j]];
      }
    }
  }
  return(0);
}

#endif
