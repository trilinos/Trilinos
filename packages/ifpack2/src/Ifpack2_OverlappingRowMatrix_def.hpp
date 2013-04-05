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
	  = std::find(ExtElements.begin(),ExtElements.end(),GID); 
        if (pos == ExtElements.end()) { 
          ExtElements.push_back(GID);
          mylist[count] = GID; 
          ++count; 
        } 
      } 
    }
    
    // Allocate & import new matrices, maps, etc.
    TmpMap      = rcp(new MapType(global_invalid,mylist(0,count),Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),A_->getComm(),A_->getNode()));
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
  this->template applyTempl<Scalar,Scalar>(X, Y, mode, alpha, beta);
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
template<class DomainScalar, class RangeScalar>
void OverlappingRowMatrix<MatrixType>::applyTempl(const Tpetra::MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
						  Tpetra::MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
						  Teuchos::ETransp mode, 
						  RangeScalar alpha,
						  RangeScalar beta) const
{
  // Note: This isn't AztecOO compliant.  But neither was Ifpack's version.
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
			     "Ifpack2::OverlappingRowMatrix::apply ERROR: X.getNumVectors() != Y.getNumVectors().");

  RangeScalar zero = Teuchos::ScalarTraits<RangeScalar>::zero();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const DomainScalar> > x_ptr = X.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<RangeScalar> >        y_ptr = Y.get2dViewNonConst();
  Y.putScalar(zero);
  size_t NumVectors = Y.getNumVectors();

  for (size_t i = 0 ; i < NumMyRowsA_ ; ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    A_->getLocalRowCopy(i,Indices_(),Values_(),Nnz);
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

 for (size_t i = 0 ; i < NumMyRowsB_ ; ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    ExtMatrix_->getLocalRowCopy(i,Indices_(),Values_(),Nnz);
    if (mode==Teuchos::NO_TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][NumMyRowsA_+i] += (RangeScalar)Values_[j] * (RangeScalar)x_ptr[k][Indices_[j]];      
    }
    else if (mode==Teuchos::TRANS){
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][NumMyRowsA_+Indices_[j]] += (RangeScalar)Values_[j] * (RangeScalar)x_ptr[k][i];
    }
    else { //mode==Teuchos::CONJ_TRANS
      for (size_t j = 0 ; j < Nnz ; ++j) 
	for (size_t k = 0 ; k < NumVectors ; ++k)
	  y_ptr[k][NumMyRowsA_+Indices_[j]] += Teuchos::ScalarTraits<RangeScalar>::conjugate((RangeScalar)Values_[j]) * (RangeScalar)x_ptr[k][i];
    }
  }
}
  


// ======================================================================
template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::importMultiVector(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
							 Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &OvX, 
							 Tpetra::CombineMode CM)
{
  this->template importMultiVectorTempl<Scalar>(X, OvX, CM);
}

// ======================================================================
template<class MatrixType>
template<class OpScalar>
void OverlappingRowMatrix<MatrixType>::importMultiVectorTempl(const Tpetra::MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
							      Tpetra::MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &OvX, 
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
  this->template exportMultiVectorTempl<Scalar>(OvX, X, CM);
}

// ======================================================================
template<class MatrixType>
template<class OpScalar>
void OverlappingRowMatrix<MatrixType>::exportMultiVectorTempl(const Tpetra::MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &OvX, 
							      Tpetra::MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
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

//==========================================================================
template<class MatrixType> 
bool OverlappingRowMatrix<MatrixType>::supportsRowViews() const
{
  return false;
}


} // namespace Ifpack2

#endif // IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP
