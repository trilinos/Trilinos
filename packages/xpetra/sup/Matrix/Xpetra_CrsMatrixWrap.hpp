// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_CRSMATRIXWRAP_HPP
#define XPETRA_CRSMATRIXWRAP_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"

#include "Xpetra_Matrix.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

/** \file Xpetra_CrsMatrixWrap.hpp

  Declarations for the class Xpetra::CrsMatrixWrap.
*/
namespace Xpetra {

  typedef std::string viewLabel_t;

/*!
  @class CrsMatrixWrap
  @brief Concrete implementation of Xpetra::Matrix.
*/
template <class Scalar, 
          class LocalOrdinal  = int, 
          class GlobalOrdinal = LocalOrdinal, 
          class Node          = Kokkos::DefaultNode::DefaultNodeType, 
          class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
class CrsMatrixWrap : public Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Matrix;
  typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsGraph;
#ifdef HAVE_XPETRA_TPETRA
  typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
#endif
  typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
  typedef Xpetra::MatrixView<LocalOrdinal, GlobalOrdinal, Node> MatrixView;

public:
  
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor specifying fixed number of entries for each row.
  CrsMatrixWrap(const RCP<const Map> &rowMap, size_t maxNumEntriesPerRow, Xpetra::ProfileType pftype = Xpetra::DynamicProfile)
    : finalDefaultView_(false)
  {
    // Set matrix data
    matrixData_ = CrsMatrixFactory::Build(rowMap, maxNumEntriesPerRow, pftype);

    // Default view
    CreateDefaultView();
  }

  //! Constructor specifying (possibly different) number of entries in each row.
  CrsMatrixWrap(const RCP<const Map> &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = Xpetra::DynamicProfile)
    : finalDefaultView_(false)
  {
    // Set matrix data
    matrixData_ = CrsMatrixFactory::Build(rowMap, NumEntriesPerRowToAlloc, pftype);
    
    // Default view
    CreateDefaultView();
  }
  
  CrsMatrixWrap(RCP<CrsMatrix> &matrix)
    : finalDefaultView_(matrix->isFillComplete())
  {
    // Set matrix data
    matrixData_ = matrix;

    // Default view
    CreateDefaultView();
  }
  
  //! Destructor
  virtual ~CrsMatrixWrap() {}

  //@}


  //! @name Insertion/Removal Methods
  //@{ 

  //! Insert matrix entries, using global IDs.
  /** All index values must be in the global space. 
      \pre \c globalRow exists as an ID in the global row map
      \pre <tt>isLocallyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isGloballyIndexed() == true</tt>

      \note If \c globalRow does not belong to the matrix on this node, then it will be communicated to the appropriate node when globalAssemble() is called (which will, at the latest, occur during the next call to fillComplete().) Otherwise, the entries will be inserted in the local matrix. 
      \note If the matrix row already contains values at the indices corresponding to values in \c cols, then the new values will be summed with the old values; this may happen at insertion or during the next call to fillComplete().
      \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where cols[i] belongs to the column map on this node will be inserted into the matrix.
  */
  void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) { 
    matrixData_->insertGlobalValues(globalRow, cols, vals);
  }

  //! Insert matrix entries, using local IDs.
  /** All index values must be in the local space.
      \pre \c localRow exists as an ID in the global row map
      \pre <tt>isGloballyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isLocallyIndexed() == true</tt>
  */
  void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
    matrixData_->insertLocalValues(localRow, cols, vals);
  }

  //@}

  //! @name Transformational Methods
  //@{ 

    /*! \brief Signal that data entry is complete, specifying domain and range maps.

    Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */ 
  void fillComplete(const RCP<const Map> &domainMap, const RCP<const Map> &rangeMap, const RCP<Teuchos::ParameterList> &params = null) {
    matrixData_->fillComplete(domainMap, rangeMap, params);

    // Update default view with the colMap because colMap can be <tt>null</tt> until fillComplete() is called.
    updateDefaultView();
  }

  /*! \brief Signal that data entry is complete. 

  Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

  \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

  \pre  <tt>isFillActive() == true<tt>
  \pre <tt>isFillComplete()() == false<tt>

  \post <tt>isFillActive() == false<tt>
  \post <tt>isFillComplete() == true<tt>
  \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
  */
  //TODO : Get ride of "Tpetra"::OptimizeOption
  void fillComplete(const RCP<ParameterList> &params = null) {
    matrixData_->fillComplete(params);

    // Update default view with the colMap because colMap can be <tt>null</tt> until fillComplete() is called.
    updateDefaultView();
  }

  //@}

  //! Returns the number of global rows in this matrix.
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumRows() const { 
    return matrixData_->getGlobalNumRows();
  }

  //! \brief Returns the number of global columns in the matrix.
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumCols() const {
    return matrixData_->getGlobalNumCols();
  }

  //! Returns the number of matrix rows owned on the calling node.
  size_t getNodeNumRows() const {
    return matrixData_->getNodeNumRows();
  }

  //! Returns the global number of entries in this matrix.
  global_size_t getGlobalNumEntries() const {
    return matrixData_->getGlobalNumEntries();
  }

  //! Returns the local number of entries in this matrix.
  size_t getNodeNumEntries() const {
    return matrixData_->getNodeNumEntries();
  }

  //! Returns the current number of entries on this node in the specified local row.
  /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
  size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    return matrixData_->getNumEntriesInLocalRow(localRow);
  }

  //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumDiags() const {
    return matrixData_->getGlobalNumDiags();
  }

  //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
  /** Undefined if isFillActive().
   */
  size_t getNodeNumDiags() const {
    return matrixData_->getNodeNumDiags();
  }

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  /** Undefined if isFillActive().
   */
  size_t getGlobalMaxNumRowEntries() const {
    return matrixData_->getGlobalMaxNumRowEntries();
  }

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  /** Undefined if isFillActive().
   */
  size_t getNodeMaxNumRowEntries() const {
    return matrixData_->getNodeMaxNumRowEntries();
  }

  //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
  bool isLocallyIndexed() const {
    return matrixData_->isLocallyIndexed();
  }

  //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
  bool isGloballyIndexed() const {
    return matrixData_->isGloballyIndexed();
  }

  //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
  bool isFillComplete() const {
    return matrixData_->isFillComplete();
  }

  //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices - (Out) Local column indices corresponding to values.
    \param Values - (Out) Matrix values.
    \param NumIndices - (Out) Number of indices.
    
    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
    returned as OrdinalTraits<size_t>::invalid().
    
    \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
  */
  void getLocalRowCopy(LocalOrdinal LocalRow, 
                       const ArrayView<LocalOrdinal> &Indices, 
                       const ArrayView<Scalar> &Values,
                       size_t &NumEntries
                       ) const {
    matrixData_->getLocalRowCopy(LocalRow, Indices, Values, NumEntries);
  }
  
  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*!
    \param GlobalRow - (In) Global row number for which indices are desired.
    \param Indices   - (Out) Global column indices corresponding to values.
    \param Values    - (Out) Row values
    \pre <tt>isLocallyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
  */
  void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const {
     matrixData_->getGlobalRowView(GlobalRow, indices, values);
  }

  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices  - (Out) Global column indices corresponding to values.
    \param Values   - (Out) Row values
    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

    Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
  */
  void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const {
     matrixData_->getLocalRowView(LocalRow, indices, values);
  }

  //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
    the zero and non-zero diagonals owned by this node. */
  void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const {
    matrixData_->getLocalDiagCopy(diag);
  }

  //! Get Frobenius norm of the matrix
  typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const {
    return matrixData_->getFrobeniusNorm();
  }

  //@}

  //! @name Advanced Matrix-vector multiplication and solve methods
  //@{

  //! Multiplies this matrix by a MultiVector.
  /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

  Both are required to have constant stride, and they are not permitted to ocupy overlapping space. No runtime checking will be performed in a non-debug build.

  This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that multiply() not be called directly; instead, use the CrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
  If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
  will be accumulated into \c Y.
  */
  //TODO virtual=0 // TODO: Add default parameters ?
//   void multiply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, Scalar alpha, Scalar beta) const {
//      matrixData_->multiply(X, Y, trans, alpha, beta);
//   }

  //@}

  //! @name Methods implementing Matrix
  //@{ 

  //! \brief Computes the sparse matrix-multivector multiplication.
  /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
    - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
  */
  void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha = ScalarTraits<Scalar>::one(),
             Scalar beta = ScalarTraits<Scalar>::zero()) const {

    matrixData_->apply(X,Y,mode,alpha,beta);
  }
  
  //! \brief Returns the Map associated with the domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
    return matrixData_->getDomainMap();
  }

  //! Returns the Map associated with the domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
    return matrixData_->getRangeMap();
  }

  //! \brief Returns the Map that describes the column distribution in this matrix.
  //! This might be <tt>null</tt> until fillComplete() is called.
  const RCP<const Map> & getColMap() const { return getColMap(Matrix::GetCurrentViewLabel()); }
  
  //! \brief Returns the Map that describes the column distribution in this matrix.
  const RCP<const Map> & getColMap(viewLabel_t viewLabel) const { 
    TEUCHOS_TEST_FOR_EXCEPTION(Matrix::operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.GetColMap(): view '" + viewLabel + "' does not exist.");
    updateDefaultView(); // If CrsMatrix::fillComplete() have been used instead of CrsMatrixWrap::fillComplete(), the default view is updated.
    return Matrix::operatorViewTable_.get(viewLabel)->GetColMap(); 
  }
  
  //@}

  //! @name Overridden from Teuchos::Describable 
  //@{
  
  /** \brief Return a simple one-line description of this object. */
  std::string description() const {
    return "Xpetra_CrsMatrixWrap.description()";
  }
  
  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { 
    //     Teuchos::EVerbosityLevel vl = verbLevel;
    //     if (vl == VERB_DEFAULT) vl = VERB_LOW;
    //     RCP<const Comm<int> > comm = this->getComm();
    //     const int myImageID = comm->getRank(),
    //       numImages = comm->getSize();
    
    //     if (myImageID == 0) out << this->description() << std::endl; 
    
    matrixData_->describe(out,verbLevel);

    // Teuchos::OSTab tab(out);
  }

  // JG: Added:

  //! Returns the CrsGraph associated with this matrix. 
  RCP<const CrsGraph> getCrsGraph() const { return matrixData_->getCrsGraph(); }

  RCP<CrsMatrix> getCrsMatrix() const {  return matrixData_; }

  //@}
  

private:

  // Default view is created after fillComplete()
  // Because ColMap might not be available before fillComplete(). 
  void CreateDefaultView() {
    
    // Create default view
    this->defaultViewLabel_ = "point";
    this->CreateView(this->GetDefaultViewLabel(), matrixData_->getRowMap(), matrixData_->getColMap());    
    
    // Set current view
    this->currentViewLabel_ = this->GetDefaultViewLabel();
  }

private:

  // The colMap can be <tt>null</tt> until fillComplete() is called. The default view of the Matrix have to be updated when fillComplete() is called. 
  // If CrsMatrix::fillComplete() have been used instead of CrsMatrixWrap::fillComplete(), the default view is updated when getColMap() is called.
  void updateDefaultView() const {
    if ((finalDefaultView_ == false) &&  matrixData_->isFillComplete() ) {
      // Update default view with the colMap
      Matrix::operatorViewTable_.get(Matrix::GetDefaultViewLabel())->SetColMap(matrixData_->getColMap());
      finalDefaultView_ = true;
    }
  }
  // The boolean finalDefaultView_ keep track of the status of the default view (= already updated or not)
  // See also CrsMatrixWrap::updateDefaultView()
  mutable bool finalDefaultView_;


  RCP<CrsMatrix> matrixData_;

}; //class Matrix

} //namespace Xpetra

#define XPETRA_CRSMATRIXWRAP_SHORT
#endif //XPETRA_CRSMATRIXWRAP_DECL_HPP

//NOTE: if CrsMatrix and VbrMatrix share a common interface for fillComplete() etc, I can move some stuff in Xpetra_Matrix.hpp
//TODO: getUnderlyingMatrix() method
