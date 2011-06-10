#ifndef CTHULHU_CRSOPERATOR_HPP
#define CTHULHU_CRSOPERATOR_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Exceptions.hpp"

#include "Cthulhu_BlockMap.hpp"
#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_CrsGraph.hpp"
#include "Cthulhu_CrsMatrix.hpp"
#include "Cthulhu_CrsMatrixFactory.hpp"

#include "Cthulhu_Operator.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

#include "Cthulhu_Debug.hpp"

/** \file Cthulhu_Operator.hpp

  Declarations for the class Cthulhu::Operator.
*/
namespace Cthulhu {

  typedef std::string viewLabel_t;

template <class Scalar, 
          class LocalOrdinal  = int, 
          class GlobalOrdinal = LocalOrdinal, 
          class Node          = Kokkos::DefaultNode::DefaultNodeType, 
          class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
class CrsOperator : public Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

  typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Cthulhu::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
  typedef Cthulhu::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsGraph;
#ifdef HAVE_CTHULHU_TPETRA
  typedef Cthulhu::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
#endif
  typedef Cthulhu::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
  typedef Cthulhu::OperatorView<LocalOrdinal, GlobalOrdinal, Node> OperatorView;

public:
  
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  CrsOperator(const RCP<const Map> &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) 
  {
    // Set matrix data
    matrixData_ = CrsMatrixFactory::Build(rowMap, maxNumEntriesPerRow, pftype);

    // Create default view
    this->defaultViewLabel_ = "point";
    CreateView(this->GetDefaultViewLabel(), matrixData_->getRowMap(), matrixData_->getColMap());    
    
    // Set current view
    this->currentViewLabel_ = this->GetDefaultViewLabel();
    this->SwitchToDefaultView(); // JG: this line is useless (for the moment).
  }
  
  //! Destructor
  virtual ~CrsOperator() {}

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
  inline void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) { 
    matrixData_->insertGlobalValues(globalRow, cols, vals);
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
  inline void fillComplete(const RCP<const Map> &domainMap, const RCP<const Map> &rangeMap, OptimizeOption os = DoOptimizeStorage) {
    matrixData_->fillComplete(domainMap, rangeMap, os);
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
  inline void fillComplete(Cthulhu::OptimizeOption os = Cthulhu::DoOptimizeStorage) {
    matrixData_->fillComplete(os);
  }

  //@}

  //! Returns the number of global rows in this matrix.
  /** Undefined if isFillActive().
   */
  inline global_size_t getGlobalNumRows() const { 
    return matrixData_->getGlobalNumRows();
  }

  //! \brief Returns the number of global columns in the matrix.
  /** Undefined if isFillActive().
   */
  inline global_size_t getGlobalNumCols() const {
    return matrixData_->getGlobalNumCols();
  }

  //! Returns the number of matrix rows owned on the calling node.
  inline size_t getNodeNumRows() const {
    return matrixData_->getNodeNumRows();
  }

  //! Returns the global number of entries in this matrix.
  inline global_size_t getGlobalNumEntries() const {
    return matrixData_->getGlobalNumEntries();
  }

  //! Returns the local number of entries in this matrix.
  inline size_t getNodeNumEntries() const {
    return matrixData_->getNodeNumEntries();
  }

  //! Returns the current number of entries on this node in the specified local row.
  /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
  inline size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    return  matrixData_->getNumEntriesInLocalRow(localRow);
  }

  //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
  /** Undefined if isFillActive().
   */
  inline global_size_t getGlobalNumDiags() const {
    return  matrixData_->getGlobalNumDiags();
  }

  //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
  /** Undefined if isFillActive().
   */
  inline size_t getNodeNumDiags() const {
    return  matrixData_->getNodeNumDiags();
  }

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  /** Undefined if isFillActive().
   */
  inline size_t getGlobalMaxNumRowEntries() const {
    return  matrixData_->getGlobalMaxNumRowEntries();
  }

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  /** Undefined if isFillActive().
   */
  inline size_t getNodeMaxNumRowEntries() const {
    return  matrixData_->getNodeMaxNumRowEntries();
  }

  //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
  inline bool isLocallyIndexed() const {
    return  matrixData_->isLocallyIndexed();
  }

  //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
  inline bool isGloballyIndexed() const {
    return  matrixData_->isGloballyIndexed();
  }

  //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
  inline bool isFillComplete() const {
    return  matrixData_->isFillComplete();
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
    virtual void getLocalRowCopy(LocalOrdinal LocalRow, 
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
  inline void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const {
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
  inline void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const {
     matrixData_->getLocalRowView(LocalRow, indices, values);
  }

  //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
    the zero and non-zero diagonals owned by this node. */
  inline void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const {
    matrixData_->getLocalDiagCopy(diag);
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
  inline void multiply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, Scalar alpha, Scalar beta) const {
     matrixData_->multiply(X, Y, trans, alpha, beta);
  }

  //@}

  //! @name Methods implementing Operator
  //@{ 

  //! \brief Computes the sparse matrix-multivector multiplication.
  /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
    - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
  */
  virtual void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                     Teuchos::ETransp mode = Teuchos::NO_TRANS,
                     Scalar alpha = ScalarTraits<Scalar>::one(),
                     Scalar beta = ScalarTraits<Scalar>::zero()) const {

    return matrixData_->apply(X,Y,mode,alpha,beta);
  }
  
  //! \brief Returns the Map associated with the domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  inline const RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
    return matrixData_->getDomainMap();
  }

  //! Returns the Map associated with the domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  inline const RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
    return matrixData_->getRangeMap();
  }

  //@}

  //! @name Overridden from Teuchos::Describable 
  //@{
  
  /** \brief Return a simple one-line description of this object. */
  std::string description() const { 
    std::ostringstream oss;
    oss << "Cthulhu_CrsOperator.description()" << std::endl;
    return oss.str();
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

  RCP<CrsMatrix> getCrsMatrix() const { CTHULHU_DEBUG_ME; return matrixData_; }

  //@}
  
private:
  RCP<CrsMatrix> matrixData_;

}; //class Operator

} //namespace Cthulhu

#define CTHULHU_CRSOPERATOR_SHORT
#endif //CTHULHU_CRSOPERATOR_DECL_HPP

//NOTE: if CrsMatrix and VbrMatrix share a common interface for fillComplete() etc, I can move some stuff in Cthulhu_Operator.hpp
//TODO: getUnderlyingMatrix() method
