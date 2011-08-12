#ifndef CTHULHU_TPETRAVBRMATRIX_HPP
#define CTHULHU_TPETRAVBRMATRIX_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_TPETRA
#error This file should be included only if HAVE_CTHULHU_TPETRA is defined.
#endif

#include "Cthulhu_VbrMatrix.hpp"

#include <Tpetra_VbrMatrix.hpp>

namespace Cthulhu {

  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::BlockSparseOps >
  class TpetraVbrMatrix : 
    public VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalOrdinal> {
  public:
  
    //! @name Constructor/Destructor Methods
    //@{

    TpetraVbrMatrix(const Teuchos::RCP<const Tpetra::VbrMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &mtx) : mtx_(mtx) {  } //TODO

    //! Destructor
    virtual ~TpetraVbrMatrix();

    //@}

    //! @name Operator Methods
    //@{

    //! Returns the Map associated with the domain of this operator.
    /*! Note that this is a point-entry map, not a block-map.
     */
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const {  return mtx_->getDomainMap(); }

    //! Returns the Map associated with the range of this operator, which must be compatible with Y.getMap().
    /*! Note that this is a point-entry map, not a block-map.
     */
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const {  return mtx_->getRangeMap(); }

    //! \brief Computes the operator-multivector application.
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{trans}} \cdot X + \beta \cdot Y\f$. However, the details of operation
      vary according to the values of \c alpha and \c beta. Specifically
      - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
      - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
    */
    inline void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                      MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                      Teuchos::ETransp trans = Teuchos::NO_TRANS,
                      Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                      Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const {  mtx_->apply(X,Y,trans,alpha,beta); }

    //! Triangular Solve -- Matrix must be triangular.
    /*! Find X such that A*X = Y.
      Both \c X and \c Y are required to have constant stride.
    */
    inline void applyInverse(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Y,
                             MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                             Teuchos::ETransp trans) const {  mtx_->applyInverse(Y,X,trans); }

    //! Indicates whether this operator supports applying the adjoint operator.
    inline bool hasTransposeApply() const {  return mtx_->hasTransposeApply(); }

    //@}

    //! @name Attribute Query Methods
    //@{

    //! Returns the block-row map.
    inline const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRowMap() const {  return mtx_->getBlockRowMap(); }

    //! Returns the block-column map.
    inline const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockColMap() const {  return mtx_->getBlockColMap(); }

    //! Returns the block-domain map.
    inline const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockDomainMap() const {  return mtx_->getBlockDomainMap(); }

    //! Returns the block-range map.
    inline const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRangeMap() const {  return mtx_->getBlockRangeMap(); }

    //! Returns the point-row map.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointRowMap() const {  return mtx_->getPointRowMap(); }

    //! Returns the point-column map.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointColMap() const {  return mtx_->getPointColMap(); }

    //! Return true if fillComplete has been called, false otherwise.
    inline bool isFillComplete() const {  return mtx_->isFillComplete(); }
    //@}

    //! @name Insertion Methods
    //@{

    //! Set the specified scalar throughout the matrix.
    /*!
      This method may be called any time (before or after fillComplete()).
    */
    inline void putScalar(Scalar s) {  mtx_->putScalar(s); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, it will be
      over-written (replaced) by the input block-entry.

      This method may be called any time (before or after fillComplete()).
    */
    inline void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry) {  mtx_->setGlobalBlockEntry(globalBlockRow, globalBlockCol, blockEntry); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The coefficients of the specified block-entry will be
      over-written (replaced) by the input block-entry.
    */
    inline void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry) {  mtx_->setLocalBlockEntry(localBlockRow, localBlockCol, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, the contents of the
      input block-entry will be added to the values that are already present.

      This method may be called any time (before or after fillComplete()).
    */
    inline void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry) {  mtx_->sumIntoGlobalBlockEntry(globalBlockRow, globalBlockCol, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The contents of the input block-entry will be added to the values that are
      already present in the matrix.
    */
    inline void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry) {  mtx_->sumIntoLocalBlockEntry(localBlockRow, localBlockCol, blockEntry); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, it will be
      over-written (replaced) by the input block-entry.

      This method may be called any time (before or after fillComplete()).
    */
    inline void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) {  mtx_->setGlobalBlockEntry(globalBlockRow, globalBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The coefficients for the specified block-entry will be
      over-written (replaced) by the input block-entry.
    */
    inline void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) {  mtx_->setLocalBlockEntry(localBlockRow, localBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, the contents of the
      input block-entry will be added to the values that are already present.

      This method may be called any time (before or after fillComplete()).
    */
    inline void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) {  mtx_->sumIntoGlobalBlockEntry(globalBlockRow, globalBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The contents of the input block-entry will be added to the values that are
      already present in the matrix.
    */
    inline void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) {  mtx_->sumIntoLocalBlockEntry(localBlockRow, localBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //@}

    //! @name Transformational Methods
    //@{

    //! Transition the matrix to the packed, optimized-storage state.
    /*!
      This method also sets the domain and range maps.
    */
    inline void fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap) {  mtx_->fillComplete(blockDomainMap, blockRangeMap); }

    //! Transition the matrix to the packed, optimized-storage state.
    /*!
      This method internally calls fillComplete(getBlockRowMap(),getBlockRowMap()).
    */
    inline void fillComplete() {  mtx_->fillComplete(); }
    //@}

    //! @name Extraction Methods
    //@{

    //! Returns a const read-only view of a block-entry.
    /*!
      The arguments numPtRows and numPtCols are set to the dimensions of the block-
      entry on output.
      The stride (LDA in Blas terminology) is equal to numPtRows.

      This method may be called any time (before or after fillComplete()), but will
      throw an exception if the specified block-entry doesn't already exist.
    */
    inline void getGlobalBlockEntryView(GlobalOrdinal globalBlockRow,
                                        GlobalOrdinal globalBlockCol,
                                        LocalOrdinal& numPtRows,
                                        LocalOrdinal& numPtCols,
                                        Teuchos::ArrayRCP<const Scalar>& blockEntry) const {  mtx_->getGlobalBlockEntryView(globalBlockRow, globalBlockCol, numPtRows, numPtCols, blockEntry); }

    //! Returns a non-const read-write view of a block-entry.
    /*! Creates the block-entry if it doesn't already exist, and if:
      - the arguments numPtRows and numPtCols are set on entry (and nonzero),
      - and if fillComplete() has not yet been called.

      Important Note: Be very careful managing the lifetime of this view.
      If fillComplete() has been called, and if you are running on a GPU,
      this view may be a copy of memory from the GPU, and your changes to the
      view won't be copied back to the GPU until your ArrayRCP is destroyed
      or set to Teuchos::null.
    */
    inline void getGlobalBlockEntryViewNonConst(GlobalOrdinal globalBlockRow,
                                                GlobalOrdinal globalBlockCol,
                                                LocalOrdinal& numPtRows,
                                                LocalOrdinal& numPtCols,
                                                Teuchos::ArrayRCP<Scalar>& blockEntry) {  mtx_->getGlobalBlockEntryViewNonConst(globalBlockRow, globalBlockCol, numPtRows, numPtCols, blockEntry); }

    //! Returns a const read-only view of a block-entry.
    /*!
      The arguments numPtRows and numPtCols are set to the dimensions of the block-
      entry on output.
      The stride (LDA in Blas terminology) is equal to numPtRows.
      Throws an exception if fillComplete() has not yet been called, or if the
      specified block-entry doesn't exist.

      This method may only be called after fillComplete() has been called, and will
      throw an exception if the specified block-entry doesn't already exist.
    */
    inline void getLocalBlockEntryView(LocalOrdinal localBlockRow,
                                       LocalOrdinal localBlockCol,
                                       LocalOrdinal& numPtRows, 
                                       LocalOrdinal& numPtCols,
                                       Teuchos::ArrayRCP<const Scalar>& blockEntry) const {  mtx_->getLocalBlockEntryView(localBlockRow, localBlockCol, numPtRows, numPtCols, blockEntry); }

    //! Returns a non-const read-write view of a block-entry.
    /*!
      The arguments numPtRows and numPtCols are set to the dimensions of the block-
      entry on output.
      The stride (LDA in Blas terminology) is equal to numPtRows.
      Throws an exception if fillComplete() has not yet been called, or if the
      specified block-entry doesn't exist.

      Important Note: Be very careful managing the lifetime of this view.
      If fillComplete() has been called, and if you are running on a GPU,
      this view may be a copy of memory from the GPU, and your changes to the
      view won't be copied back to the GPU until your ArrayRCP is destroyed
      or set to Teuchos::null.

      This method may only be called after fillComplete() has been called, and will
      throw an exception if the specified block-entry doesn't already exist.
    */
    inline void getLocalBlockEntryViewNonConst(LocalOrdinal localBlockRow,
                                               LocalOrdinal localBlockCol,
                                               LocalOrdinal& numPtRows,
                                               LocalOrdinal& numPtCols,
                                               Teuchos::ArrayRCP<Scalar>& blockEntry) {  mtx_->getLocalBlockEntryViewNonConst(localBlockRow, localBlockCol, numPtRows, numPtCols, blockEntry); }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{
    inline std::string description() const {  return mtx_->description(); }

    /** \brief Print the object with some verbosity level to a FancyOStream object.
     */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  mtx_->describable(); }
    //@}

    RCP< const Tpetra::VbrMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > getTpetra_VbrMatrix() const {  return mtx_; }

  private:
  
    const RCP< const Tpetra::VbrMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > mtx_;
  
  };//class VbrMatrix
  
}//namespace Cthulhu

#define CTHULHU_TPETRAVBRMATRIX_SHORT
#endif //CTHULHU_VBRMATRIX_DECL_HPP
