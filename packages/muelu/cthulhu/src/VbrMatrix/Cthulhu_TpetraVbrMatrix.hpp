#ifndef CTHULHU_TPETRAVBRMATRIX_HPP
#define CTHULHU_TPETRAVBRMATRIX_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_TPETRA
#error This file should be included only if HAVE_CTHULHU_TPETRA is defined.
#endif

#include "Cthulhu_Debug.hpp"
#include "Cthulhu_VbrMatrix.hpp"

#include <Tpetra_VbrMatrix.hpp>

/** \file Cthulhu_TpetraVbrMatrix.hpp

  Declarations for the class Cthulhu::VbrMatrix.
*/
namespace Cthulhu {

//! \brief VbrMatrix: Variable block row matrix.
/**
The VbrMatrix class has two significant 'states', distinguished by whether or not
storage has been optimized (packed) or not.

When the matrix is in the non-optimized-storage state, internal data
storage is in a non-contiguous data-structure that allows for
convenient insertion of data.

When the matrix is in the optimized-storage state, internal data is stored in
contiguous (packed) arrays. When in this state, existing entries may be updated
and replaced, but no new entries (indices and/or coefficients) may be inserted.
In other words, the sparsity pattern or structure of the matrix may not be
changed.

Use of the matrix as an Operator (performing matrix-vector multiplication) is
only allowed when it is in the optimized-storage state.

VbrMatrix has two constructors, one which leaves the matrix in the optimized-
storage state, and another which leaves the matrix in the non-optimized-storage
state.

When the VbrMatrix is constructed in the non-optimized-storage state, (and then
filled using methods such as setGlobalBlockEntry etc.), it can then be transformed
to the optimized-storage state by calling the method fillComplete().

Once in the optimized-storage state, the VbrMatrix can not be returned to the
non-optimized-storage state.
*/
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

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Constructor specifying the row-map and the max number of (block) non-zeros for all rows.
    /*! After this constructor completes, the VbrMatrix is in the non-packed,
      non-optimized-storage, isFillComplete()==false state.
      Block-entries (rectangular, dense submatrices) may be inserted using class
      methods such as setGlobalBlockEntry(...), declared below.
    */
    // TODO VbrMatrix(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Constructor specifying a pre-filled block-graph.
    /*! Constructing a VbrMatrix with a pre-filled graph means that the matrix will
      start out in the optimized-storage state, i.e., isFillComplete()==true.
      The graph provided to this constructor must be already filled.
      (If blkGraph->isFillComplete() != true, an exception is thrown.)

      Entries in the input BlockCrsGraph correspond to block-entries in the
      VbrMatrix. In other words, the VbrMatrix will have a block-row corresponding
      to each row in the graph, and a block-entry corresponding to each column-
      index in the graph.
    */
    //TODO: need BlockCrsGraph
    VbrMatrix(const Teuchos::RCP<const BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& blkGraph);
#endif // CTHULHU_NOT_IMPLEMENTED

    TpetraVbrMatrix(const Teuchos::RCP<const Tpetra::VbrMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &mtx) : mtx_(mtx) { CTHULHU_DEBUG_ME; } //TODO

    //! Destructor
    virtual ~TpetraVbrMatrix();

    //@}

    //! @name Advanced Mathematical operations

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Multiply this matrix by a MultiVector.
    /*! \c X is required to be post-imported, i.e., described by the column map
      of the matrix. \c Y is required to be pre-exported, i.e., described by
      the row map of the matrix.
      See also the Operator::apply method which is implemented below.
    */
    //TODO virtual
    template <class DomainScalar, class RangeScalar>
    void multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const;
#endif // CTHULHU_NOT_IMPLEMENTED
    //@}

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Triangular Solve -- Matrix must be triangular.
    /*! Find X such that A*X = Y.
      \c X is required to be post-imported, i.e., described by the column map
      of the matrix. \c Y is required to be pre-exported, i.e., described by
      the row map of the matrix.

      Both \c X and \c Y are required to have constant stride.

      Note that if the diagonal block-entries are stored, they must be triangular.
      I.e., the matrix structure must be block-triangular, and any diagonal blocks
      must be "point"-triangular, meaning that coefficients on the "wrong" side of the
      point-diagonal must be zero.
    */
    //TODO virtual
    template <class DomainScalar, class RangeScalar>
    void solve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> & Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, Teuchos::ETransp trans) const;
#endif // CTHULHU_NOT_IMPLEMENTED
    //@}

    //! @name Operator Methods
    //@{

    //! Returns the Map associated with the domain of this operator.
    /*! Note that this is a point-entry map, not a block-map.
     */
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const { CTHULHU_DEBUG_ME; return mtx_->getDomainMap(); }

    //! Returns the Map associated with the range of this operator, which must be compatible with Y.getMap().
    /*! Note that this is a point-entry map, not a block-map.
     */
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const { CTHULHU_DEBUG_ME; return mtx_->getRangeMap(); }

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
                      Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const { CTHULHU_DEBUG_ME; mtx_->apply(X,Y,trans,alpha,beta); }

    //! Triangular Solve -- Matrix must be triangular.
    /*! Find X such that A*X = Y.
      Both \c X and \c Y are required to have constant stride.
    */
    inline void applyInverse(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Y,
                             MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                             Teuchos::ETransp trans) const { CTHULHU_DEBUG_ME; mtx_->applyInverse(Y,X,trans); }

    //! Indicates whether this operator supports applying the adjoint operator.
    inline bool hasTransposeApply() const { CTHULHU_DEBUG_ME; return mtx_->hasTransposeApply(); }

    //@}

    //! @name Attribute Query Methods
    //@{

    //! Returns the block-row map.
    inline const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRowMap() const { CTHULHU_DEBUG_ME; return mtx_->getBlockRowMap(); }

    //! Returns the block-column map.
    inline const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockColMap() const { CTHULHU_DEBUG_ME; return mtx_->getBlockColMap(); }

    //! Returns the block-domain map.
    inline const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockDomainMap() const { CTHULHU_DEBUG_ME; return mtx_->getBlockDomainMap(); }

    //! Returns the block-range map.
    inline const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRangeMap() const { CTHULHU_DEBUG_ME; return mtx_->getBlockRangeMap(); }

    //! Returns the point-row map.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointRowMap() const { CTHULHU_DEBUG_ME; return mtx_->getPointRowMap(); }

    //! Returns the point-column map.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointColMap() const { CTHULHU_DEBUG_ME; return mtx_->getPointColMap(); }

    //! Return true if fillComplete has been called, false otherwise.
    inline bool isFillComplete() const { CTHULHU_DEBUG_ME; return mtx_->isFillComplete(); }
    //@}

    //! @name Insertion Methods
    //@{

    //! Set the specified scalar throughout the matrix.
    /*!
      This method may be called any time (before or after fillComplete()).
    */
    inline void putScalar(Scalar s) { CTHULHU_DEBUG_ME; mtx_->putScalar(s); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, it will be
      over-written (replaced) by the input block-entry.

      This method may be called any time (before or after fillComplete()).
    */
    inline void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->setGlobalBlockEntry(globalBlockRow, globalBlockCol, blockEntry); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The coefficients of the specified block-entry will be
      over-written (replaced) by the input block-entry.
    */
    inline void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->setLocalBlockEntry(localBlockRow, localBlockCol, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, the contents of the
      input block-entry will be added to the values that are already present.

      This method may be called any time (before or after fillComplete()).
    */
    inline void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->sumIntoGlobalBlockEntry(globalBlockRow, globalBlockCol, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The contents of the input block-entry will be added to the values that are
      already present in the matrix.
    */
    inline void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->sumIntoLocalBlockEntry(localBlockRow, localBlockCol, blockEntry); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, it will be
      over-written (replaced) by the input block-entry.

      This method may be called any time (before or after fillComplete()).
    */
    inline void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->setGlobalBlockEntry(globalBlockRow, globalBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The coefficients for the specified block-entry will be
      over-written (replaced) by the input block-entry.
    */
    inline void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->setLocalBlockEntry(localBlockRow, localBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, the contents of the
      input block-entry will be added to the values that are already present.

      This method may be called any time (before or after fillComplete()).
    */
    inline void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->sumIntoGlobalBlockEntry(globalBlockRow, globalBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The contents of the input block-entry will be added to the values that are
      already present in the matrix.
    */
    inline void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->sumIntoLocalBlockEntry(localBlockRow, localBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //@}

    //! @name Transformational Methods
    //@{

    //! Transition the matrix to the packed, optimized-storage state.
    /*!
      This method also sets the domain and range maps.
    */
    inline void fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap) { CTHULHU_DEBUG_ME; mtx_->fillComplete(blockDomainMap, blockRangeMap); }

    //! Transition the matrix to the packed, optimized-storage state.
    /*!
      This method internally calls fillComplete(getBlockRowMap(),getBlockRowMap()).
    */
    inline void fillComplete() { CTHULHU_DEBUG_ME; mtx_->fillComplete(); }
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
                                        Teuchos::ArrayRCP<const Scalar>& blockEntry) const { CTHULHU_DEBUG_ME; mtx_->getGlobalBlockEntryView(globalBlockRow, globalBlockCol, numPtRows, numPtCols, blockEntry); }

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
                                                Teuchos::ArrayRCP<Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->getGlobalBlockEntryViewNonConst(globalBlockRow, globalBlockCol, numPtRows, numPtCols, blockEntry); }

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
                                       Teuchos::ArrayRCP<const Scalar>& blockEntry) const { CTHULHU_DEBUG_ME; mtx_->getLocalBlockEntryView(localBlockRow, localBlockCol, numPtRows, numPtCols, blockEntry); }

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
                                               Teuchos::ArrayRCP<Scalar>& blockEntry) { CTHULHU_DEBUG_ME; mtx_->getLocalBlockEntryViewNonConst(localBlockRow, localBlockCol, numPtRows, numPtCols, blockEntry); }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Return a copy of the (point-entry) diagonal values.
    /*!
      Throws an exception if the input-vector's map is not the same as
      getBlockRowMap()->getPointMap().
    */
    //TODO:Vector  
    inline void getLocalDiagCopy(Cthulhu::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& diag) const { CTHULHU_DEBUG_ME; mtx_->getLocalDiagCopy(diag); }
#endif // CTHULHU_NOT_IMPLEMENTED
    //@}

    //! @name Overridden from Teuchos::Describable
    //@{
    inline std::string description() const { CTHULHU_DEBUG_ME; return mtx_->description(); }

    /** \brief Print the object with some verbosity level to a FancyOStream object.
     */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { CTHULHU_DEBUG_ME; mtx_->describable(); }
    //@}

    RCP< const Tpetra::VbrMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > getTpetra_VbrMatrix() const { CTHULHU_DEBUG_ME; return mtx_; }

  private:
  
    const RCP< const Tpetra::VbrMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > mtx_;
  
  };//class VbrMatrix
  
}//namespace Cthulhu

//----------------------------------------------------------------------------
// Description of arrays representing the VBR format:
//
// (For more description, see this URL (valid as of 5/26/2010):
// http://docs.sun.com/source/817-0086-10/prog-sparse-support.html)
// ...and of course more can be found using google...
// The old Aztec manual was a great resource for this but I can't
// find a copy of that these days...
//
//
// Here is a brief description of the 6 arrays that are required to
// represent a VBR matrix in packed (contiguous-memory-storage) format:
//
// rptr: length num_block_rows + 1
//       rptr[i]: the pt-row corresponding to the i-th block-row
//       Note: rptr is getBlockRowMap()->getNodeFirstPointInBlocks().
//
// cptr: length num_distinct_block_cols + 1
//       cptr[j]: the pt-col corresponding to the j-th block-col
//       Note: cptr is getBlockColMap()->getNodeFirstPointInBlocks().
//
// bptr: length num_block_rows + 1
//       bptr[i]: location in bindx of the first nonzero block-entry
//                of the i-th block-row
//       Note: bptr is blkGraph_->getNodeRowOffsets();
//
// bindx: length num-nonzero-block-entries
//        bindx[j]: block-col-index of j-th block-entry
//        Note: bindx is blkGraph_->getNodePackedIndices();
//
// indx: length num-nonzero-block-entries + 1
//       indx[j]: location in vals of the beginning of the j-th
//       block-entry
//
// vals: length num-nonzero-scalar-entries
//
//
// Some example look-ups:
//
// int nbr = num_block_rows;
// int total_num_block_nonzeros = bptr[nbr];
// int total_num_scalar_nonzeros = indx[num_block_nonzeros];
// 
// //get arrays for i-th block-row:
// int* bindx_i = &bindx[bptr[i]];
// double* vals_i = &val[indx[bptr[i]]];
// int num_block_nonzeros_in_row_i = bptr[i+1]-bptr[i];
// 
//----------------------------------------------------------------------------

#define CTHULHU_TPETRAVBRMATRIX_SHORT
#endif //CTHULHU_VBRMATRIX_DECL_HPP
