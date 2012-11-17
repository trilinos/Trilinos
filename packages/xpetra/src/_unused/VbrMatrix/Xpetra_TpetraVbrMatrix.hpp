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
#ifndef XPETRA_TPETRAVBRMATRIX_HPP
#define XPETRA_TPETRAVBRMATRIX_HPP

#include "Xpetra_ConfigDefs.hpp"

#ifndef HAVE_XPETRA_TPETRA
#error This file should be included only if HAVE_XPETRA_TPETRA is defined.
#endif

#include "Xpetra_VbrMatrix.hpp"

#include <Tpetra_VbrMatrix.hpp>

namespace Xpetra {

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

    //! @name Matrix Methods
    //@{

    //! Returns the Map associated with the domain of this operator.
    /*! Note that this is a point-entry map, not a block-map.
     */
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const {  return mtx_->getDomainMap(); }

    //! Returns the Map associated with the range of this operator, which must be compatible with Y.getMap().
    /*! Note that this is a point-entry map, not a block-map.
     */
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const {  return mtx_->getRangeMap(); }

    //! \brief Computes the operator-multivector application.
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{trans}} \cdot X + \beta \cdot Y\f$. However, the details of operation
      vary according to the values of \c alpha and \c beta. Specifically
      - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
      - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
    */
    void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                      MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                      Teuchos::ETransp trans = Teuchos::NO_TRANS,
                      Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                      Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const {  mtx_->apply(X,Y,trans,alpha,beta); }

    //! Triangular Solve -- Matrix must be triangular.
    /*! Find X such that A*X = Y.
      Both \c X and \c Y are required to have constant stride.
    */
    void applyInverse(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Y,
                             MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                             Teuchos::ETransp trans) const {  mtx_->applyInverse(Y,X,trans); }

    //! Indicates whether this operator supports applying the adjoint operator.
    bool hasTransposeApply() const {  return mtx_->hasTransposeApply(); }

    //@}

    //! @name Attribute Query Methods
    //@{

    //! Returns the block-row map.
    const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRowMap() const {  return mtx_->getBlockRowMap(); }

    //! Returns the block-column map.
    const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockColMap() const {  return mtx_->getBlockColMap(); }

    //! Returns the block-domain map.
    const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockDomainMap() const {  return mtx_->getBlockDomainMap(); }

    //! Returns the block-range map.
    const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRangeMap() const {  return mtx_->getBlockRangeMap(); }

    //! Returns the point-row map.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointRowMap() const {  return mtx_->getPointRowMap(); }

    //! Returns the point-column map.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointColMap() const {  return mtx_->getPointColMap(); }

    //! Return true if fillComplete has been called, false otherwise.
    bool isFillComplete() const {  return mtx_->isFillComplete(); }
    //@}

    //! @name Insertion Methods
    //@{

    //! Set the specified scalar throughout the matrix.
    /*!
      This method may be called any time (before or after fillComplete()).
    */
    void putScalar(Scalar s) {  mtx_->putScalar(s); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, it will be
      over-written (replaced) by the input block-entry.

      This method may be called any time (before or after fillComplete()).
    */
    void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry) {  mtx_->setGlobalBlockEntry(globalBlockRow, globalBlockCol, blockEntry); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The coefficients of the specified block-entry will be
      over-written (replaced) by the input block-entry.
    */
    void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry) {  mtx_->setLocalBlockEntry(localBlockRow, localBlockCol, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, the contents of the
      input block-entry will be added to the values that are already present.

      This method may be called any time (before or after fillComplete()).
    */
    void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry) {  mtx_->sumIntoGlobalBlockEntry(globalBlockRow, globalBlockCol, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The contents of the input block-entry will be added to the values that are
      already present in the matrix.
    */
    void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry) {  mtx_->sumIntoLocalBlockEntry(localBlockRow, localBlockCol, blockEntry); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, it will be
      over-written (replaced) by the input block-entry.

      This method may be called any time (before or after fillComplete()).
    */
    void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) {  mtx_->setGlobalBlockEntry(globalBlockRow, globalBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //!Copy the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The coefficients for the specified block-entry will be
      over-written (replaced) by the input block-entry.
    */
    void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) {  mtx_->setLocalBlockEntry(localBlockRow, localBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will create the specified block-entry if it doesn't already exist,
      but only if fillComplete() has not yet been called.

      If the specified block-entry already exists in the matrix, the contents of the
      input block-entry will be added to the values that are already present.

      This method may be called any time (before or after fillComplete()).
    */
    void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) {  mtx_->sumIntoGlobalBlockEntry(globalBlockRow, globalBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //!Add the contents of the input block-entry into the matrix.
    /*!
      This method will throw an exception if fillComplete() has not yet been called,
      or if the specified block-entry doesn't already exist in the matrix.

      The contents of the input block-entry will be added to the values that are
      already present in the matrix.
    */
    void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) {  mtx_->sumIntoLocalBlockEntry(localBlockRow, localBlockCol, blkRowSize, blkColSize, LDA, blockEntry); }

    //@}

    //! @name Transformational Methods
    //@{

    //! Transition the matrix to the packed, optimized-storage state.
    /*!
      This method also sets the domain and range maps.
    */
    void fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap) {  mtx_->fillComplete(blockDomainMap, blockRangeMap); }

    //! Transition the matrix to the packed, optimized-storage state.
    /*!
      This method internally calls fillComplete(getBlockRowMap(),getBlockRowMap()).
    */
    void fillComplete() {  mtx_->fillComplete(); }
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
    void getGlobalBlockEntryView(GlobalOrdinal globalBlockRow,
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
    void getGlobalBlockEntryViewNonConst(GlobalOrdinal globalBlockRow,
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
    void getLocalBlockEntryView(LocalOrdinal localBlockRow,
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
    void getLocalBlockEntryViewNonConst(LocalOrdinal localBlockRow,
                                               LocalOrdinal localBlockCol,
                                               LocalOrdinal& numPtRows,
                                               LocalOrdinal& numPtCols,
                                               Teuchos::ArrayRCP<Scalar>& blockEntry) {  mtx_->getLocalBlockEntryViewNonConst(localBlockRow, localBlockCol, numPtRows, numPtCols, blockEntry); }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{
    std::string description() const {  return mtx_->description(); }

    /** \brief Print the object with some verbosity level to a FancyOStream object.
     */
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  mtx_->describable(); }
    //@}

    RCP< const Tpetra::VbrMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > getTpetra_VbrMatrix() const {  return mtx_; }

  private:

    const RCP< const Tpetra::VbrMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > mtx_;

  };//class VbrMatrix

}//namespace Xpetra

#define XPETRA_TPETRAVBRMATRIX_SHORT
#endif //XPETRA_VBRMATRIX_DECL_HPP
