#ifndef CTHULHU_VBRMATRIX_HPP
#define CTHULHU_VBRMATRIX_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_VbrMatrix.hpp>

#include "Cthulhu_ConfigDefs.hpp"
// #include "Cthulhu_Operator.hpp"
#include "Cthulhu_BlockMap.hpp"
#include "Cthulhu_MultiVector.hpp"
// #include "Cthulhu_BlockCrsGraph.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>

#include "Cthulhu_Debug.hpp"

/** \file Cthulhu_VbrMatrix.hpp

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
class VbrMatrix { //: public Cthulhu::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
 public:

  //! @name Constructor/Destructor Methods
  //@{

  //! Destructor
  virtual ~VbrMatrix();

  //@}

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

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
  virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const =0;

  //! Returns the Map associated with the range of this operator, which must be compatible with Y.getMap().
  /*! Note that this is a point-entry map, not a block-map.
  */
  virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const =0;

  //! \brief Computes the operator-multivector application.
  /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{trans}} \cdot X + \beta \cdot Y\f$. However, the details of operation
      vary according to the values of \c alpha and \c beta. Specifically
      - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
      - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
   */
  virtual void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                    Teuchos::ETransp trans = Teuchos::NO_TRANS,
                    Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                    Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const =0;
  
  //! Triangular Solve -- Matrix must be triangular.
  /*! Find X such that A*X = Y.
      Both \c X and \c Y are required to have constant stride.
  */
  virtual void applyInverse(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Y,
                   MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                   Teuchos::ETransp trans) const =0;

  //! Indicates whether this operator supports applying the adjoint operator.
  virtual bool hasTransposeApply() const =0;

  //@}

  //! @name Attribute Query Methods
  //@{

  //! Returns the block-row map.
  virtual const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRowMap() const =0;

  //! Returns the block-column map.
  virtual const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockColMap() const =0;

  //! Returns the block-domain map.
  virtual const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockDomainMap() const =0;

  //! Returns the block-range map.
  virtual const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRangeMap() const =0;

  //! Returns the point-row map.
  virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointRowMap() const =0;

  //! Returns the point-column map.
  virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointColMap() const =0;

  //! Return true if fillComplete has been called, false otherwise.
  virtual bool isFillComplete() const =0;
  //@}

  //! @name Insertion Methods
  //@{

  //! Set the specified scalar throughout the matrix.
  /*!
    This method may be called any time (before or after fillComplete()).
  */
  virtual void putScalar(Scalar s) =0;

  //!Copy the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, it will be
    over-written (replaced) by the input block-entry.

    This method may be called any time (before or after fillComplete()).
  */
  virtual void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry) =0;

  //!Copy the contents of the input block-entry into the matrix.
  /*!
    This method will throw an exception if fillComplete() has not yet been called,
    or if the specified block-entry doesn't already exist in the matrix.

    The coefficients of the specified block-entry will be
    over-written (replaced) by the input block-entry.
  */
  virtual void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry) =0;

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, the contents of the
    input block-entry will be added to the values that are already present.

    This method may be called any time (before or after fillComplete()).
  */
  virtual void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry) =0;

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will throw an exception if fillComplete() has not yet been called,
    or if the specified block-entry doesn't already exist in the matrix.

    The contents of the input block-entry will be added to the values that are
    already present in the matrix.
  */
  virtual void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry) =0;

  //!Copy the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, it will be
    over-written (replaced) by the input block-entry.

    This method may be called any time (before or after fillComplete()).
  */
  virtual void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) =0;

  //!Copy the contents of the input block-entry into the matrix.
  /*!
    This method will throw an exception if fillComplete() has not yet been called,
    or if the specified block-entry doesn't already exist in the matrix.

    The coefficients for the specified block-entry will be
    over-written (replaced) by the input block-entry.
  */
  virtual void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) =0;

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, the contents of the
    input block-entry will be added to the values that are already present.

    This method may be called any time (before or after fillComplete()).
  */
  virtual void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) =0;

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will throw an exception if fillComplete() has not yet been called,
    or if the specified block-entry doesn't already exist in the matrix.

    The contents of the input block-entry will be added to the values that are
    already present in the matrix.
  */
  virtual void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry) =0;

  //@}

  //! @name Transformational Methods
  //@{

  //! Transition the matrix to the packed, optimized-storage state.
  /*!
    This method also sets the domain and range maps.
  */
  virtual void fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap) =0;

  //! Transition the matrix to the packed, optimized-storage state.
  /*!
    This method internally calls fillComplete(getBlockRowMap(),getBlockRowMap()).
  */
  virtual void fillComplete() =0;
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
  virtual void getGlobalBlockEntryView(GlobalOrdinal globalBlockRow,
                              GlobalOrdinal globalBlockCol,
                              LocalOrdinal& numPtRows,
                              LocalOrdinal& numPtCols,
                              Teuchos::ArrayRCP<const Scalar>& blockEntry) const =0;

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
  virtual void getGlobalBlockEntryViewNonConst(GlobalOrdinal globalBlockRow,
                                      GlobalOrdinal globalBlockCol,
                                      LocalOrdinal& numPtRows,
                                      LocalOrdinal& numPtCols,
                                      Teuchos::ArrayRCP<Scalar>& blockEntry) =0;

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
  virtual void getLocalBlockEntryView(LocalOrdinal localBlockRow,
                             LocalOrdinal localBlockCol,
                             LocalOrdinal& numPtRows, 
                             LocalOrdinal& numPtCols,
                             Teuchos::ArrayRCP<const Scalar>& blockEntry) const =0;

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
  virtual void getLocalBlockEntryViewNonConst(LocalOrdinal localBlockRow,
                                     LocalOrdinal localBlockCol,
                                     LocalOrdinal& numPtRows,
                                     LocalOrdinal& numPtCols,
                                     Teuchos::ArrayRCP<Scalar>& blockEntry) =0;

#ifdef CTHULHU_NOT_IMPLEMENTED
  //! Return a copy of the (point-entry) diagonal values.
  /*!
    Throws an exception if the input-vector's map is not the same as
    getBlockRowMap()->getPointMap().
  */
  //TODO: need Vector  
virtual void getLocalDiagCopy(Cthulhu::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& diag) const =0;
#endif // CTHULHU_NOT_IMPLEMENTED
  //@}

  //! @name Overridden from Teuchos::Describable
  //@{
  virtual std::string description() const =0;

  /** \brief Print the object with some verbosity level to a FancyOStream object.
  */
  virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0;
  //@}

#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

};//class VbrMatrix

}//namespace Cthulhu

#define CTHULHU_VBRMATRIX_SHORT
#endif //CTHULHU_VBRMATRIX_DECL_HPP
