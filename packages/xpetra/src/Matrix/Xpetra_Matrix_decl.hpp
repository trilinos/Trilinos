// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MATRIX_DECL_HPP
#define XPETRA_MATRIX_DECL_HPP

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"
#include "Xpetra_MatrixView.hpp"
#include "Xpetra_Operator.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

/** \file Xpetra_Matrix.hpp

Declarations for the class Xpetra::Matrix.
*/
namespace Xpetra {

/*!
 @class Xpetra::Matrix class.
 @brief Xpetra-specific matrix class.

 This class is specific to Xpetra and has no analogue in Epetra or Tpetra.  The main motivation for this class is to be able to access matrix data in a manner different than how it is stored.
 For example, it might be more convenient to treat ("view") a matrix stored in compressed row storage as if it were a block matrix.  The Xpetra::Matrix class is intended to manage these "views".

 <B>How to create a Matrix from an existing CrsMatrix</B>

 @code
 RCP<Xpetra::CrsMatrix> crsA;
 RCP<Xpetra::Matrix>    A  = rcp(new CrsMatrixWrap(crsA));
 @endcode

*/

typedef std::string viewLabel_t;

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class Matrix : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
  typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> CrsGraph;
  typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixFactory;
  typedef Xpetra::MatrixView<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatrixView;

 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

#ifdef HAVE_XPETRA_TPETRA
  typedef typename CrsMatrix::local_matrix_type local_matrix_type;
#endif

  //! @name Constructor/Destructor Methods
  //@{

  Matrix();

  //! Destructor
  virtual ~Matrix();

  //@}

  //! @name View management methods
  //@{
  void CreateView(viewLabel_t viewLabel, const RCP<const Map> &rowMap, const RCP<const Map> &colMap);

  // JG TODO: why this is a member function??
  void CreateView(const viewLabel_t viewLabel, const RCP<const Matrix> &A, bool transposeA = false, const RCP<const Matrix> &B = Teuchos::null, bool transposeB = false);

  //! Print all of the views associated with the Matrix.
  void PrintViews(Teuchos::FancyOStream &out) const;

  void RemoveView(const viewLabel_t viewLabel);

  const viewLabel_t SwitchToView(const viewLabel_t viewLabel);

  bool IsView(const viewLabel_t viewLabel) const;

  const viewLabel_t SwitchToDefaultView();

  const viewLabel_t &GetDefaultViewLabel() const;

  const viewLabel_t &GetCurrentViewLabel() const;

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
  virtual void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) = 0;

  //! Insert matrix entries, using local IDs.
  /** All index values must be in the local space.
      \pre \c localRow exists as an ID in the local row map
      \pre <tt>isGloballyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isLocallyIndexed() == true</tt>
  */
  virtual void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) = 0;

  //! \brief Replace matrix entries, using global IDs.
  /** All index values must be in the global space.

      \pre \c globalRow is a global row belonging to the matrix on this node.

      \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
  virtual void replaceGlobalValues(GlobalOrdinal globalRow,
                                   const ArrayView<const GlobalOrdinal> &cols,
                                   const ArrayView<const Scalar> &vals) = 0;

  //! Replace matrix entries, using local IDs.
  /** All index values must be in the local space.
      Note that if a value is not already present for the specified location in the matrix, the input value will be ignored silently.
   */
  virtual void replaceLocalValues(LocalOrdinal localRow,
                                  const ArrayView<const LocalOrdinal> &cols,
                                  const ArrayView<const Scalar> &vals) = 0;

  //! Set all matrix entries equal to scalar
  virtual void setAllToScalar(const Scalar &alpha) = 0;

  //! Scale the current values of a matrix, this = alpha*this.
  virtual void scale(const Scalar &alpha) = 0;

  //@}

  //! @name Transformational Methods
  //@{

  /*! Resume fill operations.
    After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

    resumeFill() may be called repeatedly.

    \post  <tt>isFillActive() == true<tt>
    \post  <tt>isFillComplete() == false<tt>
  */
  virtual void resumeFill(const RCP<ParameterList> &params = null) = 0;

  /*! \brief Signal that data entry is complete, specifying domain and range maps.

  Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

  \pre  <tt>isFillActive() == true<tt>
  \pre <tt>isFillComplete()() == false<tt>

  \post <tt>isFillActive() == false<tt>
  \post <tt>isFillComplete() == true<tt>
  \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
  */
  virtual void fillComplete(const RCP<const Map> &domainMap, const RCP<const Map> &rangeMap, const RCP<ParameterList> &params = null) = 0;

  /*! \brief Signal that data entry is complete.

  Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

  \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

  \pre  <tt>isFillActive() == true<tt>
  \pre <tt>isFillComplete()() == false<tt>

  \post <tt>isFillActive() == false<tt>
  \post <tt>isFillComplete() == true<tt>
  \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
  */
  // TODO : Get ride of "Tpetra"::OptimizeOption
  virtual void fillComplete(const RCP<ParameterList> &params = null) = 0;

  //@}

  //! @name Methods implementing RowMatrix
  //@{

  //! Returns the Map that describes the row distribution in this matrix.
  virtual const RCP<const Map> &getRowMap() const;

  //! Returns the Map that describes the row distribution in this matrix.
  virtual const RCP<const Map> &getRowMap(viewLabel_t viewLabel) const;

  //! \brief Returns the Map that describes the column distribution in this matrix.
  //! This might be <tt>null</tt> until fillComplete() is called.
  virtual const RCP<const Map> &getColMap() const;

  //! \brief Returns the Map that describes the column distribution in this matrix.
  virtual const RCP<const Map> &getColMap(viewLabel_t viewLabel) const;

  //! Returns the number of global rows in this matrix.
  /** Undefined if isFillActive().
   */
  virtual global_size_t getGlobalNumRows() const = 0;

  //! \brief Returns the number of global columns in the matrix.
  /** Undefined if isFillActive().
   */
  virtual global_size_t getGlobalNumCols() const = 0;

  //! Returns the number of matrix rows owned on the calling node.
  virtual size_t getLocalNumRows() const = 0;

  //! Returns the global number of entries in this matrix.
  virtual global_size_t getGlobalNumEntries() const = 0;

  //! Returns the local number of entries in this matrix.
  virtual size_t getLocalNumEntries() const = 0;

  //! Returns the current number of entries on this node in the specified local row.
  /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
  virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

  //! Returns the current number of entries in the specified global row.
  /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row is not owned by this process. */
  virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  /** Undefined if isFillActive().
   */
  virtual size_t getGlobalMaxNumRowEntries() const = 0;

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  /** Undefined if isFillActive().
   */
  virtual size_t getLocalMaxNumRowEntries() const = 0;

  //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
  virtual bool isLocallyIndexed() const = 0;

  //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
  virtual bool isGloballyIndexed() const = 0;

  //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
  virtual bool isFillComplete() const = 0;

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
                               size_t &NumEntries) const = 0;

  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*!
    \param GlobalRow - (In) Global row number for which indices are desired.
    \param Indices   - (Out) Global column indices corresponding to values.
    \param Values    - (Out) Row values
    \pre <tt>isLocallyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const = 0;

  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices  - (Out) Local column indices corresponding to values.
    \param Values   - (Out) Row values
    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

    Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const = 0;

  //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing
    the zero and non-zero diagonals owned by this node. */
  virtual void getLocalDiagCopy(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) const = 0;

  //! Get Frobenius norm of the matrix
  virtual typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const = 0;

  //! Left scale matrix using the given vector entries
  virtual void leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &x) = 0;

  //! Right scale matrix using the given vector entries
  virtual void rightScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &x) = 0;

  //! Returns true if globalConstants have been computed; false otherwise
  virtual bool haveGlobalConstants() const = 0;

  //@}

  //! Implements DistObject interface
  //{@

  //! Access function for the Tpetra::Map this DistObject was constructed with.
  virtual const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getMap() const = 0;

  // TODO: first argument of doImport/doExport should be a Xpetra::DistObject

  //! Import.
  virtual void doImport(const Matrix &source,
                        const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) = 0;

  //! Export.
  virtual void doExport(const Matrix &dest,
                        const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) = 0;

  //! Import (using an Exporter).
  virtual void doImport(const Matrix &source,
                        const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) = 0;

  //! Export (using an Importer).
  virtual void doExport(const Matrix &dest,
                        const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) = 0;

  // @}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! @name Overridden from Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  virtual std::string description() const = 0;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const = 0;
  //@}

  //! @name Overridden from Teuchos::LabeledObject
  //@{
  virtual void setObjectLabel(const std::string &objectLabel) = 0;
  //@}

  // JG: Added:

  //! Supports the getCrsGraph() call
  virtual bool hasCrsGraph() const = 0;

  //! Returns the CrsGraph associated with this matrix.
  virtual RCP<const CrsGraph> getCrsGraph() const = 0;

  // To keep the normal virtual matrix-multivector definition of apply before overloading with the region variant
  using Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply;

  //! Computes the matrix-multivector multiplication for region layout matrices
  virtual void apply(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
                     MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
                     Teuchos::ETransp mode,
                     Scalar alpha,
                     Scalar beta,
                     bool sumInterfaceValues,
                     const RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > &regionInterfaceImporter,
                     const Teuchos::ArrayRCP<LocalOrdinal> &regionInterfaceLIDs) const = 0;

  // ----------------------------------------------------------------------------------
  // "TEMPORARY" VIEW MECHANISM
  /**
   * Set fixed block size of operator (e.g., 3 for 3 DOFs per node).
   *
   * @param blksize: block size denoting how many DOFs per node are used (LocalOrdinal)
   * @param offset:  global offset allows to define operators with global indices starting from a given value "offset" instead of 0. (GlobalOrdinal, default = 0)
   * */
  void SetFixedBlockSize(LocalOrdinal blksize, GlobalOrdinal offset = 0);

  //==========================================================================

  LocalOrdinal GetFixedBlockSize() const;  // TODO: why LocalOrdinal?

  //! Returns true, if `SetFixedBlockSize` has been called before.
  bool IsFixedBlockSizeSet() const;

  //! Returns the block size of the storage mechanism, which is usually 1, except for Tpetra::BlockCrsMatrix
  virtual LocalOrdinal GetStorageBlockSize() const = 0;

  // ----------------------------------------------------------------------------------

  virtual void SetMaxEigenvalueEstimate(Scalar const &sigma);

  // ----------------------------------------------------------------------------------

  virtual Scalar GetMaxEigenvalueEstimate() const;

  // ----------------------------------------------------------------------------------
#ifdef HAVE_XPETRA_TPETRA
  virtual local_matrix_type getLocalMatrixDevice() const                    = 0;
  virtual typename local_matrix_type::HostMirror getLocalMatrixHost() const = 0;
#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif
  // ----------------------------------------------------------------------------------

  //! Compute a residual R = B - (*this) * X
  virtual void residual(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
                        const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B,
                        MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &R) const = 0;

 protected:
  Teuchos::Hashtable<viewLabel_t, RCP<MatrixView> > operatorViewTable_;  // hashtable storing the operator views (keys = view names, values = views).

  viewLabel_t defaultViewLabel_;  // label of the view associated with inital Matrix construction
  viewLabel_t currentViewLabel_;  // label of the current view

};  // class Matrix

}  // namespace Xpetra

#define XPETRA_MATRIX_SHORT
#endif  // XPETRA_MATRIX_DECL_HPP
