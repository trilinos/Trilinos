// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_CRSMATRIXWRAP_DECL_HPP
#define XPETRA_CRSMATRIXWRAP_DECL_HPP

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MultiVector_decl.hpp"
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
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class CrsMatrixWrap : public Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> CrsGraph;
#ifdef HAVE_XPETRA_TPETRA
  typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraCrsMatrix;
#endif
  typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixFactory;
  typedef Xpetra::MatrixView<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatrixView;
#ifdef HAVE_XPETRA_TPETRA
  typedef typename CrsMatrix::local_matrix_type local_matrix_type;
#endif

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor for a dynamic profile matrix (Epetra only)
  CrsMatrixWrap(const RCP<const Map> &rowMap);

  //! Constructor specifying fixed number of entries for each row.
  CrsMatrixWrap(const RCP<const Map> &rowMap,
                size_t maxNumEntriesPerRow);

  //! Constructor specifying (possibly different) number of entries in each row.
  CrsMatrixWrap(const RCP<const Map> &rowMap,
                const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc);

  //! Constructor specifying fixed number of entries for each row and column map
  CrsMatrixWrap(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, size_t maxNumEntriesPerRow);

  //! Constructor specifying fixed number of entries for each row and column map
  CrsMatrixWrap(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc);

#ifdef HAVE_XPETRA_TPETRA
  //! Constructor specifying fixed number of entries for each row and column map
  CrsMatrixWrap(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, const local_matrix_type &lclMatrix, const Teuchos::RCP<Teuchos::ParameterList> &params = null);

  //! Constructor specifying fixed number of entries for each row and column map
  CrsMatrixWrap(const local_matrix_type &lclMatrix, const RCP<const Map> &rowMap, const RCP<const Map> &colMap,
                const RCP<const Map> &domainMap = Teuchos::null, const RCP<const Map> &rangeMap = Teuchos::null,
                const Teuchos::RCP<Teuchos::ParameterList> &params = null);
#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif

  CrsMatrixWrap(RCP<CrsMatrix> matrix);

  CrsMatrixWrap(const RCP<const CrsGraph> &graph, const RCP<ParameterList> &paramList = Teuchos::null);

  CrsMatrixWrap(const RCP<const CrsGraph> &graph,
                typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::values_type &values,
                const RCP<ParameterList> &paramList = Teuchos::null);

  //! Destructor
  virtual ~CrsMatrixWrap();

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
  void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals);

  //! Insert matrix entries, using local IDs.
  /** All index values must be in the local space.
      \pre \c localRow exists as an ID in the global row map
      \pre <tt>isGloballyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isLocallyIndexed() == true</tt>
  */
  void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals);

  //! \brief Replace matrix entries, using global IDs.
  /** All index values must be in the global space.

  \pre \c globalRow is a global row belonging to the matrix on this node.

  \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
  void replaceGlobalValues(GlobalOrdinal globalRow,
                           const ArrayView<const GlobalOrdinal> &cols,
                           const ArrayView<const Scalar> &vals);

  //! Replace matrix entries, using local IDs.
  /** All index values must be in the local space.
      Note that if a value is not already present for the specified location in the matrix, the input value will be ignored silently.
  */
  void replaceLocalValues(LocalOrdinal localRow,
                          const ArrayView<const LocalOrdinal> &cols,
                          const ArrayView<const Scalar> &vals);

  //! Set all matrix entries equal to scalar
  virtual void setAllToScalar(const Scalar &alpha);

  //! Scale the current values of a matrix, this = alpha*this.
  void scale(const Scalar &alpha);

  //@}

  //! @name Transformational Methods
  //@{

  /*! Resume fill operations.
    After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

    resumeFill() may be called repeatedly.

    \post  <tt>isFillActive() == true<tt>
    \post  <tt>isFillComplete() == false<tt>
  */
  void resumeFill(const RCP<ParameterList> &params = null);

  /*! \brief Signal that data entry is complete, specifying domain and range maps.

  Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

  \pre  <tt>isFillActive() == true<tt>
  \pre <tt>isFillComplete()() == false<tt>

  \post <tt>isFillActive() == false<tt>
  \post <tt>isFillComplete() == true<tt>
  \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
  */
  void fillComplete(const RCP<const Map> &domainMap, const RCP<const Map> &rangeMap, const RCP<Teuchos::ParameterList> &params = null);

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
  void fillComplete(const RCP<ParameterList> &params = null);

  //@}

  //! Returns the number of global rows in this matrix.
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumRows() const;

  //! \brief Returns the number of global columns in the matrix.
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumCols() const;

  //! Returns the number of matrix rows owned on the calling node.
  size_t getLocalNumRows() const;

  //! Returns the global number of entries in this matrix.
  global_size_t getGlobalNumEntries() const;

  //! Returns the local number of entries in this matrix.
  size_t getLocalNumEntries() const;

  //! Returns the current number of entries on this node in the specified local row.
  /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
  size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;

  //! Returns the current number of entries in the specified global row.
  /*! Returns OrdinalTraits<size_t>::invalid() if the row is not owned by this process. */
  size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  /** Undefined if isFillActive().
   */
  size_t getGlobalMaxNumRowEntries() const;

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  /** Undefined if isFillActive().
   */
  size_t getLocalMaxNumRowEntries() const;

  //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
  bool isLocallyIndexed() const;

  //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
  bool isGloballyIndexed() const;

  //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
  bool isFillComplete() const;

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
                       size_t &NumEntries) const;

  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*!
    \param GlobalRow - (In) Global row number for which indices are desired.
    \param Indices   - (Out) Global column indices corresponding to values.
    \param Values    - (Out) Row values
    \pre <tt>isLocallyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
  */
  void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const;

  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices  - (Out) Global column indices corresponding to values.
    \param Values   - (Out) Row values
    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

    Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
  */
  void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const;

  //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the
    the zero and non-zero diagonals owned by this node. */
  void getLocalDiagCopy(Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) const;

  //! Get offsets of the diagonal entries in the matrix.
  void getLocalDiagOffsets(Teuchos::ArrayRCP<size_t> &offsets) const;

  //! Get a copy of the diagonal entries owned by this node, with local row indices, using row offsets.
  void getLocalDiagCopy(Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag, const Teuchos::ArrayView<const size_t> &offsets) const;

  //! Get Frobenius norm of the matrix
  typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const;

  //! Left scale matrix using the given vector entries
  void leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &x);

  //! Right scale matrix using the given vector entries
  void rightScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &x);

  //! Returns true if globalConstants have been computed; false otherwise
  bool haveGlobalConstants() const;

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
  // TODO virtual=0 // TODO: Add default parameters ?
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
  virtual void apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
                     Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
                     Teuchos::ETransp mode = Teuchos::NO_TRANS,
                     Scalar alpha          = ScalarTraits<Scalar>::one(),
                     Scalar beta           = ScalarTraits<Scalar>::zero()) const;

  //! Computes the matrix-multivector multiplication for region layout matrices
  virtual void apply(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
                     MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
                     Teuchos::ETransp mode,
                     Scalar alpha,
                     Scalar beta,
                     bool sumInterfaceValues,
                     const RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > &regionInterfaceImporter,
                     const Teuchos::ArrayRCP<LocalOrdinal> &regionInterfaceLIDs) const;

  //! \brief Returns the Map associated with the domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const;

  //! Returns the Map associated with the domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const;

  //! \brief Returns the Map that describes the column distribution in this matrix.
  //! This might be <tt>null</tt> until fillComplete() is called.
  const RCP<const Map> &getColMap() const;

  //! \brief Returns the Map that describes the column distribution in this matrix.
  const RCP<const Map> &getColMap(viewLabel_t viewLabel) const;

  void removeEmptyProcessesInPlace(const Teuchos::RCP<const Map> &newMap);

  //@}

  //! Implements DistObject interface
  //{@

  //! Access function for the Tpetra::Map this DistObject was constructed with.
  const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getMap() const;

  //! Import.
  void doImport(const Matrix &source,
                const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM);

  //! Export.
  void doExport(const Matrix &dest,
                const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM);

  //! Import (using an Exporter).
  void doImport(const Matrix &source,
                const Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM);

  //! Export (using an Importer).
  void doExport(const Matrix &dest,
                const Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM);

  // @}

  //! @name Overridden from Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;

  //! @name Overridden from Teuchos::LabeledObject
  //@{
  void setObjectLabel(const std::string &objectLabel);
  //@}

#ifdef HAVE_XPETRA_TPETRA
  virtual local_matrix_type getLocalMatrixDevice() const;
  virtual typename local_matrix_type::HostMirror getLocalMatrixHost() const;
#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif

  // JG: Added:

  bool hasCrsGraph() const;

  //! Returns the CrsGraph associated with this matrix.
  RCP<const CrsGraph> getCrsGraph() const;

  RCP<CrsMatrix> getCrsMatrix() const;

  //! Returns the block size of the storage mechanism, which is usually 1, except for Tpetra::BlockCrsMatrix
  LocalOrdinal GetStorageBlockSize() const;

  //! Compute a residual R = B - (*this) * X
  void residual(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
                const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B,
                MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &R) const;

  //! Expert only
  void replaceCrsMatrix(RCP<CrsMatrix> &M);

  //@}
 private:
  // Default view is created after fillComplete()
  // Because ColMap might not be available before fillComplete().
  void CreateDefaultView();

 private:
  // The colMap can be <tt>null</tt> until fillComplete() is called. The default view of the Matrix have to be updated when fillComplete() is called.
  // If CrsMatrix::fillComplete() have been used instead of CrsMatrixWrap::fillComplete(), the default view is updated when getColMap() is called.
  void updateDefaultView() const;
  // The boolean finalDefaultView_ keep track of the status of the default view (= already updated or not)
  // See also CrsMatrixWrap::updateDefaultView()
  mutable bool finalDefaultView_;

  // The underlying matrix object
  RCP<CrsMatrix> matrixData_;

};  // class CrsMatrixWrap

}  // namespace Xpetra

#define XPETRA_CRSMATRIXWRAP_SHORT
#endif  // XPETRA_CRSMATRIXWRAP_DECL_HPP

// NOTE: if CrsMatrix and VbrMatrix share a common interface for fillComplete() etc, I can move some stuff in Xpetra_Matrix.hpp
// TODO: getUnderlyingMatrix() method
