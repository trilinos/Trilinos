// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_BLOCKEDCRSMATRIX_DECL_HPP
#define XPETRA_BLOCKEDCRSMATRIX_DECL_HPP

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_BlockedMultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_BlockedVector.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"

#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MapExtractorFactory.hpp"

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"

#ifdef HAVE_XPETRA_THYRA
#include <Thyra_ProductVectorSpaceBase.hpp>
#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_LinearOpBase.hpp>
#include <Thyra_BlockedLinearOpBase.hpp>
#include <Thyra_PhysicallyBlockedLinearOpBase.hpp>
#include "Xpetra_ThyraUtils.hpp"
#endif

#include "Xpetra_VectorFactory.hpp"

/** \file Xpetra_BlockedCrsMatrix.hpp

  Declarations for the class Xpetra::BlockedCrsMatrix.
*/
namespace Xpetra {

#ifdef HAVE_XPETRA_THYRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class ThyraUtils;
#endif

typedef std::string viewLabel_t;

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class BlockedCrsMatrix : public Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

 private:
#undef XPETRA_BLOCKEDCRSMATRIX_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  /*!
   * \param rangeMaps range maps for all blocks
   * \param domainMaps domain maps for all blocks
   * \param numEntriesPerRow estimated number of entries per row in each block(!)
   */
  BlockedCrsMatrix(const Teuchos::RCP<const BlockedMap>& rangeMaps,
                   const Teuchos::RCP<const BlockedMap>& domainMaps,
                   size_t numEntriesPerRow);

  //! Constructor
  /*!
   * \param rangeMapExtractor range map extractor for all blocks
   * \param domainMapExtractor domain map extractor for all blocks
   * \param numEntriesPerRow estimated number of entries per row in each block(!)
   *
   * \note This constructor will be deprecated. Please use the constructor which takes BlockedMap objects instead.
   */
  BlockedCrsMatrix(Teuchos::RCP<const MapExtractor>& rangeMapExtractor,
                   Teuchos::RCP<const MapExtractor>& domainMapExtractor,
                   size_t numEntriesPerRow);

#ifdef HAVE_XPETRA_THYRA
  //! Constructor
  /*!
   * \param rangeMaps range maps for all blocks
   * \param domainMaps domain maps for all blocks
   * \param npr extimated number of entries per row in each block(!)
   */
  BlockedCrsMatrix(const Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> >& thyraOp, const Teuchos::RCP<const Teuchos::Comm<int> >& /* comm */);

 private:
  //! mergeMaps
  /*!
   * \param subMaps
   *
   * Merges all Xpetra::Map objects in std::vector subMaps and returns a new Xpetra::Map containing all entries.
   * Helper function only used in constructor of Xpetra_BlockedCrsMatrix for transforming a Thyra::BlockedLinearOp object
   * All GID entries are sorted and duplicates are eliminated.
   */
  Teuchos::RCP<const Map> mergeMaps(std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& subMaps);

 public:
#endif

  //! Destructor
  virtual ~BlockedCrsMatrix();

  //@}

  //! @name Insertion/Removal Methods
  //@{

  //! Insert matrix entries, using global IDs.
  /**
    Note: this routine throws for Rows() > 1 and/or Cols() > 1

    All index values must be in the global space.
    \pre \c globalRow exists as an ID in the global row map
    \pre <tt>isLocallyIndexed() == false</tt>
    \pre <tt>isStorageOptimized() == false</tt>

    \post <tt>isGloballyIndexed() == true</tt>

    \note If \c globalRow does not belong to the matrix on this node, then it
    will be communicated to the appropriate node when globalAssemble() is
    called (which will, at the latest, occur during the next call to
    fillComplete().) Otherwise, the entries will be inserted in the local
    matrix.

    \note If the matrix row already contains values at the indices
    corresponding to values in \c cols, then the new values will be summed
    with the old values; this may happen at insertion or during the next call
    to fillComplete().

    \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where
    cols[i] belongs to the column map on this node will be inserted into the
    matrix.
    */
  void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal>& cols, const ArrayView<const Scalar>& vals);

  //! Insert matrix entries, using local IDs.
  /**
     Note: this routine throws if Rows() > 1 and/or Cols() > 1

     All index values must be in the local space.
    \pre \c localRow exists as an ID in the local row map
    \pre <tt>isGloballyIndexed() == false</tt>
    \pre <tt>isStorageOptimized() == false</tt>

    \post <tt>isLocallyIndexed() == true</tt>
    */
  void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal>& cols, const ArrayView<const Scalar>& vals);

  void removeEmptyProcessesInPlace(const Teuchos::RCP<const Map>& newMap);

  //! \brief Replace matrix entries, using global IDs.
  /** All index values must be in the global space.

    \pre \c globalRow is a global row belonging to the matrix on this node.

    \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in
    this matrix row (likely because it was inserted more than once and
    fillComplete() has not been called in the interim), the behavior of this
    function is not defined. */
  void replaceGlobalValues(GlobalOrdinal globalRow,
                           const ArrayView<const GlobalOrdinal>& cols,
                           const ArrayView<const Scalar>& vals);

  //! Replace matrix entries, using local IDs.
  /** All index values must be in the local space.
    Note that if a value is not already present for the specified location in
    the matrix, the input value will be ignored silently.
    */
  void replaceLocalValues(LocalOrdinal localRow,
                          const ArrayView<const LocalOrdinal>& cols,
                          const ArrayView<const Scalar>& vals);

  //! Set all matrix entries equal to scalar
  virtual void setAllToScalar(const Scalar& alpha);

  //! Scale the current values of a matrix, this = alpha*this.
  void scale(const Scalar& alpha);

  //@}

  //! @name Transformational Methods
  //@{

  /*! Resume fill operations.
    After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

    For BlockedCrsMatrix objects we call the routine iteratively for all sub-blocks.

    resumeFill() may be called repeatedly.
    */
  void resumeFill(const RCP<ParameterList>& params = null);

  /*! \brief Signal that data entry is complete.

    Note: for blocked operators the specified domain and range maps have no meaning.
          We just call fillComplete for all underlying blocks
    */
  void fillComplete(const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, const RCP<ParameterList>& params = null);

  /*! \brief Signal that data entry is complete.

    Off-node entries are distributed (via globalAssemble()), repeated entries are
    summed, and global indices are transformed to local indices.

    \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */
  void fillComplete(const RCP<ParameterList>& params = null);

  //@}

  //! Returns the number of global rows.
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

  //! Returns the current number of entries in the specified (locally owned) global row.
  /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
  size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  /** Undefined if isFillActive().
   */
  size_t getGlobalMaxNumRowEntries() const;

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  /** Undefined if isFillActive().
   */
  size_t getLocalMaxNumRowEntries() const;

  //! \brief If matrix indices of all matrix blocks are in the local range, this function returns true. Otherwise, this function returns false.
  /** if false, then this does not automatically mean that all blocks are globally indexed. The user has to make sure, that all matrix blocks
   * are indexed in the same way (locally or globally). Otherwise the block matrix is not valid...
   */
  bool isLocallyIndexed() const;

  //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false.
  /** if false, then this does not automatically mean that all blocks are locally indexed. The user has to make sure, that all matrix blocks
   * are indexed in the same way (locally or globally). Otherwise the block matrix is not valid...
   */
  bool isGloballyIndexed() const;

  //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
  bool isFillComplete() const;

  //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices - (Out) Local column indices corresponding to values.
    \param Values - (Out) Matrix values.
    \param NumIndices - (Out) Number of indices.

    Note: A std::runtime_error exception is thrown if either \c Indices or \c
    Values is not large enough to hold the data associated with row \c
    LocalRow. If \c LocalRow is not valid for this node, then \c Indices and
    \c Values are unchanged and \c NumIndices is returned as
    OrdinalTraits<size_t>::invalid().

    \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
    */
  virtual void getLocalRowCopy(LocalOrdinal LocalRow,
                               const ArrayView<LocalOrdinal>& Indices,
                               const ArrayView<Scalar>& Values,
                               size_t& NumEntries) const;

  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*!
    \param GlobalRow - (In) Global row number for which indices are desired.
    \param Indices   - (Out) Global column indices corresponding to values.
    \param Values    - (Out) Row values
    \pre <tt>isLocallyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
    */
  void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal>& indices, ArrayView<const Scalar>& values) const;

  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices  - (Out) Global column indices corresponding to values.
    \param Values   - (Out) Row values
    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

    Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
    */
  void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal>& indices, ArrayView<const Scalar>& values) const;

  //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
  /*! Returns a distributed Vector object partitioned according to this
    matrix's row map, containing the
    the zero and non-zero diagonals owned by this node. */
  void getLocalDiagCopy(Vector& diag) const;

  //! Left scale matrix using the given vector entries
  void leftScale(const Vector& x);

  //! Right scale matrix using the given vector entries
  void rightScale(const Vector& x);

  //! Get Frobenius norm of the matrix
  virtual typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const;

  //! Returns true if globalConstants have been computed; false otherwise
  virtual bool haveGlobalConstants() const;

  //@}

  //! @name Advanced Matrix-vector multiplication and solve methods
  //@{

  //! Multiplies this matrix by a MultiVector.
  /*! \c X is required to be post-imported, i.e., described by the column map
   * of the matrix. \c Y is required to be pre-exported, i.e., described by the
   * row map of the matrix.

   Both are required to have constant stride, and they are not permitted to
   ocupy overlapping space. No runtime checking will be performed in a
   non-debug build.

   This method is templated on the scalar type of MultiVector objects, allowing
   this method to be applied to MultiVector objects of arbitrary type. However,
   it is recommended that multiply() not be called directly; instead, use the
   CrsMatrixMultiplyOp, as it will handle the import/exprt operations required
   to apply a matrix with non-trivial communication needs.

   If \c beta is equal to zero, the operation will enjoy overwrite semantics
   (\c Y will be overwritten with the result of the multiplication). Otherwise,
   the result of the multiplication will be accumulated into \c Y.
   */

  //@}

  //! @name Methods implementing Matrix
  //@{

  //! sparse matrix-multivector multiplication for the region layout matrices (currently no blocked implementation)
  virtual void apply(const MultiVector& X, MultiVector& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta, bool sumInterfaceValues,
                     const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter,
                     const Teuchos::ArrayRCP<LocalOrdinal>& regionInterfaceLIDs) const;

  //! \brief Computes the sparse matrix-multivector multiplication.
  /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
    - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
    */
  virtual void apply(const MultiVector& X, MultiVector& Y,
                     Teuchos::ETransp mode = Teuchos::NO_TRANS,
                     Scalar alpha          = ScalarTraits<Scalar>::one(),
                     Scalar beta           = ScalarTraits<Scalar>::zero()) const;

  //! \brief Returns the Map associated with the full domain of this operator.
  RCP<const Map> getFullDomainMap() const;

  //! \brief Returns the BlockedMap associated with the domain of this operator.
  RCP<const BlockedMap> getBlockedDomainMap() const;

  //! \brief Returns the Map associated with the domain of this operator.
  const RCP<const Map> getDomainMap() const;

  //! \brief Returns the Map associated with the i'th block domain of this operator.
  RCP<const Map> getDomainMap(size_t i) const;

  //! \brief Returns the Map associated with the i'th block domain of this operator.
  RCP<const Map> getDomainMap(size_t i, bool bThyraMode) const;

  //! Returns the Map associated with the full range of this operator.
  RCP<const Map> getFullRangeMap() const;

  //! \brief Returns the BlockedMap associated with the range of this operator.
  RCP<const BlockedMap> getBlockedRangeMap() const;

  //! Returns the Map associated with the range of this operator.
  const RCP<const Map> getRangeMap() const;

  //! Returns the Map associated with the i'th block range of this operator.
  RCP<const Map> getRangeMap(size_t i) const;

  //! Returns the Map associated with the i'th block range of this operator.
  RCP<const Map> getRangeMap(size_t i, bool bThyraMode) const;

  //! Returns map extractor class for range map
  RCP<const MapExtractor> getRangeMapExtractor() const;

  //! Returns map extractor for domain map
  RCP<const MapExtractor> getDomainMapExtractor() const;

  //@}

  //! Special multiplication routine (for BGS/Jacobi smoother)
  //{@

  /*! \brief Computes the sparse matrix-multivector multiplication (plus linear combination with input/result vector)
   *
   *  Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exception:
   *  - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
   *  - calculates result only for blocked row "row"
   *  - useful for BGS/Jacobi smoother in MueLu: there we have to calculate the residual for the current block row
   *    we can skip the MatVec calls in all other block rows
   */
  virtual void bgs_apply(
      const MultiVector& X,                                 ///< Vector to be multiplied by matrix (input)
      MultiVector& Y,                                       ///< result vector
      size_t row,                                           ///< Index of block row to be treated
      Teuchos::ETransp mode = Teuchos::NO_TRANS,            ///< Transpose mode
      Scalar alpha          = ScalarTraits<Scalar>::one(),  ///< scaling factor for result of matrix-vector product
      Scalar beta           = ScalarTraits<Scalar>::zero()  ///< scaling factor for linear combination with result vector
  ) const;

  //@}

  //! Implements DistObject interface
  //{@

  //! Access function for the Tpetra::Map this DistObject was constructed with.
  const Teuchos::RCP<const Map> getMap() const;

  //! Import.
  void doImport(const Matrix& source, const Import& importer, CombineMode CM);

  //! Export.
  void doExport(const Matrix& dest, const Import& importer, CombineMode CM);

  //! Import (using an Exporter).
  void doImport(const Matrix& source, const Export& exporter, CombineMode CM);

  //! Export (using an Importer).
  void doExport(const Matrix& dest, const Export& exporter, CombineMode CM);

  // @}

  //! @name Overridden from Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;

  //! @name Overridden from Teuchos::LabeledObject
  //@{
  void setObjectLabel(const std::string& objectLabel);
  //@}

  //! Supports the getCrsGraph() call
  bool hasCrsGraph() const;

  //! Returns the CrsGraph associated with this matrix.
  RCP<const CrsGraph> getCrsGraph() const;

  //@}

  //! @name Block matrix access
  //@{

  virtual bool isDiagonal() const;

  /// number of row blocks
  virtual size_t Rows() const;

  /// number of column blocks
  virtual size_t Cols() const;

  /// return unwrap 1x1 blocked operators
  Teuchos::RCP<Matrix> getCrsMatrix() const;

  /// helper routine recursively returns the first inner-most non-null matrix block from a (nested) blocked operator
  Teuchos::RCP<Matrix> getInnermostCrsMatrix();

  /// return block (r,c)
  Teuchos::RCP<Matrix> getMatrix(size_t r, size_t c) const;

  /// set matrix block
  // void setMatrix(size_t r, size_t c, Teuchos::RCP<CrsMatrix> mat) {
  void setMatrix(size_t r, size_t c, Teuchos::RCP<Matrix> mat);

  /// merge BlockedCrsMatrix blocks in a CrsMatrix
  // NOTE: This is a rather expensive operation, since all blocks are copied
  // into a new big CrsMatrix
  Teuchos::RCP<Matrix> Merge() const;
  //@}

  typedef typename CrsMatrix::local_matrix_type local_matrix_type;
  /// \brief Access the underlying local Kokkos::CrsMatrix object
  local_matrix_type getLocalMatrixDevice() const;
  /// \brief Access the underlying local Kokkos::CrsMatrix object
  typename local_matrix_type::HostMirror getLocalMatrixHost() const;

#ifdef HAVE_XPETRA_THYRA
  Teuchos::RCP<Thyra::BlockedLinearOpBase<Scalar> > getThyraOperator();
#endif
  //! Returns the block size of the storage mechanism
  LocalOrdinal GetStorageBlockSize() const;

  //! Compute a residual R = B - (*this) * X
  void residual(const MultiVector& X,
                const MultiVector& B,
                MultiVector& R) const;

 private:
  /** \name helper functions */
  //@{

  /// Add a Xpetra::CrsMatrix to another: B = B*scalarB + A*scalarA
  /**
   * Note, that this routine works only correctly if A only has entries which are empty (zero) in B.
   * We use the insertGlobalValues routine for inserting the new values from A in B. The sumIntoGlobalValues
   * routine is not implemented in Xpetra (and would not extend the graph of B for new entries).
   * Here we need something to catch the exceptions of a future implementation of sumIntoGlobalValues that
   * then adds the remaining new entries with insertGlobal Values.
   *
   * This routine is private and used only by Merge. Since the blocks in BlockedCrsMatrix are seperated,
   * this routine works for merging a BlockedCrsMatrix.
   */
  void Add(const Matrix& A, const Scalar scalarA, Matrix& B, const Scalar scalarB) const;

  //@}

  // Default view is created after fillComplete()
  // Because ColMap might not be available before fillComplete().
  void CreateDefaultView();

 private:
  bool is_diagonal_;                             ///< If we're diagonal, a bunch of the extraction stuff should work
  Teuchos::RCP<const MapExtractor> domainmaps_;  ///< full domain map together with all partial domain maps
  Teuchos::RCP<const MapExtractor> rangemaps_;   ///< full range map together with all partial domain maps

  std::vector<Teuchos::RCP<Matrix> > blocks_;  ///< row major matrix block storage
#ifdef HAVE_XPETRA_THYRA
  Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> > thyraOp_;  ///< underlying thyra operator
#endif
  bool bRangeThyraMode_;   ///< boolean flag, which is true, if BlockedCrsMatrix has been created using Thyra-style numbering for sub blocks, i.e. all GIDs of submaps are contiguous and start from 0.
  bool bDomainThyraMode_;  ///< boolean flag, which is true, if BlockedCrsMatrix has been created using Thyra-style numbering for sub blocks, i.e. all GIDs of submaps are contiguous and start from 0.
};

}  // namespace Xpetra

#define XPETRA_BLOCKEDCRSMATRIX_SHORT
#endif /* XPETRA_BLOCKEDCRSMATRIX_DECL_HPP */
