// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDCRSMATRIX_DECL_HPP
#define TPETRA_BLOCKEDCRSMATRIX_DECL_HPP

/// \file Tpetra_BlockedCrsMatrix_decl.hpp
/// \brief Declaration of the Tpetra::BlockedCrsMatrix class.
///
/// Direct port of Xpetra::BlockedCrsMatrix.  A BlockedCrsMatrix is an N x M
/// "matrix of matrices": an operator whose blocks are themselves
/// Tpetra::RowMatrix objects (typically Tpetra::CrsMatrix, possibly nested
/// Tpetra::BlockedCrsMatrix).  It derives from Tpetra::RowMatrix so that it can
/// be handed to any Tpetra consumer expecting a row matrix / operator.  The
/// map bookkeeping lives in two Tpetra::MapExtractor objects (range and
/// domain), and the block sub-vectors are moved with Tpetra::BlockedMultiVector.
///
/// This is the Tpetra-native core.  The Thyra constructor and
/// getThyraOperator() stay in the Xpetra wrapper (Thyra is not a Tpetra
/// dependency); "Thyra-style" numbering (bThyraMode) itself is plain
/// Tpetra::Map/Import arithmetic and is fully supported here.

#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"

#include "Tpetra_BlockedCrsMatrix_fwd.hpp"
#include "Tpetra_MapExtractor_fwd.hpp"
#include "Tpetra_BlockedMap_fwd.hpp"
#include "Tpetra_BlockedMultiVector_fwd.hpp"
#include "Tpetra_BlockedVector_fwd.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_RowGraph_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_Export_fwd.hpp"

#include "Tpetra_RowMatrix_decl.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_Vector_decl.hpp"

namespace Tpetra {

/// \class BlockedCrsMatrix
/// \brief Block operator ("matrix of matrices") built on Tpetra::RowMatrix blocks.
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class BlockedCrsMatrix : public ::Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using scalar_type         = Scalar;
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;

 private:
  using base_type          = ::Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Map                = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using MultiVector        = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Vector             = ::Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using RowMatrix          = ::Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using CrsMatrix          = ::Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using CrsGraph           = ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using RowGraph           = ::Tpetra::RowGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using Operator           = ::Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Import             = ::Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
  using Export             = ::Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>;
  using MapExtractor       = ::Tpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedMap         = ::Tpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedMultiVector = ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedVector      = ::Tpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

 public:
  //! The Kokkos-view types inherited from Tpetra::RowMatrix (used by overrides).
  using typename base_type::global_inds_host_view_type;
  using typename base_type::local_inds_host_view_type;
  using typename base_type::values_host_view_type;
  using typename base_type::nonconst_global_inds_host_view_type;
  using typename base_type::nonconst_local_inds_host_view_type;
  using typename base_type::nonconst_values_host_view_type;
  using typename base_type::mag_type;
  using typename base_type::impl_scalar_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor from range/domain BlockedMap objects.
  /*!
   * \param rangeMaps range maps for all blocks
   * \param domainMaps domain maps for all blocks
   * \param numEntriesPerRow estimated number of entries per row in each block
   */
  BlockedCrsMatrix(const Teuchos::RCP<const BlockedMap>& rangeMaps,
                   const Teuchos::RCP<const BlockedMap>& domainMaps,
                   size_t numEntriesPerRow);

  //! Constructor from range/domain MapExtractor objects.
  BlockedCrsMatrix(Teuchos::RCP<const MapExtractor>& rangeMapExtractor,
                   Teuchos::RCP<const MapExtractor>& domainMapExtractor,
                   size_t numEntriesPerRow);

  //! Destructor.
  virtual ~BlockedCrsMatrix();

  //@}

  //! @name Insertion/Removal Methods
  //@{

  //! Insert matrix entries, using global IDs (throws unless 1x1).
  void insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal>& cols, const Teuchos::ArrayView<const Scalar>& vals);

  //! Insert matrix entries, using local IDs (throws unless 1x1).
  void insertLocalValues(LocalOrdinal localRow, const Teuchos::ArrayView<const LocalOrdinal>& cols, const Teuchos::ArrayView<const Scalar>& vals);

  void removeEmptyProcessesInPlace(const Teuchos::RCP<const Map>& newMap);

  //! Replace matrix entries, using global IDs (throws unless 1x1).
  void replaceGlobalValues(GlobalOrdinal globalRow,
                           const Teuchos::ArrayView<const GlobalOrdinal>& cols,
                           const Teuchos::ArrayView<const Scalar>& vals);

  //! Replace matrix entries, using local IDs (throws unless 1x1).
  void replaceLocalValues(LocalOrdinal localRow,
                          const Teuchos::ArrayView<const LocalOrdinal>& cols,
                          const Teuchos::ArrayView<const Scalar>& vals);

  //! Set all matrix entries equal to scalar.
  virtual void setAllToScalar(const Scalar& alpha);

  //! Scale the current values of a matrix, this = alpha*this.
  void scale(const Scalar& alpha);

  //@}

  //! @name Transformational Methods
  //@{

  //! Resume fill operations (iteratively for all sub-blocks).
  void resumeFill(const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  //! Signal that data entry is complete (calls fillComplete for all underlying blocks).
  void fillComplete(const Teuchos::RCP<const Map>& domainMap, const Teuchos::RCP<const Map>& rangeMap, const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  //! Signal that data entry is complete.
  void fillComplete(const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  //@}

  //! @name Methods implementing RowMatrix
  //@{

  //! The communicator over which this matrix is distributed.
  Teuchos::RCP<const Teuchos::Comm<int>> getComm() const override;

  //! The (full range) Map that describes the distribution of rows over processes.
  Teuchos::RCP<const Map> getRowMap() const override;

  //! The (full domain) Map that describes the distribution of columns over processes.
  Teuchos::RCP<const Map> getColMap() const override;

  //! The RowGraph associated with this matrix (only supported for 1x1).
  Teuchos::RCP<const RowGraph> getGraph() const override;

  //! Returns the number of global rows.
  global_size_t getGlobalNumRows() const override;

  //! Returns the number of global columns in the matrix.
  global_size_t getGlobalNumCols() const override;

  //! Returns the number of matrix rows owned on the calling node.
  size_t getLocalNumRows() const override;

  //! Returns the number of columns connected to the locally owned rows.
  size_t getLocalNumCols() const override;

  //! The index base for global indices in this matrix.
  GlobalOrdinal getIndexBase() const override;

  //! Returns the global number of entries in this matrix.
  global_size_t getGlobalNumEntries() const override;

  //! Returns the local number of entries in this matrix.
  size_t getLocalNumEntries() const override;

  //! Returns the current number of entries in the specified (locally owned) global row.
  size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const override;

  //! Returns the current number of entries on this node in the specified local row.
  size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const override;

  //! Returns the maximum number of entries across all rows on all nodes.
  size_t getGlobalMaxNumRowEntries() const override;

  //! The number of degrees of freedom per mesh point.
  LocalOrdinal getBlockSize() const override;

  //! Returns the maximum number of entries across all rows on this node.
  size_t getLocalMaxNumRowEntries() const override;

  //! Whether this matrix has a well-defined column Map.
  bool hasColMap() const override;

  //! Whether matrix indices (of all blocks) are locally indexed.
  bool isLocallyIndexed() const override;

  //! Whether matrix indices (of all blocks) are globally indexed.
  bool isGloballyIndexed() const override;

  //! Whether fillComplete() has been called on all blocks.
  bool isFillComplete() const override;

  //! Whether this object implements getLocalRowView() and getGlobalRowView().
  bool supportsRowViews() const override;

  //! Extract a copy of a global row (only supported for 1x1).
  void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                        nonconst_global_inds_host_view_type& Indices,
                        nonconst_values_host_view_type& Values,
                        size_t& NumEntries) const override;

  //! Extract a copy of a local row (only supported for 1x1).
  void getLocalRowCopy(LocalOrdinal LocalRow,
                       nonconst_local_inds_host_view_type& Indices,
                       nonconst_values_host_view_type& Values,
                       size_t& NumEntries) const override;

  //! Extract a const view of a global row (only supported for 1x1).
  void getGlobalRowView(GlobalOrdinal GlobalRow,
                        global_inds_host_view_type& indices,
                        values_host_view_type& values) const override;

  //! Extract a const view of a local row (1x1 delegate, or diagonal path).
  void getLocalRowView(LocalOrdinal LocalRow,
                       local_inds_host_view_type& indices,
                       values_host_view_type& values) const override;

  //! Get a copy of the diagonal entries owned by this node.
  void getLocalDiagCopy(Vector& diag) const override;

  //! Left scale matrix using the given vector entries.
  void leftScale(const Vector& x) override;

  //! Right scale matrix using the given vector entries.
  void rightScale(const Vector& x) override;

  //! Get Frobenius norm of the matrix.
  mag_type getFrobeniusNorm() const override;

  //@}

  //! @name Methods implementing Operator
  //@{

  //! sparse matrix-multivector multiplication for the region layout matrices (no blocked implementation).
  virtual void apply(const MultiVector& X, MultiVector& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta, bool sumInterfaceValues,
                     const Teuchos::RCP<Import>& regionInterfaceImporter,
                     const Teuchos::ArrayRCP<LocalOrdinal>& regionInterfaceLIDs) const;

  //! Computes the sparse matrix-multivector multiplication Y = alpha*A^mode*X + beta*Y.
  void apply(const MultiVector& X, MultiVector& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const override;

  //! Returns the Map associated with the (full) domain of this operator.
  Teuchos::RCP<const Map> getDomainMap() const override;

  //! Returns the Map associated with the (full) range of this operator.
  Teuchos::RCP<const Map> getRangeMap() const override;

  //@}

  //! @name Block matrix map access
  //@{

  //! Returns the Map associated with the full domain of this operator.
  Teuchos::RCP<const Map> getFullDomainMap() const;

  //! Returns the BlockedMap associated with the domain of this operator.
  Teuchos::RCP<const BlockedMap> getBlockedDomainMap() const;

  //! Returns the Map associated with the i'th block domain of this operator.
  Teuchos::RCP<const Map> getDomainMap(size_t i) const;

  //! Returns the Map associated with the i'th block domain of this operator.
  Teuchos::RCP<const Map> getDomainMap(size_t i, bool bThyraMode) const;

  //! Returns the Map associated with the full range of this operator.
  Teuchos::RCP<const Map> getFullRangeMap() const;

  //! Returns the BlockedMap associated with the range of this operator.
  Teuchos::RCP<const BlockedMap> getBlockedRangeMap() const;

  //! Returns the Map associated with the i'th block range of this operator.
  Teuchos::RCP<const Map> getRangeMap(size_t i) const;

  //! Returns the Map associated with the i'th block range of this operator.
  Teuchos::RCP<const Map> getRangeMap(size_t i, bool bThyraMode) const;

  //! Returns map extractor class for range map.
  Teuchos::RCP<const MapExtractor> getRangeMapExtractor() const;

  //! Returns map extractor for domain map.
  Teuchos::RCP<const MapExtractor> getDomainMapExtractor() const;

  //@}

  //! @name Special multiplication routine (for BGS/Jacobi smoother)
  //@{

  //! Sparse matrix-multivector multiplication computed only for blocked row "row".
  virtual void bgs_apply(const MultiVector& X,
                         MultiVector& Y,
                         size_t row,
                         Teuchos::ETransp mode = Teuchos::NO_TRANS,
                         Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
                         Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //@}

  //! @name "DistObject"-like interface (delegate-or-throw, matching Xpetra)
  //@{

  //! Import (delegates to block (0,0) for 1x1, else throws).
  void doImport(const RowMatrix& source, const Import& importer, CombineMode CM);

  //! Export (delegates to block (0,0) for 1x1, else throws).
  void doExport(const RowMatrix& dest, const Import& importer, CombineMode CM);

  //! Import using an Exporter (delegates to block (0,0) for 1x1, else throws).
  void doImport(const RowMatrix& source, const Export& exporter, CombineMode CM);

  //! Export using an Importer (delegates to block (0,0) for 1x1, else throws).
  void doExport(const RowMatrix& dest, const Export& exporter, CombineMode CM);

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  std::string description() const override;

  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const override;

  //@}

  //! @name Overridden from Teuchos::LabeledObject
  //@{
  void setObjectLabel(const std::string& objectLabel) override;
  //@}

  //! Supports the getCrsGraph() call (only for 1x1).
  bool hasCrsGraph() const;

  //! Returns the CrsGraph associated with this matrix (only for 1x1).
  Teuchos::RCP<const CrsGraph> getCrsGraph() const;

  //! @name Block matrix access
  //@{

  virtual bool isDiagonal() const;

  //! Number of row blocks.
  virtual size_t Rows() const;

  //! Number of column blocks.
  virtual size_t Cols() const;

  //! Return unwrapped 1x1 blocked operators.
  Teuchos::RCP<RowMatrix> getCrsMatrix() const;

  //! Recursively returns the first inner-most non-null matrix block.
  Teuchos::RCP<RowMatrix> getInnermostCrsMatrix();

  //! Return block (r,c).
  Teuchos::RCP<RowMatrix> getMatrix(size_t r, size_t c) const;

  //! Set matrix block (r,c).
  void setMatrix(size_t r, size_t c, Teuchos::RCP<RowMatrix> mat);

  //! Merge BlockedCrsMatrix blocks into a single CrsMatrix.
  Teuchos::RCP<CrsMatrix> Merge() const;
  //@}

  //! Returns the block size of the storage mechanism.
  LocalOrdinal GetStorageBlockSize() const;

  //! Compute a residual R = B - (*this) * X.
  void residual(const MultiVector& X,
                const MultiVector& B,
                MultiVector& R) const;

 private:
  /** \name helper functions */
  //@{

  //! Add a Tpetra::RowMatrix to a Tpetra::CrsMatrix: B = B*scalarB + A*scalarA.
  /*!
   * Works only correctly if A only has entries which are empty (zero) in B.
   * Used only by Merge; since the blocks in a BlockedCrsMatrix are disjoint,
   * this suffices for merging.
   */
  void Add(const RowMatrix& A, const Scalar scalarA, CrsMatrix& B, const Scalar scalarB) const;

  //! Translate the Thyra-style column GIDs of a block into Xpetra-style GIDs.
  /*!
   * Ported from Xpetra::MapUtils::transformThyra2XpetraGIDs; pure Tpetra::Map
   * arithmetic plus Teuchos::reduceAll communication.  Only used by the
   * Thyra-mode branch of Merge().
   */
  static Teuchos::RCP<Map> transformThyra2XpetraGIDs(const Map& input,
                                                     const Map& nonOvlInput,
                                                     const Map& nonOvlReferenceInput);

  //@}

 private:
  bool is_diagonal_;                             ///< If we're diagonal, a bunch of the extraction stuff should work
  Teuchos::RCP<const MapExtractor> domainmaps_;  ///< full domain map together with all partial domain maps
  Teuchos::RCP<const MapExtractor> rangemaps_;   ///< full range map together with all partial range maps

  std::vector<Teuchos::RCP<RowMatrix>> blocks_;  ///< row major matrix block storage

  bool bRangeThyraMode_;   ///< true if built using Thyra-style numbering for range sub blocks
  bool bDomainThyraMode_;  ///< true if built using Thyra-style numbering for domain sub blocks
};

}  // namespace Tpetra

#endif  // TPETRA_BLOCKEDCRSMATRIX_DECL_HPP
