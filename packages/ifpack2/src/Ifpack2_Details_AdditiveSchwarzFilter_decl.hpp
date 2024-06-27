
// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_ADDITIVE_SCHWARZ_FILTER_DECL_HPP
#define IFPACK2_ADDITIVE_SCHWARZ_FILTER_DECL_HPP

/*!
\class AdditiveSchwarzFilter
\brief Wraps a Tpetra::CrsMatrix or Ifpack2::OverlappingRowMatrix in a filter that removes off-process edges, may reorder rows/columns, and may remove singletons (rows with no connections to other rows).

This is essentially a fusion of LocalFilter, AdditiveSchwarzFilter and SingletonFilter. These three filters are used together in AdditiveSchwarz.
*/

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Details_RowMatrix.hpp"
#include "Ifpack2_OverlappingRowMatrix.hpp"

namespace Ifpack2
{
namespace Details
{
  template<typename MatrixType>
  class AdditiveSchwarzFilter :
    public Ifpack2::Details::RowMatrix<MatrixType> {
public:
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename Kokkos::ArithTraits<scalar_type>::val_type impl_scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef typename MatrixType::global_inds_host_view_type global_inds_host_view_type;
  typedef typename MatrixType::local_inds_host_view_type local_inds_host_view_type;
  typedef typename MatrixType::values_host_view_type values_host_view_type;

  typedef typename MatrixType::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef typename MatrixType::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename MatrixType::nonconst_values_host_view_type nonconst_values_host_view_type;

  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;
  typedef Tpetra::MultiVector<scalar_type,
                              local_ordinal_type,
                              global_ordinal_type,
                              node_type> mv_type;
  typedef typename crs_matrix_type::device_type device_type;
  typedef typename crs_matrix_type::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space> policy_type;
  typedef Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<2>> policy_2d_type;
  typedef Ifpack2::OverlappingRowMatrix<row_matrix_type> overlapping_matrix_type;
  typedef typename crs_matrix_type::local_matrix_device_type local_matrix_type;
  typedef typename local_matrix_type::size_type size_type;
  typedef typename local_matrix_type::row_map_type::non_const_type row_map_type;
  typedef typename local_matrix_type::index_type entries_type;
  typedef typename local_matrix_type::values_type values_type;
  typedef typename row_map_type::HostMirror host_row_map_type;
  typedef typename entries_type::HostMirror host_entries_type;
  typedef typename values_type::HostMirror host_values_type;
  typedef typename local_matrix_type::HostMirror host_local_matrix_type;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::AdditiveSchwarzFilter: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.  The constructor can take either a RowMatrix or a CrsMatrix just fine.");

  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;

  typedef typename row_matrix_type::mag_type mag_type;

  //! \name Constructor & destructor methods
  //@{

  /// \brief Constructor.
  ///
  /// \param A [in] The matrix to which to apply the filter.
  //                Must be a Tpetra::CrsMatrix or an Ifpack2::OverlappingRowMatrix that wraps a Tpetra::CrsMatrix.
  /// \param perm [in] Forward permutation of A's rows and columns.
  /// \param reverseperm [in] Reverse permutation of A's rows and columns.
  /// \param filterSingletons [in] If true, remove rows that have no local neighbors.
  ///
  /// It must make sense to apply the given permutation to both the
  /// rows and columns.  This means that the row and column Maps must
  /// have the same numbers of entries on all processes, and must have
  /// the same order of GIDs on all processes.
  ///
  /// perm[i] gives the where OLD index i shows up in the NEW
  /// ordering.  revperm[i] gives the where NEW index i shows up in
  /// the OLD ordering.  Note that perm is actually the "inverse
  /// permutation," in Zoltan2 terms.
  AdditiveSchwarzFilter(const Teuchos::RCP<const row_matrix_type>& A,
                 const Teuchos::ArrayRCP<local_ordinal_type>& perm,
                 const Teuchos::ArrayRCP<local_ordinal_type>& reverseperm,
                 bool filterSingletons);

  //! Update the filtered matrix in response to A's values changing (A's structure isn't allowed to change, though)
  void updateMatrixValues();

  Teuchos::RCP<const crs_matrix_type> getFilteredMatrix() const;
  
  //! Destructor.
  virtual ~AdditiveSchwarzFilter ();

  //@}
  //! \name Matrix query methods
  //@{

  //! The matrix's communicator.
  virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! Returns the Map that describes the row distribution in this matrix.
  virtual Teuchos::RCP<const map_type> getRowMap() const;

  //! Returns the Map that describes the column distribution in this matrix.
  virtual Teuchos::RCP<const map_type> getColMap() const;

  //! Returns the Map that describes the domain distribution in this matrix.
  virtual Teuchos::RCP<const map_type> getDomainMap() const;

  //! Returns the Map that describes the range distribution in this matrix.
  virtual Teuchos::RCP<const map_type> getRangeMap() const;

  //! Returns the RowGraph associated with this matrix.
  virtual Teuchos::RCP<const Tpetra::RowGraph<local_ordinal_type,global_ordinal_type,node_type> > getGraph() const;

  //! The number of degrees of freedom per mesh point.
  virtual local_ordinal_type getBlockSize() const;

  //! Returns the number of global rows in this matrix.
  virtual global_size_t getGlobalNumRows() const;

  //! \brief Returns the number of global columns in this matrix.
  virtual global_size_t getGlobalNumCols() const;

  //! Returns the number of rows owned on the calling node.
  virtual size_t getLocalNumRows() const;

  //! Returns the number of columns needed to apply the forward operator on this node, i.e., the number of elements listed in the column map.
  virtual size_t getLocalNumCols() const;

  //! Returns the index base for global indices for this matrix.
  virtual global_ordinal_type getIndexBase() const;

  //! Returns the global number of entries in this matrix.
  virtual global_size_t getGlobalNumEntries() const;

  //! Returns the local number of entries in this matrix.
  virtual size_t getLocalNumEntries() const;

  /// \brief The current number of entries in this matrix, stored on
  ///   the calling process, in the row whose global index is \c globalRow.
  ///
  /// \return The number of entries, or
  ///   Teuchos::OrdinalTraits<size_t>::invalid() if the specified row
  ///   is not owned by the calling process.
  virtual size_t getNumEntriesInGlobalRow (global_ordinal_type globalRow) const;

  /// \brief The current number of entries in this matrix, stored on
  ///   the calling process, in the row whose local index is \c globalRow.
  ///
  /// \return The number of entries, or
  ///   Teuchos::OrdinalTraits<size_t>::invalid() if the specified row
  ///   is not owned by the calling process.
  virtual size_t getNumEntriesInLocalRow (local_ordinal_type localRow) const;

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  virtual size_t getGlobalMaxNumRowEntries() const;

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  virtual size_t getLocalMaxNumRowEntries() const;

  //! \brief Indicates whether this matrix has a well-defined column map.
  virtual bool hasColMap() const;

  //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
  virtual bool isLocallyIndexed() const;

  //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
  virtual bool isGloballyIndexed() const;

  //! Returns \c true if fillComplete() has been called.
  virtual bool isFillComplete() const;

  //! Returns \c true if RowViews are supported.
  virtual bool supportsRowViews() const;

  //@}

  //! @name Extraction Methods
  //@{

  //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
  /*!
    \param GlobalRow  - (In) Global row number for which indices are desired.
    \param Indices    - (Out) Global column indices corresponding to values.
    \param Values     - (Out) Matrix values.
    \param NumEntries - (Out) Number of indices.

    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices and \c Values are unchanged and \c NumIndices is
    returned as Teuchos::OrdinalTraits<size_t>::invalid().
  */
  virtual void
  getGlobalRowCopy (global_ordinal_type GlobalRow,
                   nonconst_global_inds_host_view_type &Indices,
                   nonconst_values_host_view_type &Values,
                   size_t& NumEntries) const;
  //! Extract a list of entries in a specified local row of the graph. Put into storage allocated by calling routine.
  /*!
    \param LocalRow   - (In) Local row number for which indices are desired.
    \param Indices    - (Out) Local column indices corresponding to values.
    \param Values     - (Out) Matrix values.
    \param NumIndices - (Out) Number of indices.

    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is
    returned as Teuchos::OrdinalTraits<size_t>::invalid().
  */
  virtual void
  getLocalRowCopy (local_ordinal_type LocalRow,
                   nonconst_local_inds_host_view_type &Indices,
                   nonconst_values_host_view_type &Values,
                   size_t& NumEntries) const;

  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*!
    \param GlobalRow - (In) Global row number for which indices are desired.
    \param Indices   - (Out) Global column indices corresponding to values.
    \param Values    - (Out) Row values
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>
    \pre <tt>isLocallyIndexed() == false</tt>
    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void
  getGlobalRowView (global_ordinal_type GlobalRow,
                    global_inds_host_view_type &indices,
                    values_host_view_type &values) const;
  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices  - (Out) Global column indices corresponding to values.
    \param Values   - (Out) Row values
    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInDropRow(LocalRow)</tt>

    Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void
  getLocalRowView (local_ordinal_type LocalRow,
                   local_inds_host_view_type & indices,
                   values_host_view_type & values) const;
  //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the
    the zero and non-zero diagonals owned by this node. */
  virtual void getLocalDiagCopy(Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &diag) const;

  //@}

  //! \name Mathematical Methods
  //@{

  /**
   * \brief Scales the RowMatrix on the left with the Vector x.
   *
   * This matrix will be scaled such that A(i,j) = x(i)*A(i,j)
   * where i denoes the global row number of A and
   * j denotes the global column number of A.
   *
   * \param x A vector to left scale this matrix.
   */
  virtual void leftScale(const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x);

  /**
   * \brief Scales the RowMatrix on the right with the Vector x.
   *
   * This matrix will be scaled such that A(i,j) = x(j)*A(i,j)
   * where i denoes the global row number of A and
   * j denotes the global column number of A.
   *
   * \param x A vector to right scale this matrix.
   */
  virtual void rightScale(const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x);

  //! Returns the Frobenius norm of the matrix.
  /** Computes and returns the Frobenius norm of the matrix, defined as:
      \f$ \|A\|_F = \sqrt{\sum_{i,j} \|\a_{ij}\|^2} \f$
  */
  virtual mag_type getFrobeniusNorm() const;

  /// \brief \f$ Y := \beta Y + \alpha Op(A) X \f$,
  ///   where Op(A) is either A, \f$A^T\f$, or \f$A^H\f$.
  ///
  /// Apply the reordered version of the matrix (or its transpose or
  /// conjugate transpose) to the given multivector X, producing Y.
  ///   - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y,
  ///     so that any values in \c Y (including NaNs) are ignored.
  ///   - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the
  ///     matrix, so that any values in \c X (including NaNs) are
  ///     ignored.
  ///
  /// This method assumes that X and Y are in the reordered order.
  ///
  /// If hasTransposeApply() returns false, then the only valid value
  /// of \c mode is Teuchos::NO_TRANS (the default).  Otherwise, it
  /// accepts the following values:
  ///   - mode = Teuchos::NO_TRANS: Op(A) is the reordered version of A.
  ///   - mode = Teuchos::TRANS: Op(A) is the reordered version of the
  ///     transpose of A.
  ///   - mode = Teuchos::CONJ_TRANS: Op(A) is reordered version of
  ///     the conjugate transpose of A.
  virtual void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Whether apply() can apply the transpose or conjugate transpose.
  virtual bool hasTransposeApply() const;

  //! Solve singleton rows of OverlappingY, then filter and permute OverlappingB to get the reduced B.
  void CreateReducedProblem(
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& OverlappingB,
      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>&       OverlappingY,
      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>&       ReducedReorderedB) const;

  //! Scatter ReducedY back to non-singleton rows of OverlappingY, according to the reordering.
  void UpdateLHS(
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& ReducedReorderedY,
      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>&       OverlappingY) const;

  void setup(const Teuchos::RCP<const row_matrix_type>& A_unfiltered,
      const Teuchos::ArrayRCP<local_ordinal_type>& perm,
      const Teuchos::ArrayRCP<local_ordinal_type>& reverseperm,
      bool filterSingletons);

  /// \brief Fill the entries/values in the local matrix, and from that construct A_.
  ///
  /// Requires that the structure of A_unfiltered_ has not changed since this was constructed.
  /// Used by both the constructor and updateMatrixValues().
  void fillLocalMatrix(local_matrix_type localMatrix);

  //@}

private:

  /// \brief Whether map1 is fitted to map2.
  ///
  /// \param map1 [in] The first map.
  /// \param map2 [in] The second map.
  ///
  /// For Tpetra::Map instances map1 and map2, we say that map1 is
  /// <i>fitted</i> to map2 (on the calling process), when the initial
  /// indices of map1 (on the calling process) are the same and in the
  /// same order as those of map2.  "Fittedness" is entirely a local
  /// (per MPI process) property.  The predicate "map1 is fitted to
  /// map2?" is <i>not</i> symmetric.  For example, map2 may have more
  /// entries than map1.
  static bool
  mapPairIsFitted (const map_type& map1, const map_type& map2);

  /// \brief Whether the domain Map of A is fitted to its column Map,
  ///   and the range Map of A is fitted to its row Map.
  ///
  // If both pairs of Maps of the original matrix A are fitted on this
  // process, then this process can use a fast "view" implementation.
  static bool
  mapPairsAreFitted (const row_matrix_type& A);

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> A_unfiltered_;
  //! Filtered and reordered matrix.
  Teuchos::RCP<crs_matrix_type> A_;
  //! Permutation: Original to reordered. Size is number of local rows in A_unfiltered_. Rows which were filtered from A_ map to -1 (Invalid).
  Kokkos::DualView<local_ordinal_type*, device_type> perm_;
  //! Permutation: Reordered to original. Size is number of local rows in A_.
  Kokkos::DualView<local_ordinal_type*, device_type> reverseperm_;
  local_ordinal_type numSingletons_;
  //! List of singleton rows (indices are rows in A_unfiltered_)
  Kokkos::DualView<local_ordinal_type*, device_type> singletons_;
  //! Diagonal values of singletons (used for solving those rows)
  Kokkos::DualView<impl_scalar_type*, device_type> singletonDiagonals_;

  //! Whether singletons (disconnected rows) are filtered out of the local matrix
  // (apply and solve still process those rows separately)
  bool FilterSingletons_;

  //! Row, col, domain and range map of this locally filtered matrix (it's square and non-distributed)
  Teuchos::RCP<const map_type> localMap_;
};

}}  //namespace Ifpack2::Details

#endif

