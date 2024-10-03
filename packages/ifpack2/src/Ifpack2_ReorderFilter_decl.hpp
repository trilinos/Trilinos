// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_REORDERFILTER_DECL_HPP
#define IFPACK2_REORDERFILTER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Details_RowMatrix.hpp"
#include <type_traits>

namespace Ifpack2 {

/*!
\class ReorderFilter
\brief Wraps a Tpetra::RowMatrix in a filter that reorders local rows and columns.
\tparam MatrixType A specialization of Tpetra::RowMatrix.

This class is used in AdditiveSchwarz to reorder (if required by the
user) the localized matrix.  As the localized matrix is defined on a
serial communicator only, all maps are trivial (as all elements reside
on the same process).  This class does not attemp to define properly
reordered maps, hence it should not be used for distributed matrices.

To improve the performance of Ifpack2::AdditiveSchwarz, some
operations are not performed in the construction phase (like for
instance the computation of the 1-norm and infinite-norm, of check
whether the reordered matrix is lower/upper triangular or not).
*/
template<class MatrixType>
class ReorderFilter :
    virtual public Ifpack2::Details::RowMatrix<MatrixType> {
public:
  typedef typename MatrixType::scalar_type scalar_type;
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

  static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::ReorderFilter: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.  The constructor can take either a RowMatrix or a CrsMatrix just fine.");

  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;

  typedef typename row_matrix_type::mag_type mag_type;

  //! \name Constructor & destructor methods
  //@{

  /// \brief Constructor.
  ///
  /// \param A [in] The matrix to which to apply the filter.
  /// \param perm [in] Forward permutation of A's rows and columns.
  /// \param reverseperm [in] Reverse permutation of A's rows and columns.
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
  ReorderFilter (const Teuchos::RCP<const row_matrix_type>& A,
                 const Teuchos::ArrayRCP<local_ordinal_type>& perm,
                 const Teuchos::ArrayRCP<local_ordinal_type>& reverseperm);

  //! Destructor.
  virtual ~ReorderFilter ();

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

  //! The number of degrees of freedom per mesh point.
  virtual local_ordinal_type getBlockSize () const;

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

  //! Permute multivector: original-to-reordered
  virtual void permuteOriginalToReordered(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &originalX,
                                          Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &reorderedY) const;

  template <class DomainScalar, class RangeScalar>
  void
  permuteOriginalToReorderedTempl (const Tpetra::MultiVector<DomainScalar,local_ordinal_type,global_ordinal_type,node_type> &originalX,
                                       Tpetra::MultiVector<RangeScalar,local_ordinal_type,global_ordinal_type,node_type> &reorderedY) const;

  //! Permute multivector: reordered-to-original
  virtual void permuteReorderedToOriginal(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &reorderedX,
                                          Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &originalY) const;

  template <class DomainScalar, class RangeScalar>
  void permuteReorderedToOriginalTempl(const Tpetra::MultiVector<DomainScalar,local_ordinal_type,global_ordinal_type,node_type> &reorderedX,
                                       Tpetra::MultiVector<RangeScalar,local_ordinal_type,global_ordinal_type,node_type> &originalY) const;
  //@}

private:
  //! Pointer to the matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> A_;
  //! Permutation: Original to reordered
  Teuchos::ArrayRCP<local_ordinal_type> perm_;
  //! Permutation: Reordered to original
  Teuchos::ArrayRCP<local_ordinal_type> reverseperm_;

  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable nonconst_local_inds_host_view_type Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable nonconst_values_host_view_type Values_;
};// class ReorderFilter

}// namespace Ifpack2

#endif /* IFPACK2_REORDERFILTER_DECL_HPP */
