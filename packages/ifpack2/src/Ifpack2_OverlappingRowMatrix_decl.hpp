// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_OVERLAPPINGROWMATRIX_DECL_HPP
#define IFPACK2_OVERLAPPINGROWMATRIX_DECL_HPP

#include "Ifpack2_Details_RowMatrix.hpp"
#include "Tpetra_CrsMatrix_decl.hpp" // only need the declaration here
#include "Tpetra_Import_decl.hpp"
#include "Tpetra_Map_decl.hpp"
#include <type_traits>

namespace Ifpack2 {

/// \class OverlappingRowMatrix
/// \brief Sparse matrix (Tpetra::RowMatrix subclass) with ghost rows.
/// \tparam MatrixType Tpetra::RowMatrix specialization.
template<class MatrixType>
class OverlappingRowMatrix :
    virtual public Ifpack2::Details::RowMatrix<MatrixType> {
public:
  //! \name Typedefs
  //@{
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
  typedef typename MatrixType::global_inds_host_view_type global_inds_host_view_type;
  typedef typename MatrixType::local_inds_host_view_type local_inds_host_view_type;
  typedef typename MatrixType::values_host_view_type values_host_view_type;

  typedef typename MatrixType::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef typename MatrixType::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename MatrixType::nonconst_values_host_view_type nonconst_values_host_view_type;

  using row_matrix_type = Tpetra::RowMatrix<scalar_type, local_ordinal_type,
					    global_ordinal_type, node_type>;
  using crs_matrix_type = Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
					    global_ordinal_type, node_type>;

  // device typedefs
  typedef typename MatrixType::node_type::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef typename MatrixType::local_inds_device_view_type local_inds_device_view_type;
  typedef typename MatrixType::global_inds_device_view_type global_inds_device_view_type;
  typedef typename MatrixType::values_device_view_type values_device_view_type;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::OverlappingRowMatrix: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.  The constructor can take either a RowMatrix or a CrsMatrix just fine.");

  typedef typename row_matrix_type::mag_type mag_type;

  //@}
  //! \name Constructors and destructor
  //@{

  /// Standard constructor.
  ///
  /// \param A [in] The input matrix.  Currently this class requires
  ///   that A be a Tpetra::CrsMatrix instance with the same first
  ///   four template parameters as MatrixType, and with a default
  ///   fifth template parameter.  Furthermore, A must have a
  ///   nonoverlapping row Map and must be distributed over more than
  ///   one MPI process.
  ///
  /// \param overlapLevel [in] The number of levels of overlap.
  OverlappingRowMatrix (const Teuchos::RCP<const row_matrix_type>& A,
                        const int overlapLevel);

  //! Destructor
  ~OverlappingRowMatrix () = default;

  //@}
  //! @name Matrix query methods
  //@{

  //! The communicator over which the matrix is distributed.
  virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;


  //! The Map that describes the distribution of rows over processes.
  virtual Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> >
  getRowMap () const;

  //! The Map that describes the distribution of columns over processes.
  virtual Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> >
  getColMap () const;

  /// \brief The Map that describes the domain of this matrix.
  ///
  /// The domain is the distribution of valid input vectors of apply().
  virtual Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> >
  getDomainMap () const;

  /// \brief The Map that describes the range of this matrix.
  ///
  /// The domain is the distribution of valid output vectors of apply().
  virtual Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> >
  getRangeMap () const;

  //! This matrix's graph.
  virtual Teuchos::RCP<const Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type> >
  getGraph () const;

  //! The global number of rows in this matrix.
  virtual global_size_t getGlobalNumRows () const;

  //! The global number of columns in this matrix.
  virtual global_size_t getGlobalNumCols () const;

  //! The number of rows owned by the calling process.
  virtual size_t getLocalNumRows () const;

  /// \brief The number of columns owned by the calling process.
  ///
  /// This is the number of columns needed to apply the forward
  /// operator on the calling process, that is, the number of elements
  /// listed in the column Map on the calling process.
  virtual size_t getLocalNumCols () const;

  //! The index base for global indices for this matrix.
  virtual global_ordinal_type getIndexBase () const;

  //! The global number of entries in this matrix.
  virtual global_size_t getGlobalNumEntries () const;

  //! The number of entries in this matrix owned by the calling process.
  virtual size_t getLocalNumEntries () const;

  /// \brief The number of entries in the given global row that are
  ///   owned by the calling process.
  ///
  /// \param globalRow [in] Global index of the row.
  ///
  /// \return Teuchos::OrdinalTraits<size_t>::invalid() if the
  ///   specified row (either in the input matrix or in the overlap
  ///   matrix) is not owned by the calling process, else the number
  ///   of entries in that row that are owned by the calling process.
  virtual size_t getNumEntriesInGlobalRow (global_ordinal_type globalRow) const;

  /// \brief The number of entries in the given local row that are
  ///   owned by the calling process.
  ///
  /// \param globalRow [in] Local index of the row.
  ///
  /// \return Teuchos::OrdinalTraits<size_t>::invalid() if the
  ///   specified row (either in the input matrix or in the overlap
  ///   matrix) is not owned by the calling process, else the number
  ///   of entries in that row that are owned by the calling process.
  virtual size_t getNumEntriesInLocalRow (local_ordinal_type localRow) const;

  //! The maximum number of entries in any row on any process.
  virtual size_t getGlobalMaxNumRowEntries () const;

  //! The maximum number of entries in any row on the calling process.
  virtual size_t getLocalMaxNumRowEntries() const;
  
  //! The number of degrees of freedom per mesh point.
  virtual local_ordinal_type getBlockSize () const;

  //! Whether this matrix has a column Map.
  virtual bool hasColMap() const;

  //! Whether this matrix is locally indexed.
  virtual bool isLocallyIndexed () const;

  //! Whether this matrix is globally indexed.
  virtual bool isGloballyIndexed () const;

  //! \c true if fillComplete() has been called, else \c false.
  virtual bool isFillComplete() const;

  //! \c true if row views are supported, else \c false.
  virtual bool supportsRowViews() const;

  //@}
  //! @name Extraction methods
  //@{

  //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
  /*!
    \param LocalRow - (In) Global row number for which indices are desired.
    \param Indices - (Out) Global column indices corresponding to values.
    \param Values - (Out) Matrix values.
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
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices - (Out) Local column indices corresponding to values.
    \param Values - (Out) Matrix values.
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
    \pre <tt>isLocallyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

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
    \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

    Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
  */
 virtual void
  getLocalRowView (local_ordinal_type LocalRow,
                   local_inds_host_view_type & indices,
                   values_host_view_type & values) const;

  //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the
    the zero and non-zero diagonals owned by this node. */
  virtual
  void getLocalDiagCopy (Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &diag) const;

  //@}
  //! \name Mathematical methods
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
  virtual void
  leftScale (const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x);

  /**
   * \brief Scales the RowMatrix on the right with the Vector x.
   *
   * This matrix will be scaled such that A(i,j) = x(j)*A(i,j)
   * where i denoes the global row number of A and
   * j denotes the global column number of A.
   *
   * \param x A vector to right scale this matrix.
   */
  virtual void
  rightScale (const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x);

  //! Returns the Frobenius norm of the matrix.
  /** Computes and returns the Frobenius norm of the matrix, defined as:
      \f$ \|A\|_F = \sqrt{\sum_{i,j} \|\a_{ij}\|^2} \f$
  */
  virtual mag_type
  getFrobeniusNorm () const;

  //! \brief Computes the operator-multivector application.
  /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
    vary according to the values of \c alpha and \c beta. Specifically
    - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
    - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.

    This is analagous to the *Multiply* function in Ifpack, not the *Apply*
  */
  virtual void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Whether this operator's apply() method can apply the adjoint (transpose).
  virtual bool hasTransposeApply() const;

  virtual void
  importMultiVector (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
                     Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &OvX,
                     Tpetra::CombineMode CM = Tpetra::INSERT);

  virtual void
  exportMultiVector (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &OvX,
                     Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
                     Tpetra::CombineMode CM = Tpetra::ADD);

  std::string description() const;

  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const;

  Teuchos::RCP<const crs_matrix_type> getUnderlyingMatrix() const;

  Teuchos::RCP<const crs_matrix_type> getExtMatrix() const;

  Kokkos::View<size_t*, typename OverlappingRowMatrix<MatrixType>::device_type> getExtHaloStarts() const;
  typename Kokkos::View<size_t*, typename OverlappingRowMatrix<MatrixType>::device_type>::HostMirror getExtHaloStartsHost() const;

  void doExtImport();

private:
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
  typedef Tpetra::Export<local_ordinal_type, global_ordinal_type, node_type> export_type;
  typedef Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type> row_graph_type;
  typedef Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> vector_type;

  //! The input matrix to the constructor.
  Teuchos::RCP<const crs_matrix_type> A_;

  Tpetra::global_size_t NumGlobalRows_;
  Tpetra::global_size_t NumGlobalNonzeros_;
  size_t MaxNumEntries_;
  int OverlapLevel_;

  // Wrapper matrix objects
  Teuchos::RCP<const map_type> RowMap_;
  Teuchos::RCP<const map_type> ColMap_;
  Teuchos::RCP<const import_type> Importer_;

  //! The matrix containing only the overlap rows.
  Teuchos::RCP<crs_matrix_type> ExtMatrix_;
  Teuchos::RCP<const map_type>        ExtMap_;
  Teuchos::RCP<const import_type>     ExtImporter_;
  Kokkos::View<size_t*, device_type>  ExtHaloStarts_;
  typename Kokkos::View<size_t*, device_type>::HostMirror ExtHaloStarts_h;

  //! Graph of the matrix (as returned by getGraph()).
  Teuchos::RCP<const row_graph_type> graph_;
  //! Used in apply(), to avoid allocation each time.
  mutable nonconst_local_inds_host_view_type Indices_;
  //! Used in apply(), to avoid allocation each time.
  mutable nonconst_values_host_view_type Values_;


}; // class OverlappingRowMatrix

} // namespace Ifpack2

#endif // IFPACK2_OVERLAPPINGROWMATRIX_DECL_HPP
