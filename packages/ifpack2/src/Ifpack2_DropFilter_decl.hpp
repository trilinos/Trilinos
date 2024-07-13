// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DROPFILTER_DECL_HPP
#define IFPACK2_DROPFILTER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Ifpack2_Details_RowMatrix.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <type_traits>

namespace Ifpack2 {

/// \class DropFilter
/// \brief Filter based on matrix entries
/// \tparam MatrixType Tpetra::RowMatrix specialization
///
/// This class filters a matrix, by dropping all elements whose
/// absolute value is below a specified threshold.  The matrix should
/// be "localized," which means that it has no remote entries (e.g.,
/// is a LocalFilter).
///
/// A typical use is as follows:
/// \code
/// Teuchos::RCP<Tpetra::RowMatrix<> > A;
/// // first localize the matrix
/// Ifpack2::LocalFilter<Tpetra::RowMatrix<> > LocalA(A);
/// // drop all elements below this value
/// double DropTol = 1e-5;
/// // now create the matrix, elements below DropTol are
/// // not included in calls to getLocalRowCopy() and apply()
/// // and Apply()
/// Ifpack2::DropFilter<Tpetra::RowMatrix<> > DropA(LocalA,DropTol)
/// \endcode
template<class MatrixType>
class DropFilter :
    virtual public Ifpack2::Details::RowMatrix<MatrixType> {
public:
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  typedef typename MatrixType::global_inds_host_view_type global_inds_host_view_type;
  typedef typename MatrixType::local_inds_host_view_type local_inds_host_view_type;
  typedef typename MatrixType::values_host_view_type values_host_view_type;

  typedef typename MatrixType::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef typename MatrixType::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename MatrixType::nonconst_values_host_view_type nonconst_values_host_view_type;

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  typedef typename row_matrix_type::mag_type mag_type;


  static_assert(std::is_same<MatrixType, Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::value, "Ifpack2::DropFilter: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.  The constructor can take either a RowMatrix or a CrsMatrix just fine.");

  //! \name Constructor & destructor methods
  //@{

  //! Constructor.
  explicit DropFilter(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix,
                      magnitudeType DropTol);
  //! Destructor.
  virtual ~DropFilter();

  //@}

  //! \name Matrix Query Methods
  //@{

  //! Returns the communicator.
  virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;


  //! Returns the Map that describes the row distribution in this matrix.
  virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const;

  //! \brief Returns the Map that describes the column distribution in this matrix.
  virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const;

  //! Returns the Map that describes the domain distribution in this matrix.
  virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const;

  //! \brief Returns the Map that describes the range distribution in this matrix.
  virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const;

  //! Returns the RowGraph associated with this matrix.
  virtual Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const;

  //! Returns the number of global rows in this matrix.
  virtual global_size_t getGlobalNumRows() const;

  //! \brief Returns the number of global columns in this matrix.
  virtual global_size_t getGlobalNumCols() const;

  //! Returns the number of rows owned on the calling node.
  virtual size_t getLocalNumRows() const;

  //! Returns the number of columns needed to apply the forward operator on this node, i.e., the number of elements listed in the column map.
  virtual size_t getLocalNumCols() const;

  //! Returns the index base for global indices for this matrix.
  virtual GlobalOrdinal getIndexBase() const;

  //! Returns the global number of entries in this matrix.
  virtual global_size_t getGlobalNumEntries() const;

  //! Returns the local number of entries in this matrix.
  virtual size_t getLocalNumEntries() const;

  //! \brief Returns the current number of entries on this node in the specified global row.
  /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
  virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

  //! Returns the current number of entries on this node in the specified local row.
  /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
  virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  virtual size_t getGlobalMaxNumRowEntries() const;

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  virtual size_t getLocalMaxNumRowEntries() const;
  
  //! The number of degrees of freedom per mesh point.
  virtual LocalOrdinal getBlockSize () const;

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
    \param DropRow - (In) Global row number for which indices are desired.
    \param Indices - (Out) Global column indices corresponding to values.
    \param Values - (Out) Matrix values.
    \param NumEntries - (Out) Number of indices.

    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices and \c Values are unchanged and \c NumIndices is
    returned as Teuchos::OrdinalTraits<size_t>::invalid().
  */
  virtual void
  getGlobalRowCopy (GlobalOrdinal GlobalRow,
                    nonconst_global_inds_host_view_type &Indices,
                    nonconst_values_host_view_type &Values,
                    size_t& NumEntries) const;

  //! Extract a list of entries in a specified local row of the graph. Put into storage allocated by calling routine.
  /*!
    \param DropRow - (In) Drop row number for which indices are desired.
    \param Indices - (Out) Drop column indices corresponding to values.
    \param Values - (Out) Matrix values.
    \param NumIndices - (Out) Number of indices.

    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c DropRow. If \c DropRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is
    returned as Teuchos::OrdinalTraits<size_t>::invalid().
  */
  virtual void
  getLocalRowCopy (LocalOrdinal DropRow,
                   nonconst_local_inds_host_view_type &Indices,
                   nonconst_values_host_view_type &Values,
                   size_t& NumEntries) const;

  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*!
    \param GlobalRow - (In) Global row number for which indices are desired.
    \param Indices   - (Out) Global column indices corresponding to values.
    \param Values    - (Out) Row values
    \pre <tt>isDroplyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void
  getGlobalRowView (GlobalOrdinal GlobalRow,
                    global_inds_host_view_type &indices,
                    values_host_view_type &values) const;

  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*!
    \param DropRow - (In) Drop row number for which indices are desired.
    \param Indices  - (Out) Global column indices corresponding to values.
    \param Values   - (Out) Row values
    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInDropRow(DropRow)</tt>

    Note: If \c DropRow does not belong to this node, then \c indices is set to null.
  */
  virtual void
  getLocalRowView (LocalOrdinal LocalRow,
                   local_inds_host_view_type & indices,
                   values_host_view_type & values) const;

  //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the
    the zero and non-zero diagonals owned by this node. */
  virtual void getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const;

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
  virtual void leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x);

  /**
   * \brief Scales the RowMatrix on the right with the Vector x.
   *
   * This matrix will be scaled such that A(i,j) = x(j)*A(i,j)
   * where i denoes the global row number of A and
   * j denotes the global column number of A.
   *
   * \param x A vector to right scale this matrix.
   */
  virtual void rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x);

  //! Returns the Frobenius norm of the matrix.
  /** Computes and returns the Frobenius norm of the matrix, defined as:
      \f$ \|A\|_F = \sqrt{\sum_{i,j} \|\a_{ij}\|^2} \f$
  */
  virtual mag_type getFrobeniusNorm() const;

  //! \brief Computes the operator-multivector application.
  /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
    vary according to the values of \c alpha and \c beta. Specifically
    - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
    - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
  */
  virtual void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                     Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                     Teuchos::ETransp mode = Teuchos::NO_TRANS,
                     Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                     Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Indicates whether this operator supports applying the adjoint operator.
  virtual bool hasTransposeApply() const;

  //@}

private:

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  //! Drop Tolerance
  magnitudeType DropTol_;
  //! Number of rows in the local matrix.
  size_t NumRows_;
  //! Number of nonzeros in the local matrix.
  size_t NumNonzeros_;
  //! Maximum number of nonzero entries in a row for the filtered matrix.
  size_t MaxNumEntries_;
  //! Maximum number of nonzero entries in a row for the unfiltered matrix.
  size_t MaxNumEntriesA_;
  //! NumEntries_[i] contains the nonzero entries in row `i'.
  std::vector<size_t> NumEntries_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable nonconst_local_inds_host_view_type Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time
  mutable nonconst_values_host_view_type Values_;

};

} // namespace Ifpack2

#endif /* IFPACK2_DROPFILTER_DECL_HPP */
