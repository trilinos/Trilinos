/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_SPARSITYFILTER_DECL_HPP
#define IFPACK2_SPARSITYFILTER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Details_RowMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Ifpack2 {

/// \class SparsityFilter
/// \brief Drop entries of a matrix, based on the sparsity pattern.
/// \tparam MatrixType A specialization of Tpetra::RowMatrix.
///
/// This class takes an existing Tpetra::RowMatrix, and wraps it in a
/// matrix interface that drops entries using the following criteria:
///
/// 1. Drop elements beyond a certain distance from the diagonal, and / or
/// 2. Keep only the N largest entries in each row.
///
/// Here is a typical use case:
/// \code
/// Teuchos::RCP<Tpetra::RowMatrix<> > A = ...;
/// // First filter out all entries in columns
/// // which are not in A's domain Map.
/// Ifpack2::LocalFilter<Tpetra::RowMatrix<> > A_local (A);
/// // Drop all but the largest 20 elements in each row.
/// const size_t maxEntriesPerRow = 20;
/// // Exclude elements in each row whose local column index
/// // is more than 10 greater or less than the local column
/// // index of the diagonal element in that row.
/// const int maxBw = 10;
/// // Now create the sparsity filter. Elements dropped are
/// // not included in calls to getLocalRowCopy() and apply().
/// Ifpack2::SparsityFilter<Tpetra::RowMatrix<> > A_drop (A_local, maxEntriesPerRow, maxBw);
/// \endcode
///
/// This class currently only works if the matrix is the result of a
/// LocalFilter, that is, if the row and column Maps are the same.
template<class MatrixType>
class SparsityFilter :
    virtual public Ifpack2::Details::RowMatrix<MatrixType> {
public:
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> row_matrix_type;
  typedef typename row_matrix_type::mag_type mag_type;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::SparsityFilter: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.");


  //! \name Constructor & destructor methods
  //@{

  //! Constructor.
  explicit SparsityFilter (const Teuchos::RCP<const row_matrix_type>& Matrix,
                           size_t AllowedNumEntries,
                           LocalOrdinal AllowedBandwidth = -Teuchos::ScalarTraits<LocalOrdinal>::one());
  //! Destructor.
  virtual ~SparsityFilter();

  //@}
  //! \name Matrix Query Methods
  //@{

  //! Returns the communicator.
  virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! Returns the underlying node.
  virtual Teuchos::RCP<Node> getNode() const;

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
  virtual size_t getNodeNumRows() const;

  //! Returns the number of columns needed to apply the forward operator on this node, i.e., the number of elements listed in the column map.
  virtual size_t getNodeNumCols() const;

  //! Returns the index base for global indices for this matrix.
  virtual GlobalOrdinal getIndexBase() const;

  //! Returns the global number of entries in this matrix.
  virtual global_size_t getGlobalNumEntries() const;

  //! Returns the local number of entries in this matrix.
  virtual size_t getNodeNumEntries() const;

  //! \brief Returns the current number of entries on this node in the specified global row.
  /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
  virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

  //! Returns the current number of entries on this node in the specified local row.
  /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
  virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  virtual size_t getGlobalMaxNumRowEntries() const;

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  virtual size_t getNodeMaxNumRowEntries() const;

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
  virtual void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                                const Teuchos::ArrayView<GlobalOrdinal> &Indices,
                                const Teuchos::ArrayView<Scalar> &Values,
                                size_t &NumEntries) const;

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
  virtual void getLocalRowCopy(LocalOrdinal DropRow,
                               const Teuchos::ArrayView<LocalOrdinal> &Indices,
                               const Teuchos::ArrayView<Scalar> &Values,
                               size_t &NumEntries) const ;

  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*!
    \param GlobalRow - (In) Global row number for which indices are desired.
    \param Indices   - (Out) Global column indices corresponding to values.
    \param Values    - (Out) Row values
    \pre <tt>isDroplyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void getGlobalRowView(GlobalOrdinal GlobalRow,
                                Teuchos::ArrayView<const GlobalOrdinal> &indices,
                                Teuchos::ArrayView<const Scalar> &values) const;

  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*!
    \param DropRow - (In) Drop row number for which indices are desired.
    \param Indices  - (Out) Global column indices corresponding to values.
    \param Values   - (Out) Row values
    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInDropRow(DropRow)</tt>

    Note: If \c DropRow does not belong to this node, then \c indices is set to null.
  */
  virtual void getLocalRowView(LocalOrdinal DropRow,
                               Teuchos::ArrayView<const LocalOrdinal> &indices,
                               Teuchos::ArrayView<const Scalar> &values) const;

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

  //! Maximum allowed entries per row.
  size_t AllowedNumEntries_;
  //! Maximum allowed bandwidth.
  LocalOrdinal AllowedBandwidth_;
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
  mutable Teuchos::Array<LocalOrdinal> Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time
  mutable Teuchos::Array<Scalar> Values_;

};// class SparsityFilter

}// namespace Ifpack2

#endif /* IFPACK2_SPARSITYFILTER_DECL_HPP */
