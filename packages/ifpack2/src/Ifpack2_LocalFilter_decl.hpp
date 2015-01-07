/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_LOCALFILTER_DECL_HPP
#define IFPACK2_LOCALFILTER_DECL_HPP

#include <Ifpack2_ConfigDefs.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <vector>


namespace Ifpack2 {

/// \class LocalFilter
/// \brief Access only local rows and columns of a sparse matrix.
/// \tparam MatrixType A specialization of either Tpetra::RowMatrix
///   (preferred) or Tpetra::CrsMatrix (which is a subclass of
///   Tpetra::RowMatrix).
///
/// \section Ifpack2_LocalFilter_Summary Summary
///
/// LocalFilter makes a global matrix "appear" local.  Solving a
/// linear system with the LocalFilter of a matrix A, if possible, is
/// equivalent to applying one step of nonoverlapping additive Schwarz
/// to A.  This class is an implementation detail of Ifpack2.  It is
/// not meant for users.  It may be useful to Ifpack2 developers who
/// want to implement a "local" preconditioner (i.e., that only works
/// within a single MPI process), and would like the preconditioner to
/// work for a global matrix.  For example, ILUT uses LocalFilter to
/// extract and compute the incomplete factorizaton of the "local
/// matrix."  This makes the input to ILUT appear like a square
/// matrix, assuming that the domain and range Maps have the same
/// number of entries on each process.  This in turn makes the result
/// of ILUT suitable for solving linear systems.
///
/// \section Ifpack2_LocalFilter_How How does it work?
///
/// \subsection Ifpack2_LocalFilter_How_Local View of the local rows and columns
///
/// LocalFilter provides a view of only the local rows and columns of
/// an existing sparse matrix.  "Local rows" are those rows in the row
/// Map that are owned by the range Map on the calling process.
/// "Local columns" are those columns in the column Map that are owned
/// by the domain Map on the calling process.  The view's communicator
/// contains only the local process (in MPI terms,
/// <tt>MPI_COMM_SELF</tt>), so each process will have its own
/// distinct view of its local part of the matrix.
///
/// \subsection Ifpack2_LocalFilter_How_DiagBlock Square diagonal block of the original matrix
///
/// If the following conditions hold on the Maps of a matrix A, then
/// the view resulting from a LocalFilter of A will be that of a
/// square diagonal block of A:
///   - Domain and range Maps are the same
///   - On every process, the row Map's entries are the same as the
///     range Map's entries, possibly followed by additional remote
///     entries
///   - On every process, the column Map's entries are the same as the
///     domain Map's entries, possibly followed by additional remote
///     entries
///
/// These conditions commonly hold for a Tpetra::CrsMatrix constructed
/// without a column Map, after fillComplete() has been called on it,
/// with default domain and range Maps.
///
/// \subsection Ifpack2_LocalFilter_How_Ind Remapping of global indices
///
/// The global indices of the view's domain and range Maps will be
/// different than those in the original matrix's domain and range
/// Maps.  The global indices of the new (domain, range) Map will be
/// remapped to a consecutive sequence, corresponding exactly to the
/// local indices of the original (domain, range) Map.  This ensures
/// that the local indices of the old and new Maps match.
///
/// \subsection Ifpack2_LocalFilter_How_Copy Not necessarily a copy
///
/// This class does not necessarily copy the entire sparse matrix.  It
/// may choose instead just to "filter out" the nonlocal entries.  The
/// effect on the LocalFilter of changing the entries or structure of
/// the original matrix is undefined.  LocalFilter may even make
/// different implementation choices on different MPI processes.  It
/// may do so because it does not invoke MPI communication outside of
/// the calling process.
///
/// \section Ifpack2_LocalFilter_Examples Examples
///
/// Here is an example of how to apply a LocalFilter to an existing
/// Tpetra sparse matrix:
/// \code
/// #include "Ifpack2_LocalFilter.hpp"
/// // ...
/// using Teuchos::RCP;
/// typedef Tpetra::RowMatrix<double> row_matrix_type;
/// typedef Tpetra::CrsMatrix<double> crs_matrix_type;
///
/// RCP<crs_matrix_type> A = ...;
/// // ... fill the entries of A ...
/// A->FillComplete ();
///
/// Ifpack2::LocalFilter<row_matrix_type> A_local (A);
/// \endcode
///
/// Here is an example of how to apply a LocalFilter, using \c A and
/// \c A_local from the example above.
/// \code
/// typedef Tpetra::Vector<double> vec_type;
///
/// // Make vectors x and y suitable for A->apply()
/// vec_type x (A.domainMap ());
/// vec_type y (A.rangeMap ());
/// // ... fill x ...
/// A->apply (x, y);
///
/// // Reinterpret x and y as local views suitable for A_local->apply()
/// RCP<const vec_type> x_local = x.offsetView (A_local.getDomainMap (), 0);
/// RCP<vec_type> y_local = y.offsetViewNonConst (A_local.getRangeMap (), 0);
///
/// // Apply the locally filtered version of A
/// A_local.apply (*x_local, *y_local);
/// \endcode
template<class MatrixType>
class LocalFilter :
    virtual public Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                     typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type>,
    virtual public Teuchos::Describable
{
public:
  //! \name Typedefs
  //@{

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! Preserved only for backwards compatibility.  Please use "scalar_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::scalar_type Scalar;


  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "local_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::local_ordinal_type LocalOrdinal;


  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "global_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::global_ordinal_type GlobalOrdinal;


  //! The Node type used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! Preserved only for backwards compatibility.  Please use "node_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::node_type Node;


  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Preserved only for backwards compatibility.  Please use "magnitude_type".
  TEUCHOS_DEPRECATED typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitudeType;

  //! Type of the Tpetra::RowMatrix specialization that this class uses.
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;

  //! Type of the Tpetra::Map specialization that this class uses.
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;

  typedef typename row_matrix_type::mag_type mag_type;

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  //! A one-line description of this object.
  virtual std::string description () const;

  //! Print the object to the given output stream.
  virtual void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;

  //@}
  //! \name Constructor and destructor
  //@{

  /// \brief Constructor
  ///
  /// \param A [in] The sparse matrix to which to apply the local filter.
  ///
  /// This class will <i>not</i> modify the input matrix.
  explicit LocalFilter (const Teuchos::RCP<const row_matrix_type>& A);

  //! Destructor
  virtual ~LocalFilter();

  //@}
  //! \name Matrix Query Methods
  //@{

  //! Returns the communicator.
  virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! Returns the underlying Node object.
  virtual Teuchos::RCP<node_type> getNode() const;

  //! Returns the Map that describes the row distribution in this matrix.
  virtual Teuchos::RCP<const map_type> getRowMap() const;

  //! Returns the Map that describes the column distribution in this matrix.
  virtual Teuchos::RCP<const map_type> getColMap() const;

  //! Returns the Map that describes the domain distribution in this matrix.
  virtual Teuchos::RCP<const map_type> getDomainMap() const;

  //! Returns the Map that describes the range distribution in this matrix.
  virtual Teuchos::RCP<const map_type> getRangeMap() const;

  //! The (locally filtered) matrix's graph.
  virtual Teuchos::RCP<const Tpetra::RowGraph<local_ordinal_type,global_ordinal_type,node_type> >
  getGraph () const;

  //! The number of global rows in this matrix.
  virtual global_size_t getGlobalNumRows() const;

  //! The number of global columns in this matrix.
  virtual global_size_t getGlobalNumCols() const;

  //! The number of rows owned on the calling process.
  virtual size_t getNodeNumRows() const;

  //! The number of columns in the (locally filtered) matrix.
  virtual size_t getNodeNumCols() const;

  //! Returns the index base for global indices for this matrix.
  virtual global_ordinal_type getIndexBase() const;

  //! Returns the global number of entries in this matrix.
  virtual global_size_t getGlobalNumEntries() const;

  //! Returns the local number of entries in this matrix.
  virtual size_t getNodeNumEntries() const;

  /// \brief The current number of entries on this node in the specified global row.
  ///
  /// \return <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt> if
  ///   the specified row is not owned by this process, otherwise the
  ///   number of entries in that row on this process.
  virtual size_t getNumEntriesInGlobalRow (global_ordinal_type globalRow) const;

  /// \brief The current number of entries on this node in the specified local row.
  ///
  /// \return <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt> if
  ///   the specified local row is not valid on this process,
  ///   otherwise the number of entries in that row on this process.
  virtual size_t getNumEntriesInLocalRow (local_ordinal_type localRow) const;

  //! The number of global diagonal entries, based on global row/column index comparisons.
  virtual global_size_t getGlobalNumDiags() const;

  //! The number of local diagonal entries, based on global row/column index comparisons.
  virtual size_t getNodeNumDiags() const;

  //! The maximum number of entries across all rows/columns on all processes.
  virtual size_t getGlobalMaxNumRowEntries() const;

  //! The maximum number of entries across all rows/columns on this process.
  virtual size_t getNodeMaxNumRowEntries() const;

  //! Whether this matrix has a well-defined column Map.
  virtual bool hasColMap() const;

  //! Whether this matrix is lower triangular.
  virtual bool isLowerTriangular() const;

  //! Whether this matrix is upper triangular.
  virtual bool isUpperTriangular() const;

  //! Whether the underlying sparse matrix is locally (opposite of globally) indexed.
  virtual bool isLocallyIndexed() const;

  //! Whether the underlying sparse matrix is globally (opposite of locally) indexed.
  virtual bool isGloballyIndexed() const;

  //! Returns \c true if fillComplete() has been called.
  virtual bool isFillComplete() const;

  //! Returns \c true if RowViews are supported.
  virtual bool supportsRowViews() const;

  //@}
  //! \name Methods for getting the entries in a row.
  //@{

  /// \brief Get the entries in the given row, using global indices.
  ///
  /// \param GlobalRow [i] Global index of the row for which indices are desired.
  /// \param Indices [out] Global indices of the entries in the given row.
  /// \param Values [out] Values of the entries in the given row.
  /// \param NumEntries [out] Number of entries in the row.
  ///
  /// This method throws an exception if either output array is not
  /// large enough to hold the data associated with row \c
  /// GlobalRow. If \c GlobalRow does not belong to the calling
  /// process, then \c Indices and \c Values are unchanged and
  /// \c NumIndices is <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>
  /// on output.
  virtual void
  getGlobalRowCopy (global_ordinal_type GlobalRow,
                    const Teuchos::ArrayView<global_ordinal_type> &Indices,
                    const Teuchos::ArrayView<scalar_type> &Values,
                    size_t &NumEntries) const;

  /// \brief Get the entries in the given row, using local indices.
  ///
  /// \param LocalRow [i] Local index of the row for which indices are desired.
  /// \param Indices [out] Local indices of the entries in the given row.
  /// \param Values [out] Values of the entries in the given row.
  /// \param NumEntries [out] Number of entries in the row.
  ///
  /// This method throws an exception if either output array is not
  /// large enough to hold the data associated with row \c
  /// LocalRow. If \c LocalRow does not belong to the calling
  /// process, then \c Indices and \c Values are unchanged and
  /// \c NumIndices is <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>
  /// on output.
  virtual void
  getLocalRowCopy (local_ordinal_type LocalRow,
                   const Teuchos::ArrayView<local_ordinal_type> &Indices,
                   const Teuchos::ArrayView<scalar_type> &Values,
                   size_t &NumEntries) const ;

  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*!
    \param GlobalRow [in] Global row number for which indices are desired.
    \param Indices [out] Global column indices corresponding to values.
    \param Values [out] Row values

    \pre <tt>isLocallyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void
  getGlobalRowView (global_ordinal_type GlobalRow,
                    Teuchos::ArrayView<const global_ordinal_type> &indices,
                    Teuchos::ArrayView<const scalar_type> &values) const;

  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*!
    \param LocalRow [in] Local row number for which indices are desired.
    \param Indices [out] Local column indices corresponding to values.
    \param Values [out] Row values

    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

    Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void
  getLocalRowView (local_ordinal_type LocalRow,
                   Teuchos::ArrayView<const local_ordinal_type> &indices,
                   Teuchos::ArrayView<const scalar_type> &values) const;

  /// \brief Get the diagonal entries of the (locally filtered) matrix.
  ///
  /// \param diag [in/out] On input: a Tpetra::Vector whose Map is the
  ///   same as the local filter's row Map, which in turn is the same
  ///   as the original matrix's row Map.  On output: filled with the
  ///   diagonal entries owned by the calling process.
  virtual void
  getLocalDiagCopy (Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> &diag) const;

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

  /// \brief The Frobenius norm of the (locally filtered) matrix.
  ///
  /// This method may return a different value on each process,
  /// because this method computes the norm of the locally filtered
  /// matrix, which may be different on each process.
  ///
  /// The Frobenius norm of a matrix \f$A\f$ is defined as
  /// \f$\|A\|_F = \sqrt{\sum_{i,j} \|A_{ij}\|^2}\f$.
  virtual mag_type getFrobeniusNorm() const;

  /// \brief Compute Y = beta*Y + alpha*A_local*X.
  ///
  /// If <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so
  /// that any values in \c Y (including NaNs) are ignored.  If
  /// <tt>alpha == 0</tt>, apply() must short-circuit the operator, so
  /// that any values in \c X (including NaNs) are ignored.
  virtual void
  apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> &X,
         Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> &Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Whether this operator supports applying the transpose or conjugate transpose.
  virtual bool hasTransposeApply() const;

  //@}
  //! \name Deprecated routines to be removed at some point in the future.
  //@{

  //! Deprecated. Get a persisting const view of the entries in a specified global row of this matrix.
  /*!
    \param GlobalRow - (In) Global row from which to retrieve matrix entries.
    \param Indices - (Out) Indices for the global row.
    \param Values - (Out) Values for the global row.

    Note: If \c GlobalRow does not belong to this node, then \c Indices and \c Values are set to <tt>Teuchos::null</t>>.

    \pre isLocallyIndexed()==false
  */
  TPETRA_DEPRECATED
  virtual void getGlobalRowView(global_ordinal_type GlobalRow,
                                Teuchos::ArrayRCP<const global_ordinal_type> &indices,
                                Teuchos::ArrayRCP<const scalar_type>        &values) const;

  //! Deprecated. Get a persisting const view of the entries in a specified local row of this matrix.
  /*!
    \param LocalRow - (In) Local row from which to retrieve matrix entries.
    \param Indices - (Out) Indices for the local row.
    \param Values - (Out) Values for the local row.

    Note: If \c LocalRow is not valid for this node, then \c Indices and \c Values are set to <tt>Teuchos::null</tt>.

    \pre isGloballyIndexed()==false
  */
  TPETRA_DEPRECATED
  virtual void getLocalRowView(local_ordinal_type LocalRow,
                               Teuchos::ArrayRCP<const local_ordinal_type> &indices,
                               Teuchos::ArrayRCP<const scalar_type>       &values) const;
  //@}


private:
  //! Type of a read-only interface to a graph.
  typedef Tpetra::RowGraph<local_ordinal_type,
                           global_ordinal_type,
                           node_type> row_graph_type;
  //! Special case of apply() for when X and Y do not alias one another.
  void
  applyNonAliased (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
                   Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
                   Teuchos::ETransp mode,
                   scalar_type alpha,
                   scalar_type beta) const;

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
  Teuchos::RCP<const row_matrix_type> A_;

  //! Row Map of the locally filtered matrix.
  Teuchos::RCP<const map_type> localRowMap_;

  //! Domain Map of the locally filtered matrix.
  Teuchos::RCP<const map_type> localDomainMap_;

  //! Range Map of the locally filtered matrix.
  Teuchos::RCP<const map_type> localRangeMap_;

  //! Number of nonzeros in the local matrix.
  size_t NumNonzeros_;
  //! Maximum number of nonzero entries in a row for the filtered matrix.
  size_t MaxNumEntries_;
  //! Maximum number of nonzero entries in a row for the unfiltered matrix.
  size_t MaxNumEntriesA_;
  //! NumEntries_[i] contains the nonzero entries in row `i'.
  std::vector<size_t> NumEntries_;

  //! Temporary array used in getLocalRowCopy().
  mutable Teuchos::Array<local_ordinal_type> localIndices_;

  //! Temporary array used in getLocalRowCopy().
  mutable Teuchos::Array<scalar_type> Values_;
};// class LocalFilter

}// namespace Ifpack2

#endif /* IFPACK2_LOCALFILTER_DECL_HPP */
