// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_ROWMATRIX_HPP
#define TPETRA_ROWMATRIX_HPP

#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_RowGraph.hpp"

namespace Tpetra {
  //
  // Forward declarations.  The "doxygen" bit simply tells Doxygen
  // (our automatic documentation generation system) to skip forward
  // declarations.
  //
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  class Map;

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class Vector;
#endif // DOXYGEN_SHOULD_SKIP_THIS

  //! \brief A pure virtual interface for row-partitioned matrices.
  /*!
     This class is templated on \c Scalar, \c LocalOrdinal, \c GlobalOrdinal and \c Node.
     The \c LocalOrdinal type, if omitted, defaults to \c int.
     The \c GlobalOrdinal type defaults to the \c LocalOrdinal type.
     The \c Node type defaults to the default node in Kokkos.
   */
  template <class Scalar,
            class LocalOrdinal = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node = Kokkos::DefaultNode::DefaultNodeType>
  class RowMatrix :
    virtual public Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
  public:
    //! @name Typedefs
    //@{

    //! Same as the \c Scalar typedef of Operator.
    typedef Scalar        scalar_type;
    //! Same as the \c LocalOrdinal typedef of Operator.
    typedef LocalOrdinal  local_ordinal_type;
    //! Same as the \c GlobalOrdinal typedef of Operator.
    typedef GlobalOrdinal global_ordinal_type;
    //! Same as the \c Node typedef of Operator.
    typedef Node          node_type;

    //@}
    //! @name Destructor
    //@{

    //! Destructor.
    virtual ~RowMatrix();

    //@}
    //! @name Matrix query methods
    //@{

    //! The communicator over which this matrix is distributed.
    virtual const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const = 0;

    //! The Kokkos Node instance.
    virtual Teuchos::RCP<Node> getNode() const = 0;

    //! The Map that describes the distribution of rows over processes.
    virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const = 0;

    //! The Map that describes the distribution of columns over processes.
    virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const = 0;

    //! The RowGraph associated with this matrix.
    virtual Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const = 0;

    //! The global number of rows of this matrix.
    virtual global_size_t getGlobalNumRows() const = 0;

    //! The global number of columns of this matrix.
    virtual global_size_t getGlobalNumCols() const = 0;

    //! The number of rows owned by the calling process.
    virtual size_t getNodeNumRows() const = 0;

    /// \brief The number of columns needed to apply the forward operator on this node.
    ///
    /// This is the same as the number of elements listed in the
    /// column Map.  It is <i>not</i> necessarily the same as the
    /// number of domain Map elements owned by the calling process.
    virtual size_t getNodeNumCols() const = 0;

    //! The index base for global indices in this matrix.
    virtual GlobalOrdinal getIndexBase() const = 0;

    //! The global number of stored (structurally nonzero) entries.
    virtual global_size_t getGlobalNumEntries() const = 0;

    //! The local number of stored (structurally nonzero) entries.
    virtual size_t getNodeNumEntries() const = 0;

    /// \brief The current number of entries on the calling process in the specified global row.
    ///
    /// Note that if the row Map is overlapping, then the calling
    /// process might not necessarily store all the entries in the
    /// row.  Some other process might have the rest of the entries.
    ///
    /// \return <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt> if
    ///   the specified global row does not belong to this graph, else
    ///   the number of entries.
    virtual size_t getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const = 0;

    /// \brief The current number of entries on the calling process in the specified local row.
    ///
    /// Note that if the row Map is overlapping, then the calling
    /// process might not necessarily store all the entries in the
    /// row.  Some other process might have the rest of the entries.
    ///
    /// \return <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt> if
    ///   the specified local row is not valid for this graph, else
    ///   the number of entries.
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

    //! The number of global diagonal entries, based on global row/column index comparisons.
    virtual global_size_t getGlobalNumDiags() const = 0;

    //! The number of local diagonal entries, based on global row/column index comparisons.
    virtual size_t getNodeNumDiags() const = 0;

    //! The maximum number of entries across all rows/columns on all nodes.
    virtual size_t getGlobalMaxNumRowEntries() const = 0;

    //! The maximum number of entries across all rows/columns on this node.
    virtual size_t getNodeMaxNumRowEntries() const = 0;

    //! Whether this matrix has a well-defined column map.
    virtual bool hasColMap() const = 0;

    //! Whether this matrix is lower triangular.
    virtual bool isLowerTriangular() const = 0;

    //! Whether this matrix is upper triangular.
    virtual bool isUpperTriangular() const = 0;

    /// \brief Whether matrix indices are locally indexed.
    ///
    /// A RowMatrix may store column indices either as global indices
    /// (of type <tt>GlobalOrdinal</tt>), or as local indices (of type
    /// <tt>LocalOrdinal</tt>).  In some cases (for example, if the
    /// column Map has not been computed), it is not possible to
    /// switch from global to local indices without extra work.
    /// Furthermore, some operations only work for one or the other
    /// case.
    virtual bool isLocallyIndexed() const = 0;

    /// \brief Whether matrix indices are globally indexed.
    ///
    /// A RowMatrix may store column indices either as global indices
    /// (of type <tt>GlobalOrdinal</tt>), or as local indices (of type
    /// <tt>LocalOrdinal</tt>).  In some cases (for example, if the
    /// column Map has not been computed), it is not possible to
    /// switch from global to local indices without extra work.
    /// Furthermore, some operations only work for one or the other
    /// case.
    virtual bool isGloballyIndexed() const = 0;

    //! Whether fillComplete() has been called.
    virtual bool isFillComplete() const = 0;

    //! Whether this object implements getLocalRowView() and getGlobalRowView().
    virtual bool supportsRowViews() const = 0;

    //@}
    //! @name Extraction Methods
    //@{

    /// \brief Get a copy of the given global row's entries.
    ///
    /// This method only gets the entries in the given row that are
    /// stored on the calling process.  Note that if the matrix has an
    /// overlapping row Map, it is possible that the calling process
    /// does not store all the entries in that row.
    ///
    /// \param GlobalRow [in] Global index of the row.
    /// \param Indices [out] Global indices of the columns
    ///   corresponding to values.
    /// \param Values [out] Matrix values.
    /// \param NumEntries [out] Number of stored entries on the
    ///   calling process; length of Indices and Values.
    ///
    /// This method throws <tt>std::runtime_error</tt> if either
    /// Indices or Values is not large enough to hold the data
    /// associated with row GlobalRow. If GlobalRow does not belong to
    /// the calling process, then the method sets NumIndices to
    /// <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>, and does
    /// not modify Indices or Values.
    virtual void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      const Teuchos::ArrayView<GlobalOrdinal> &Indices,
                      const Teuchos::ArrayView<Scalar> &Values,
                      size_t &NumEntries) const = 0;

    /// \brief Get a copy of the given local row's entries.
    ///
    /// This method only gets the entries in the given row that are
    /// stored on the calling process.  Note that if the matrix has an
    /// overlapping row Map, it is possible that the calling process
    /// does not store all the entries in that row.
    ///
    /// \param LocalRow [in] Local index of the row.
    /// \param Indices [out] Local indices of the columns
    ///   corresponding to values.
    /// \param Values [out] Matrix values.
    /// \param NumEntries [out] Number of stored entries on the
    ///   calling process; length of Indices and Values.
    ///
    /// This method throws <tt>std::runtime_error</tt> if either
    /// Indices or Values is not large enough to hold the data
    /// associated with row LocalRow. If LocalRow does not belong to
    /// the calling process, then the method sets NumIndices to
    /// <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>, and does
    /// not modify Indices or Values.
    virtual void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     const Teuchos::ArrayView<LocalOrdinal> &Indices,
                     const Teuchos::ArrayView<Scalar> &Values,
                     size_t &NumEntries) const = 0;

    /// \brief Get a constant, nonpersisting, globally indexed view of
    ///   the given row of the matrix.
    ///
    /// The returned views of the column indices and values are not
    /// guaranteed to persist beyond the lifetime of <tt>this</tt>.
    /// Furthermore, some RowMatrix implementations allow changing the
    /// values, or the indices and values.  Any such changes
    /// invalidate the returned views.
    ///
    /// This method only gets the entries in the given row that are
    /// stored on the calling process.  Note that if the matrix has an
    /// overlapping row Map, it is possible that the calling process
    /// does not store all the entries in that row.
    ///
    /// \pre <tt>isGloballyIndexed () && supportsRowViews ()</tt>
    /// \post <tt>indices.size () == getNumEntriesInGlobalRow (GlobalRow)</tt>
    ///
    /// \param GlobalRow [in] Global index of the row.
    /// \param Indices [out] Global indices of the columns
    ///   corresponding to values.
    /// \param Values [out] Matrix values.
    ///
    /// If \c GlobalRow does not belong to this node, then \c indices
    /// is set to \c null.
    virtual void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      ArrayView<const GlobalOrdinal> &indices,
                      ArrayView<const Scalar> &values) const = 0;

    /// \brief Get a constant, nonpersisting, locally indexed view of
    ///   the given row of the matrix.
    ///
    /// The returned views of the column indices and values are not
    /// guaranteed to persist beyond the lifetime of <tt>this</tt>.
    /// Furthermore, some RowMatrix implementations allow changing the
    /// values, or the indices and values.  Any such changes
    /// invalidate the returned views.
    ///
    /// This method only gets the entries in the given row that are
    /// stored on the calling process.  Note that if the matrix has an
    /// overlapping row Map, it is possible that the calling process
    /// does not store all the entries in that row.
    ///
    /// \pre <tt>isLocallyIndexed () && supportsRowViews ()</tt>
    /// \post <tt>indices.size () == getNumEntriesInGlobalRow (LocalRow)</tt>
    ///
    /// \param LocalRow [in] Local index of the row.
    /// \param Indices [out] Local indices of the columns
    ///   corresponding to values.
    /// \param Values [out] Matrix values.
    ///
    /// If \c LocalRow does not belong to this node, then \c indices
    /// is set to \c null.
    virtual void
    getLocalRowView (LocalOrdinal LocalRow,
                     ArrayView<const LocalOrdinal> &indices,
                     ArrayView<const Scalar> &values) const = 0;

    /// \brief Get a copy of the diagonal entries, distributed by the row Map.
    ///
    /// On input, the Vector's Map must be the same as the row Map of the matrix.
    /// (That is, <tt>this->getRowMap ()->isSameAs (* (diag.getMap ())) == true</tt>.)
    ///
    /// On return, the entries of \c diag are filled with the diagonal
    /// entries of the matrix stored on this process.  Note that if
    /// the row Map is overlapping, multiple processes may own the
    /// same diagonal element.  You may combine these overlapping
    /// diagonal elements by doing an Export from the row Map Vector
    /// to a range Map Vector.
    virtual void getLocalDiagCopy (Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const = 0;

    //@}
    //! \name Mathematical methods
    //@{

    /**
     * \brief Scale the RowMatrix on the left with the given Vector x.
     *
     * On return, for all entries i,j in the matrix, \f$A(i,j) = x(i)*A(i,j)\f$.
     */
    virtual void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) = 0;

    /**
     * \brief Scale the RowMatrix on the right with the given Vector x.
     *
     * On return, for all entries i,j in the matrix, \f$A(i,j) = x(j)*A(i,j)\f$.
     */
    virtual void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) = 0;

    /// \brief The Frobenius norm of the matrix.
    ///
    /// This method computes and returns the Frobenius norm of the
    /// matrix.  The Frobenius norm \f$\|A\|_F\f$ for the matrix
    /// \f$A\f$ is defined as
    /// \f$\|A\|_F = \sqrt{ \sum_{i,j} |\a_{ij}|^2 }\f$.
    /// It has the same value as the Euclidean norm of a vector made
    /// by stacking the columns of \f$A\f$.
    virtual typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const = 0;
    //@}
  }; // class RowMatrix

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~RowMatrix() {
  }

} // namespace Tpetra

#endif
