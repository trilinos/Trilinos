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

#ifndef TPETRA_ROWMATRIX_DECL_HPP
#define TPETRA_ROWMATRIX_DECL_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_RowGraph_fwd.hpp"
#include "Tpetra_Packable.hpp"
#include "Tpetra_SrcDistObject.hpp"
#include "Teuchos_Describable.hpp"

namespace Tpetra {
namespace Classes {
  /// \class RowMatrix
  /// \brief A read-only, row-oriented interface to a sparse matrix.
  ///
  /// \tparam Scalar The type of the entries in the sparse matrix.
  ///   Same as the \c Scalar typedef of Operator.
  /// \tparam LocalOrdinal The type of local indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam GlobalOrdinal The type of global indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam Node The Kokkos Node type.  See the documentation of Map
  ///   for requirements.
  ///
  /// RowMatrix provides a read-only, row-oriented interface to view
  /// the entries of a sparse matrix, which is distributed over one or
  /// more (MPI) processes.  "Read only" means that you may view the
  /// entries, but not modify them.  "Row-oriented" means that you may
  /// view all the entries owned by the calling process, in any
  /// particular row of the matrix that is owned by the calling
  /// process.
  ///
  /// CrsMatrix implements this interface, but also lets you add or
  /// modify entries of the sparse matrix.  Ifpack2 provides other
  /// implementations of RowMatrix, which do useful things like
  /// wrapping an existing matrix to view only certain desired
  /// entries.
  template <class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
            class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
            class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
            class Node = ::Tpetra::Details::DefaultTypes::node_type>
  class RowMatrix :
    virtual public Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
    virtual public SrcDistObject,
    public Packable<char, LocalOrdinal> {
  public:
    //! @name Typedefs
    //@{

    //! The type of the entries in the sparse matrix.
    typedef Scalar        scalar_type;
    //! The type of local indices.
    typedef LocalOrdinal  local_ordinal_type;
    //! The type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! The Kokkos Node type.
    typedef Node          node_type;

    /// \brief Type of a norm result.
    ///
    /// This is usually the same as the type of the magnitude
    /// (absolute value) of <tt>Scalar</tt>, but may differ for
    /// certain <tt>Scalar</tt> types.
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type mag_type;

    //@}
    //! @name Destructor
    //@{

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~RowMatrix();

    //@}
    //! @name Matrix query methods
    //@{

    //! The communicator over which this matrix is distributed.
    virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const = 0;

    //! The Kokkos Node instance.
    virtual Teuchos::RCP<Node> getNode() const = 0;

    //! The Map that describes the distribution of rows over processes.
    virtual Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const = 0;

    //! The Map that describes the distribution of columns over processes.
    virtual Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const = 0;

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

    /// \brief Number of diagonal entries in the matrix's graph, over
    ///   all processes in the matrix's communicator.
    ///
    /// \warning DO NOT CALL THIS METHOD!  This method is DEPRECATED
    ///   and will DISAPPEAR VERY SOON per #2630.
    ///
    /// \pre Subclasses reserve the right to impose preconditions on
    ///   the matrix's state.
    virtual global_size_t TPETRA_DEPRECATED getGlobalNumDiags () const = 0;

    /// \brief Number of diagonal entries in the matrix's graph, on
    ///   the calling process.
    ///
    /// \warning DO NOT CALL THIS METHOD!  This method is DEPRECATED
    ///   and will DISAPPEAR VERY SOON per #2630.
    ///
    /// \pre Subclasses reserve the right to impose preconditions on
    ///   the matrix's state.
    virtual size_t TPETRA_DEPRECATED getNodeNumDiags () const = 0;

    /// \brief Maximum number of entries in any row of the matrix,
    ///   over all processes.
    ///
    /// \pre Subclasses reserve the right to impose preconditions on
    ///   the matrix's state.
    ///
    /// This method only uses the matrix's graph.  Explicitly stored
    /// zeros count as "entries."
    virtual size_t getGlobalMaxNumRowEntries () const = 0;

    /// \brief Maximum number of entries in any row of the matrix,
    ///   on this process.
    ///
    /// \pre Subclasses reserve the right to impose preconditions on
    ///   the matrix's state.
    ///
    /// This method only uses the matrix's graph.  Explicitly stored
    /// zeros count as "entries."
    virtual size_t getNodeMaxNumRowEntries () const = 0;

    //! Whether this matrix has a well-defined column Map.
    virtual bool hasColMap () const = 0;

    /// \brief Whether the matrix is locally lower triangular.
    ///
    /// \warning DO NOT CALL THIS METHOD!  This method is DEPRECATED
    ///   and will DISAPPEAR VERY SOON per #2630.
    ///
    /// \pre Subclasses reserve the right to impose preconditions on
    ///   the matrix's state.
    ///
    /// \note This is entirely a local property.  That means this
    ///   method may return different results on different processes.
    virtual bool TPETRA_DEPRECATED isLowerTriangular () const = 0;

    /// \brief Whether the matrix is locally upper triangular.
    ///
    /// \warning DO NOT CALL THIS METHOD!  This method is DEPRECATED
    ///   and will DISAPPEAR VERY SOON per #2630.
    ///
    /// \pre Subclasses reserve the right to impose preconditions on
    ///   the matrix's state.
    ///
    /// \note This is entirely a local property.  That means this
    ///   method may return different results on different processes.
    virtual bool TPETRA_DEPRECATED isUpperTriangular () const = 0;

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
                      Teuchos::ArrayView<const GlobalOrdinal> &indices,
                      Teuchos::ArrayView<const Scalar> &values) const = 0;

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
                     Teuchos::ArrayView<const LocalOrdinal>& indices,
                     Teuchos::ArrayView<const Scalar>& values) const = 0;

    /// \brief Get a constant, nonpersisting, locally indexed view of
    ///   the given row of the matrix, using "raw" pointers instead of
    ///   Teuchos::ArrayView.
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
    /// \post <tt>numEnt == getNumEntriesInGlobalRow (LocalRow)</tt>
    ///
    /// \param lclRow [in] Local index of the row.
    /// \param numEnt [out] Number of entries in the row that are
    ///   stored on the calling process.
    /// \param lclColInds [out] Local indices of the columns
    ///   corresponding to values.
    /// \param vals [out] Matrix values.
    ///
    /// \return Error code; zero on no error.
    virtual LocalOrdinal
    getLocalRowViewRaw (const LocalOrdinal lclRow,
                        LocalOrdinal& numEnt,
                        const LocalOrdinal*& lclColInds,
                        const Scalar*& vals) const;

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
    virtual void getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) const = 0;

    //@}
    //! \name Mathematical methods
    //@{

    /// \brief Scale the matrix on the left with the given Vector.
    ///
    /// On return, for all entries i,j in the matrix,
    /// \f$A(i,j) = x(i)*A(i,j)\f$.
    virtual void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) = 0;

    /// \brief Scale the matrix on the right with the given Vector.
    ///
    /// On return, for all entries i,j in the matrix,
    /// \f$A(i,j) = x(j)*A(i,j)\f$.
    virtual void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) = 0;

    /// \brief The Frobenius norm of the matrix.
    ///
    /// This method computes and returns the Frobenius norm of the
    /// matrix.  The Frobenius norm \f$\|A\|_F\f$ for the matrix
    /// \f$A\f$ is defined as
    /// \f$\|A\|_F = \sqrt{ \sum_{i,j} |\a_{ij}|^2 }\f$.
    /// It has the same value as the Euclidean norm of a vector made
    /// by stacking the columns of \f$A\f$.
    virtual mag_type getFrobeniusNorm() const = 0;

    /// \brief Return a new RowMatrix which is the result of <tt>beta*this + alpha*A</tt>.
    ///
    /// The new RowMatrix is actually a CrsMatrix (which see).  Note
    /// that RowMatrix is a read-only interface (not counting the left
    /// and right scale methods), so it is impossible to implement an
    /// in-place add using just that interface.
    ///
    /// For brevity, call this matrix B, and the result matrix C.  C's
    /// row Map will be identical to B's row Map.  It is correct,
    /// though less efficient, for A and B not to have the same row
    /// Maps.  We could make C's row Map the union of the two row Maps
    /// in that case.  However, we don't want row Maps to grow for a
    /// repeated sequence of additions with matrices with different
    /// row Maps.  Furthermore, the fact that the user called this
    /// method on B, rather than on A, suggests a preference for using
    /// B's distribution.  The most reasonable thing to do, then, is
    /// to use B's row Map for C.
    ///
    /// A and B must have identical or congruent communicators.  This
    /// method must be called as a collective over B's communicator.
    ///
    /// The parameters are optional and may be null.  Here are the
    /// parameters that this function accepts:
    ///   - "Call fillComplete" (\c bool): If true, call fillComplete on
    ///     the result matrix C.  This is \c true by default.
    ///   - "Constructor parameters" (sublist): If provided, give these
    ///     parameters to C's constructor.
    ///   - "fillComplete parameters" (sublist): If provided, and if
    ///     "Call fillComplete" is true, then give these parameters to
    ///     C's fillComplete call.
    ///
    /// It is not strictly necessary that a RowMatrix always have a
    /// domain and range Map.  For example, a CrsMatrix does not have
    /// a domain and range Map until after its first fillComplete
    /// call.  Neither A nor B need to have a domain and range Map in
    /// order to call add().  If at least one of them has a domain and
    /// range Map, you need not supply a domain and range Map to this
    /// method.  If you ask this method to call fillComplete on C (it
    /// does by default), it will supply the any missing domain or
    /// range Maps from either B's or A's (in that order) domain and
    /// range Maps.  If neither A nor B have a domain and range Map,
    /// and if you ask this method to call fillComplete, then you
    /// <i>must</i> supply both a domain Map and a range Map to this
    /// method.
    ///
    /// This method comes with a default implementation, since the
    /// RowMatrix interface suffices for implementing it.  Subclasses
    /// (like CrsMatrix) may override this implementation, for example
    /// to improve its performance, given additional knowledge about
    /// the subclass.  Subclass implementations may need to do a
    /// dynamic cast on A in order to know its type.
    virtual Teuchos::RCP<RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    add (const Scalar& alpha,
         const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Scalar& beta,
         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap = Teuchos::null,
         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap = Teuchos::null,
         const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;
    //@}
    //! \name Implementation of Packable interface
    //@{
  private:
    bool
    packRow (char* const numEntOut,
             char* const valOut,
             char* const indOut,
             const size_t numEnt,
             const LocalOrdinal lclRow) const;

    // TODO (mfh 25 Jan 2015) Could just make this "protected" and let
    // CrsMatrix use it, since it's exactly the same there.
    void
    allocatePackSpace (Teuchos::Array<char>& exports,
                       size_t& totalNumEntries,
                       const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs) const;

    /// \brief Pack this object's data for an Import or Export.
    ///
    /// \warning To be called only by the default implementation of
    ///   pack() (see below).
    void
    packImpl (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
              Teuchos::Array<char>& exports,
              const Teuchos::ArrayView<size_t>& numPacketsPerLID,
              size_t& constantNumPackets,
              Distributor& distor) const;


  public:
    /// \brief Pack this object's data for an Import or Export.
    ///
    /// \warning To be called only by the packAndPrepare method of
    ///   appropriate classes of DistObject.
    ///
    /// Subclasses may override this method to speed up or otherwise
    /// improve the implementation by exploiting more specific details
    /// of the subclass.
    virtual void
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
          Teuchos::Array<char>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets,
          Distributor& distor) const;
    //@}
  }; // class RowMatrix
} // namespace Classes
} // namespace Tpetra

#endif // TPETRA_ROWMATRIX_DECL_HPP

