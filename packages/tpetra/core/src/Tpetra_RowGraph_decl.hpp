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

#ifndef TPETRA_ROWGRAPH_DECL_HPP
#define TPETRA_ROWGRAPH_DECL_HPP

#include "Tpetra_RowGraph_fwd.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Packable.hpp"
#include "Teuchos_Describable.hpp"

namespace Tpetra {

  /// \class RowGraph
  /// \brief An abstract interface for graphs accessed by rows.
  ///
  /// This class is to CrsGraph, what RowMatrix is to CrsMatrix.
  /// CrsGraph is an implementation of RowGraph.
  ///
  /// \tparam LocalOrdinal The type of local indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam GlobalOrdinal The type of global indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam Node The Kokkos Node type.  See the documentation of Map
  ///   for requirements.
  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class RowGraph :
    virtual public Teuchos::Describable,
    public Packable<GlobalOrdinal, LocalOrdinal> {
  public:
    //! \name Typedefs
    //@{
    //! The type of local indices in the graph.
    typedef LocalOrdinal  local_ordinal_type;
    //! The type of global indices in the graph.
    typedef GlobalOrdinal global_ordinal_type;
    //! The Kokkos Node type.
    typedef Node          node_type;
    //@}

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~RowGraph() {};

    //! @name Graph query methods
    //@{

    //! The communicator over which this graph is distributed.
    virtual Teuchos::RCP<const Teuchos::Comm<int> >
    getComm () const = 0;

    //! The Kokkos Node instance with which this object was created.
    virtual Teuchos::RCP<Node> getNode () const = 0;

    //! The Map that describes this graph's distribution of rows over processes.
    virtual Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
    getRowMap () const = 0;

    //! The Map that describes this graph's distribution of columns over processes.
    virtual Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
    getColMap () const = 0;

    //! The Map associated with the domain of this graph.
    virtual Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
    getDomainMap () const = 0;

    //! The Map associated with the range of this graph.
    virtual Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
    getRangeMap () const = 0;

    //! This graph's Import object.
    virtual Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >
    getImporter () const = 0;

    //! This graph's Export object.
    virtual Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> >
    getExporter () const = 0;

    //! Returns the number of global rows in the graph.
    virtual global_size_t getGlobalNumRows() const = 0;

    //! \brief Returns the number of global columns in the graph.
    virtual global_size_t getGlobalNumCols() const = 0;

    //! Returns the number of rows owned on the calling node.
    virtual size_t getNodeNumRows() const = 0;

    //! Returns the number of columns connected to the locally owned rows of this graph.
    virtual size_t getNodeNumCols() const = 0;

    //! Returns the index base for global indices for this graph.
    virtual GlobalOrdinal getIndexBase() const = 0;

    //! Returns the global number of entries in the graph.
    virtual global_size_t getGlobalNumEntries() const = 0;

    //! Returns the local number of entries in the graph.
    virtual size_t getNodeNumEntries() const = 0;

    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
    virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

    /// \brief Number of diagonal entries over all processes in the
    ///   graph's communicator.
    ///
    /// \warning DO NOT CALL THIS METHOD!  This method is DEPRECATED
    ///   and will DISAPPEAR VERY SOON per #2630.
    virtual global_size_t TPETRA_DEPRECATED getGlobalNumDiags () const = 0;

    /// \brief Number of diagonal entries on the calling process.
    ///
    /// \warning DO NOT CALL THIS METHOD!  This method is DEPRECATED
    ///   and will DISAPPEAR VERY SOON per #2630.
    virtual size_t TPETRA_DEPRECATED getNodeNumDiags () const = 0;

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    virtual size_t getGlobalMaxNumRowEntries() const = 0;

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    virtual size_t getNodeMaxNumRowEntries() const = 0;

    //! Whether the graph has a well-defined column Map.
    virtual bool hasColMap() const = 0;

    /// \brief Whether the graph is locally lower triangular.
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

    /// \brief Whether the graph is locally upper triangular.
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

    //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
    virtual bool isLocallyIndexed() const = 0;

    //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
    virtual bool isGloballyIndexed() const = 0;

    //! Whether fillComplete() has been called (without an intervening resumeFill()).
    virtual bool isFillComplete() const = 0;

    //@}
    //! @name Extraction Methods
    //@{

    //! Extract a list of entries in a specified global row of the graph. Put into pre-allocated storage.
    /*!
      \param LocalRow - (In) Global row number for which indices are desired.
      \param Indices - (Out) Global column indices corresponding to values.
      \param NumIndices - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown if \c Indices is not large enough to hold the column indices associated
      with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices is unchanged and \c NumIndices is
      returned as Teuchos::OrdinalTraits<size_t>::invalid().
    */
    virtual void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      const Teuchos::ArrayView<GlobalOrdinal> &Indices,
                      size_t &NumIndices) const = 0;

    //! Extract a list of entries in a specified local row of the graph. Put into storage allocated by calling routine.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices - (Out) Local column indices corresponding to values.
      \param NumIndices - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown if \c Indices is not large enough to hold the column indices associated
      with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices is unchanged and \c NumIndices is
      returned as Teuchos::OrdinalTraits<size_t>::invalid().
    */
    virtual void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     const Teuchos::ArrayView<LocalOrdinal> &Indices,
                     size_t &NumIndices) const = 0;

    /// \brief Whether this class implements getLocalRowView() and
    ///   getGlobalRowView().
    ///
    /// If subclasses override the default (trivial) implementation of
    /// getLocalRowView() and getGlobalRowView(), then they need to
    /// override this method as well.
    virtual bool supportsRowViews () const {
      return false;
    }

    /// \brief Get a constant, nonpersisting, locally indexed view of
    ///   the given row of the graph.
    ///
    /// The returned views of the column indices are not guaranteed to
    /// persist beyond the lifetime of <tt>this</tt>.  Furthermore,
    /// some RowGraph implementations allow changing the values, or
    /// the indices and values.  Any such changes invalidate the
    /// returned views.
    ///
    /// This method only gets the entries in the given row that are
    /// stored on the calling process.  Note that if the graph has an
    /// overlapping row Map, it is possible that the calling process
    /// does not store all the entries in that row.
    ///
    /// \pre <tt>isLocallyIndexed () && supportsRowViews ()</tt>
    /// \post <tt>indices.size () == getNumEntriesInGlobalRow (LocalRow)</tt>
    ///
    /// \param lclRow [in] Local index of the row.
    /// \param lclColInds [out] Local indices of the columns in the
    ///   row.  If the given row is not a valid row index on the
    ///   calling process, then the result has no entries (its size is
    ///   zero).
    ///
    /// Subclasses are expected to implement this method.  We would
    /// have made this method pure virtual, but that would have broken
    /// backwards compatibility, since we added the method at least
    /// one major release after introducing this class.
    virtual void
    getLocalRowView (const LocalOrdinal lclRow,
                     Teuchos::ArrayView<const LocalOrdinal>& lclColInds) const;

    /// \brief Get a const, non-persisting view of the given global
    ///   row's global column indices, as a Teuchos::ArrayView.
    ///
    /// \param gblRow [in] Global index of the row.
    /// \param gblColInds [out] Global column indices in the row.  If
    ///   the given row is not a valid row index on the calling
    ///   process, then the result has no entries (its size is zero).
    ///
    /// \pre <tt>! isLocallyIndexed()</tt>
    /// \post <tt>gblColInds.size() == getNumEntriesInGlobalRow(gblRow)</tt>
    ///
    /// Subclasses are expected to implement this method.  We would
    /// have made this method pure virtual, but that would have broken
    /// backwards compatibility, since we added the method at least
    /// one major release after introducing this class.
    virtual void
    getGlobalRowView (const GlobalOrdinal gblRow,
                      Teuchos::ArrayView<const GlobalOrdinal>& gblColInds) const;

    //@}
    //! \name Implementation of Packable interface
    //@{

    //! Pack this object's data for Import or Export.
    virtual void
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
          Teuchos::Array<GlobalOrdinal>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets,
          Distributor& distor) const;
    //@}
  }; // class RowGraph
} // namespace Tpetra

#endif // TPETRA_ROWGRAPH_DECL_HPP
