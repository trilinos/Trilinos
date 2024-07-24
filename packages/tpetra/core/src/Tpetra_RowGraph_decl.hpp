// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    typedef LocalOrdinal local_ordinal_type;
    //! The type of global indices in the graph.
    typedef GlobalOrdinal global_ordinal_type;
    //! The Kokkos Node type.
    typedef Node node_type;
    //@}

    typedef typename
        Kokkos::View<LocalOrdinal *, typename Node::device_type>::const_type
        local_inds_device_view_type;
    typedef typename local_inds_device_view_type::HostMirror::const_type
        local_inds_host_view_type;
    typedef typename local_inds_device_view_type::HostMirror
        nonconst_local_inds_host_view_type;


    typedef typename
        Kokkos::View<GlobalOrdinal *, typename Node::device_type>::const_type
        global_inds_device_view_type;
    typedef typename global_inds_device_view_type::HostMirror::const_type
        global_inds_host_view_type;
    typedef typename global_inds_device_view_type::HostMirror
        nonconst_global_inds_host_view_type;

    typedef typename 
        Kokkos::View<const size_t*, typename Node::device_type>::const_type
        row_ptrs_device_view_type;
    typedef typename row_ptrs_device_view_type::HostMirror::const_type
        row_ptrs_host_view_type;

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~RowGraph() {};

    //! @name Graph query methods
    //@{

    //! The communicator over which this graph is distributed.
    virtual Teuchos::RCP<const Teuchos::Comm<int> >
    getComm () const = 0;


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
    virtual size_t getLocalNumRows() const = 0;

    //! Returns the number of columns connected to the locally owned rows of this graph.
    virtual size_t getLocalNumCols() const = 0;

    //! Returns the index base for global indices for this graph.
    virtual GlobalOrdinal getIndexBase() const = 0;

    //! Returns the global number of entries in the graph.
    virtual global_size_t getGlobalNumEntries() const = 0;

    //! Returns the local number of entries in the graph.
    virtual size_t getLocalNumEntries() const = 0;

    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
    virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    virtual size_t getGlobalMaxNumRowEntries() const = 0;

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    virtual size_t getLocalMaxNumRowEntries() const = 0;

    //! Whether the graph has a well-defined column Map.
    virtual bool hasColMap() const = 0;

    //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
    virtual bool isLocallyIndexed() const = 0;

    //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
    virtual bool isGloballyIndexed() const = 0;

    //! Whether fillComplete() has been called (without an intervening resumeFill()).
    virtual bool isFillComplete() const = 0;


    //@}
    //! @name Access to entries in a row
    //@{

    /// \brief Get a copy of the global column indices in a given row
    ///   of the graph.
    ///
    /// Given the global index of a row of the graph, get a copy of
    /// all the global column indices in that row that the calling
    /// process stores.
    ///
    /// \param gblRow [in] Global index of the row.
    /// \param gblColInds [in/out] On output: All the global column
    ///   indices in that row on the calling process.
    /// \param numColInds [out] Number of indices in the row on the
    ///   calling process.
    ///
    /// \pre <tt>getRowMap()->isNodeGlobalElement(gblRow)<tt> is <tt>true</tt>.
    /// \pre <tt>gblColInds.size() >= getNumEntriesInGlobalRow(gblRow)</tt> is <tt>true</tt>.

    virtual void
    getGlobalRowCopy (const GlobalOrdinal gblRow,
                      nonconst_global_inds_host_view_type& gblColInds,
                      size_t& numColInds) const = 0;

    /// \brief Get a copy of the local column indices in a given row
    ///   of the graph.
    ///
    /// Given the local index of a row of the graph, get a copy of
    /// all the local column indices in that row that the calling
    /// process stores.
    ///
    /// \param lclRow [in] Local index of the row.
    /// \param lclColInds [in/out] On output: All the local column
    ///   indices in that row on the calling process.
    /// \param numColInds [out] Number of indices in the row on the
    ///   calling process.
    ///
    /// \pre <tt>hasColMap()</tt> is <tt>true</tt>.
    /// \pre <tt>getRowMap()->isNodeLocalElement(lclRow)<tt> is <tt>true</tt>.
    /// \pre <tt>lclColInds.size() >= getNumEntriesInLocalRow(lclRow)</tt> is <tt>true</tt>.
    virtual void
    getLocalRowCopy (const LocalOrdinal lclRow,
                     nonconst_local_inds_host_view_type & lclColInds,
                     size_t& numColInds) const = 0;

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
                     local_inds_host_view_type & lclColInds) const = 0;

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
                      global_inds_host_view_type& gblColInds) const = 0;

    //@}
    //! \name Implementation of Packable interface
    //@{

    //! Pack this object's data for Import or Export.
    virtual void
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
          Teuchos::Array<GlobalOrdinal>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets) const;
    //@}
  }; // class RowGraph
} // namespace Tpetra

#endif // TPETRA_ROWGRAPH_DECL_HPP
