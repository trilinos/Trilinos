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

#ifndef TPETRANEW_ROWGRAPH_HPP
#define TPETRANEW_ROWGRAPH_HPP

#include "TpetraNew_RowGraph_fwd.hpp"
#include "TpetraNew_Map_fwd.hpp"
#include "TpetraNew_Import_fwd.hpp"
#include "TpetraNew_Export_fwd.hpp"
#include "Tpetra_Packable.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos {
  template<class T> class ArrayView; // forward declaration
  template<class T> class Array; // forward declaration  
  template<class OrdinalType> class Comm; // forward declaration
} // namespace Teuchos

namespace Tpetra {
  class Distributor; // forward declaration
} // namespace Tpetra

namespace TpetraNew {

  /// \typedef global_size_t
  ///
  /// The maximum global number of entries that a RowGraph or
  /// RowMatrix may have.  This typedef exists in case an individual
  /// MPI process has 32-bit size_t.
  using global_size_t = unsigned long long;

  /// \class RowGraph
  /// \brief An abstract interface for graphs accessed by rows.
  ///
  /// This class is to CrsGraph, what RowMatrix is to CrsMatrix.
  /// CrsGraph is an implementation of RowGraph.
  class RowGraph :
    virtual public Teuchos::Describable,
    public ::Tpetra::Packable<
      ::Tpetra::Details::DefaultTypes::global_ordinal_type,
      ::Tpetra::Details::DefaultTypes::local_ordinal_type>
  {
  public:
    //! \name Typedefs
    //@{

    //! The type of local indices.
    using local_ordinal_type =
      ::Tpetra::Details::DefaultTypes::local_ordinal_type;
    //! The type of global indices.
    using global_ordinal_type =
      ::Tpetra::Details::DefaultTypes::global_ordinal_type;
    //! The Kokkos execution space.
    using execution_space =
      ::Tpetra::Details::DefaultTypes::execution_space;
    //! The Kokkos memory space.
    using memory_space = ::Tpetra::Details::DefaultTypes::memory_space;
    //! The Kokkos device type.
    using device_type = ::Tpetra::Details::DefaultTypes::device_type;

    //@}

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~RowGraph() {};

    //! @name Graph query methods
    //@{

    //! The communicator over which this graph is distributed.
    virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm () const = 0;

    //! The Map that describes this graph's distribution of rows over processes.
    virtual Teuchos::RCP<const Map> getRowMap () const = 0;

    //! The Map that describes this graph's distribution of columns over processes.
    virtual Teuchos::RCP<const Map> getColMap () const = 0;

    //! The Map associated with the domain of this graph.
    virtual Teuchos::RCP<const Map> getDomainMap () const = 0;

    //! The Map associated with the range of this graph.
    virtual Teuchos::RCP<const Map> getRangeMap () const = 0;

    //! This graph's Import object.
    virtual Teuchos::RCP<const Import> getImporter () const = 0;

    //! This graph's Export object.
    virtual Teuchos::RCP<const Export> getExporter () const = 0;

    /// \brief The global number of rows in the graph.
    ///
    /// \note This used to return Tpetra::global_size_t, generally an
    ///   unsigned type.  Now it returns global_ordinal_type, which
    ///   may be signed.
    virtual global_ordinal_type getGlobalNumRows () const = 0;

    /// \brief The global number of columns in the graph.
    ///
    /// \note This used to return Tpetra::global_size_t, generally an
    ///   unsigned type.  Now it returns global_ordinal_type, which
    ///   may be signed.
    virtual global_ordinal_type getGlobalNumCols () const = 0;

    /// \brief The number of rows that live on the calling process.
    ///
    /// \note This used to return size_t, generally an unsigned type.
    ///   Now it returns local_ordinal_type, which may be signed.
    virtual local_ordinal_type getMyNumRows () const = 0;

    /// \brief The number of rows that live on the calling process.    
    ///
    /// \warning Don't call this; call getMyNumRows() instead.
    virtual local_ordinal_type TPETRA_DEPRECATED
    getNodeNumRows () const = 0;

    /// \brief The number of columns connected to the locally owned
    ///   rows of this graph.
    ///
    /// \note This used to return size_t, generally an unsigned type.
    ///   Now it returns local_ordinal_type, which may be signed.
    virtual local_ordinal_type getMyNumCols () const = 0;    

    /// \brief The number of columns connected to the locally owned rows of this graph.
    ///
    /// \warning Don't call this; call getMyNumCols() instead.
    virtual local_ordinal_type TPETRA_DEPRECATED
    getNodeNumCols () const = 0;

    //! Returns the index base for global indices for this graph.
    virtual global_ordinal_type getIndexBase() const = 0;

    //! The global number of entries in the graph.
    virtual global_size_t getGlobalNumEntries () const = 0;

    //! The number of entries in the graph on the calling process.
    virtual size_t getMyNumEntries () const = 0;

    /// \brief The local number of entries in the graph.
    ///
    /// \warning Don't call this; call getMyNumEntries() instead.
    virtual size_t TPETRA_DEPRECATED getNodeNumEntries () const = 0;

    /// \brief Current number of entries on this process in the row,
    ///   specified by global index.
    ///
    /// \return Teuchos::OrdinalTraits<local_ordinal_type>::invalid()
    ///   if the specified row does not live on this process, else the
    ///   number of entries.
    virtual local_ordinal_type
    getNumEntriesInGlobalRow (const global_ordinal_type globalRow) const = 0;

    /// \brief Current number of entries on this process in the row,
    ///   specified by local index.
    ///
    /// \return Teuchos::OrdinalTraits<local_ordinal_type>::invalid()
    ///   if the specified row does not live on this process, else the
    ///   number of entries.
    virtual local_ordinal_type
    getNumEntriesInLocalRow (const local_ordinal_type localRow) const = 0;

    /// \brief The maximum number of graph entries in any row of the
    ///   graph, over all processes in the graph's communicator.
    ///
    /// This returns local_ordinal_type because no row of the graph
    ///   may have more than local_ordinal_type entries.  (This
    ///   follows from the definition of column Map.)
    virtual global_ordinal_type getGlobalMaxNumRowEntries () const = 0;

    /// \brief The maximum number of entries in any row of the graph,
    ///   on this process.
    virtual local_ordinal_type getMyMaxNumRowEntries () const = 0;

    /// \brief The maximum number of entries in any row of the graph,
    ///   on this process.
    ///
    /// \warning Don't call this; call getMyMaxNumRowEntries() instead.
    virtual local_ordinal_type TPETRA_DEPRECATED
    getNodeMaxNumRowEntries () const = 0;

    //! Whether the graph has a well-defined column Map.
    virtual bool hasColMap() const = 0;

    //! Whether the graph stores its column indices as local indices.
    virtual bool isLocallyIndexed () const = 0;

    //! Whether the graph stores its column indices as global indices.
    virtual bool isGloballyIndexed () const = 0;

    /// \brief Whether fillComplete() has been called (without an
    ///   intervening resumeFill()).
    ///
    /// See documentation of CrsMatrix::fillComplete.
    virtual bool isFillComplete () const = 0;

    //@}
    //! @name Extraction Methods
    //@{

    /// \brief Given a global index of a row, copy all column indices
    ///   in that row that live on the calling process into
    ///   user-provided storage.  Store the column indices as global
    ///   indices, regardless of how the graph stores them.
    ///
    /// \param gblRow [in] Global index of row.
    ///
    /// \param gblColInds [out] All column indices in the row on the
    ///   calling process.  If \c gblColInds is not long enough to
    ///   hold all the indices, then the contents of this array are
    ///   unspecified.
    ///
    /// \return The number of column indices in the row on the calling
    ///   process.  If \c gblRow does not live on the calling process,
    ///   then return zero (as there are no column indices in this row
    ///   on the calling process).
    virtual local_ordinal_type
    getGlobalRowCopy (const global_ordinal_type gblRow,
                      const Teuchos::ArrayView<global_ordinal_type>& gblColInds) const = 0;

    /// \brief Given a local index of a row, copy all column indices
    ///   in that row that live on the calling process into
    ///   user-provided storage.  Store the column indices as local
    ///   indices, regardless of how the graph stores them.
    ///
    /// \pre <tt>hasColMap()</tt> is <tt>true</tt>.
    ///
    /// \param gblRow [in] Global index of row.
    ///
    /// \param lclColInds [out] All column indices in the row on the
    ///   calling process.  If \c lclColInds is not long enough to
    ///   hold all the indices, then the contents of this array are
    ///   unspecified.
    ///
    /// \return The number of column indices in the row on the calling
    ///   process.  If \c lclRow does not live on the calling process,
    ///   then return zero (as there are no column indices in this row
    ///   on the calling process).
    virtual local_ordinal_type
    getLocalRowCopy (const local_ordinal_type lclRow,
                     const Teuchos::ArrayView<local_ordinal_type>& lclColInds) = 0;

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
    getLocalRowView (const local_ordinal_type lclRow,
                     Teuchos::ArrayView<const local_ordinal_type>& lclColInds) const;

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
    getGlobalRowView (const global_ordinal_type gblRow,
                      Teuchos::ArrayView<const global_ordinal_type>& gblColInds) const;

    //@}
    //! \name Implementation of Packable interface
    //@{

    /// \brief Pack this object's data for Import or Export.
    ///
    /// This is a fall-back implementation that relies only on the
    /// RowGraph interface.  If you want your subclass of RowGraph to
    /// thread-parallelize or otherwise optimize this operation for
    /// your data structure, you must override this method with your
    /// own implementation.
    virtual void
    pack (const Teuchos::ArrayView<const local_ordinal_type>& exportLIDs,
          Teuchos::Array<global_ordinal_type>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets,
          ::Tpetra::Distributor& distor) const;
    //@}
  }; // class RowGraph

} // namespace TpetraNew

#endif // TPETRANEW_ROWGRAPH_HPP
