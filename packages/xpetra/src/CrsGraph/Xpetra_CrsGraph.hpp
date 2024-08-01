// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_CRSGRAPH_HPP
#define XPETRA_CRSGRAPH_HPP

#include <Teuchos_ParameterList.hpp>

#include <Teuchos_Describable.hpp>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DistObject.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_Map.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include <Kokkos_StaticCrsGraph.hpp>
#endif

namespace Xpetra {

using Teuchos::ParameterList;

struct RowInfo {
  size_t localRow;
  size_t allocSize;
  size_t numEntries;
  size_t offset1D;
};

enum ELocalGlobal {
  LocalIndices,
  GlobalIndices
};

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class CrsGraph
  : /*public RowGraph<>,*/ public DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

  //! @name Constructor/Destructor Methods
  //@

  //! Destructor.
  virtual ~CrsGraph() {}

  //@}

  //! @name Insertion/Removal Methods
  //@{

  //! Insert global indices into the graph.
  virtual void insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &indices) = 0;

  //! Insert local indices into the graph.
  virtual void insertLocalIndices(const LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &indices) = 0;

  //! Remove all graph indices from the specified local row.
  virtual void removeLocalIndices(LocalOrdinal localRow) = 0;

  //! Allocates the 1D pointer arrays of the graph
  virtual void allocateAllIndices(size_t numNonZeros, ArrayRCP<size_t> &rowptr, ArrayRCP<LocalOrdinal> &colind) = 0;

  //! Sets the 1D pointer arrays of the graph.
  virtual void setAllIndices(const ArrayRCP<size_t> &rowptr, const ArrayRCP<LocalOrdinal> &colind) = 0;

  //! Gets the 1D pointer arrays of the graph.
  virtual void getAllIndices(ArrayRCP<const size_t> &rowptr, ArrayRCP<const LocalOrdinal> &colind) const = 0;

  //@}

  //! @name Transformational Methods
  //@{

  //! Signal that data entry is complete, specifying domain and range maps.
  virtual void fillComplete(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rangeMap, const RCP<ParameterList> &params = null) = 0;

  //! Signal that data entry is complete.
  virtual void fillComplete(const RCP<ParameterList> &params = null) = 0;

  //! Expert version of fillComplete
  virtual void
  expertStaticFillComplete(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap,
                           const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rangeMap,
                           const RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > &importer = null,
                           const RCP<const Export<LocalOrdinal, GlobalOrdinal, Node> > &exporter = null,
                           const RCP<Teuchos::ParameterList> &params                             = null) = 0;

  //@}

  //! @name Methods implementing RowGraph.
  //@{

  //! Returns the communicator.
  virtual RCP<const Comm<int> > getComm() const = 0;

  //! Returns the Map that describes the row distribution in this graph.
  virtual RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getRowMap() const = 0;

  //! Returns the Map that describes the column distribution in this graph.
  virtual RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getColMap() const = 0;

  //! Returns the Map associated with the domain of this graph.
  virtual RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const = 0;

  //! Returns the Map associated with the domain of this graph.
  virtual RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const = 0;

  //! Returns the importer associated with this graph.
  virtual RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > getImporter() const = 0;

  //! Returns the exporter associated with this graph.
  virtual RCP<const Export<LocalOrdinal, GlobalOrdinal, Node> > getExporter() const = 0;

  //! Returns the number of global rows in the graph.
  virtual global_size_t getGlobalNumRows() const = 0;

  //! Returns the number of global columns in the graph.
  virtual global_size_t getGlobalNumCols() const = 0;

  //! Returns the number of graph rows owned on the calling node.
  virtual size_t getLocalNumRows() const = 0;

  //! Returns the number of columns connected to the locally owned rows of this graph.
  virtual size_t getLocalNumCols() const = 0;

  //! Returns the index base for global indices for this graph.
  virtual GlobalOrdinal getIndexBase() const = 0;

  //! Returns the global number of entries in the graph.
  virtual global_size_t getGlobalNumEntries() const = 0;

  //! Returns the local number of entries in the graph.
  virtual size_t getLocalNumEntries() const = 0;

  //! Returns the current number of entries on this node in the specified global row.
  virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

  //! Returns the current number of entries on this node in the specified local row.
  virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

  //! Returns the current number of allocated entries for this node in the specified global row .
  virtual size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

  //! Returns the current number of allocated entries on this node in the specified local row.
  virtual size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const = 0;

  //! Maximum number of entries in all rows over all processes.
  virtual size_t getGlobalMaxNumRowEntries() const = 0;

  //! Maximum number of entries in all rows owned by the calling process.
  virtual size_t getLocalMaxNumRowEntries() const = 0;

  //! Whether the graph has a column Map.
  virtual bool hasColMap() const = 0;

  //! Whether column indices are stored using local indices on the calling process.
  virtual bool isLocallyIndexed() const = 0;

  //! Whether column indices are stored using global indices on the calling process.
  virtual bool isGloballyIndexed() const = 0;

  //! Whether fillComplete() has been called and the graph is in compute mode.
  virtual bool isFillComplete() const = 0;

  //! Returns true if storage has been optimized.
  virtual bool isStorageOptimized() const = 0;

  //! Return a const, nonpersisting view of global indices in the given row.
  virtual void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &Indices) const = 0;

  //! Return a const, nonpersisting view of local indices in the given row.
  virtual void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices) const = 0;

  //! Force the computation of global constants if we don't have them
  virtual void computeGlobalConstants() = 0;

  //@}

  //! @name Tpetra-specific routines
  //@{
#ifdef HAVE_XPETRA_TPETRA
  typedef typename node_type::execution_space execution_space;
  typedef typename node_type::device_type device_type;
  typedef Kokkos::StaticCrsGraph<LocalOrdinal, Kokkos::LayoutLeft, device_type, void, size_t> local_graph_type;

  /// \brief Get the local graph.
  ///
  /// \warning THIS IS AN EXPERT MODE FUNCTION.  THIS IS AN
  ///   IMPLEMENTATION DETAIL.  DO NOT CALL THIS FUNCTION!!!
  ///
  /// This is only a valid representation of the local graph if the
  /// (global) graph is fill complete.
  virtual typename local_graph_type::HostMirror getLocalGraphHost() const = 0;
  virtual local_graph_type getLocalGraphDevice() const                    = 0;

  //! Get offsets of the diagonal entries in the matrix.
  virtual void getLocalDiagOffsets(const Kokkos::View<size_t *, device_type, Kokkos::MemoryUnmanaged> &offsets) const = 0;

#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  virtual std::string description() const = 0;

  //! Print the object with some verbosity level to an FancyOStream object.
  virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const = 0;

  //@}

};  // CrsGraph class

}  // namespace Xpetra

#define XPETRA_CRSGRAPH_SHORT
#endif  // XPETRA_CRSGRAPH_HPP
