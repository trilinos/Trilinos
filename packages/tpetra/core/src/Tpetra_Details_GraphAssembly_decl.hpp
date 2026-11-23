// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_GRAPHASSEMBLY_DECL_HPP
#define TPETRA_DETAILS_GRAPHASSEMBLY_DECL_HPP

/// \file Tpetra_Details_GraphAssembly_decl.hpp
///
/// Declaration of Tpetra::Details::GraphAssembly.

#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Geometry.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {

/// \class GraphAssembly
/// \brief Assemble the graph of a finite-element operator with Kokkos kernels,
///   from the element connectivity of the locally-owned elements.
///
/// \tparam LocalOrdinal  Type of local indices; matches CrsGraph.
/// \tparam GlobalOrdinal Type of global indices; matches CrsGraph.
/// \tparam Node          Kokkos Node type; matches CrsGraph.
///
/// This is the implementation of Tpetra::assembleFECrsGraph and
/// Tpetra::assembleFECrsGraphWithMaps; most users should call those instead.
///
/// The maps have the same meaning as the corresponding arguments of the
/// FECrsGraph constructors.  The result is <i>identical</i> to what FECrsGraph
/// produces when each element's (row, column) couplings are inserted with
/// insertGlobalIndices between beginAssembly() and endAssembly(): the same
/// column map (same global indices in the same order), and the same local row
/// offsets and local column indices, in the same order, for both the
/// owned+shared graph and the owned graph.
///
/// Algorithm:
///   1. The Geometry's blocks are combined into one element-to-row and one
///      element-to-column structure, translating every index to a local index
///      exactly once.  Rows are local indices of the owned+shared row map;
///      columns are local indices of the "assembly column map": the
///      owned+shared column map if given, else the owned+shared domain map if
///      given, else the owned+shared row map.  Every index must be found in its
///      map (an exception is thrown otherwise).
///   2. The local owned+shared graph is assembled on device (one team per row,
///      using a scratch hash table to deduplicate columns).
///   3. In parallel, a fused Export + fillComplete merges each owned row with
///      the contributions of other processes.
///   4. The final column map is built exactly as FECrsGraph builds it (if no
///      column map was given), and all column indices are translated to it and
///      sorted.
///
/// Typical usage:
/// \code
/// Tpetra::Details::GraphAssembly<LO, GO, NT> assembly(geometry, ownedRowMap, ownedPlusSharedRowMap);
/// assembly.build();
/// RCP<CrsGraph> ownedPlusSharedGraph = assembly.getOwnedPlusSharedGraph();
/// RCP<CrsGraph> ownedGraph           = assembly.getGraph();
/// \endcode
template <class LocalOrdinal  = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node          = ::Tpetra::Details::DefaultTypes::node_type>
class GraphAssembly {
 public:
  //! @name Typedefs
  //@{
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;

  using map_type       = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using crs_graph_type = ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using import_type    = ::Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
  using export_type    = ::Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>;
  using geometry_type  = ::Tpetra::Geometry<LocalOrdinal, GlobalOrdinal, Node>;

  using device_type     = typename Node::device_type;
  using execution_space = typename device_type::execution_space;
  //@}

  //! @name Constructor/destructor
  //@{

  /// \brief Constructor.
  ///
  /// \param geometry [in] The element connectivity of the locally-owned
  ///   elements.  See Tpetra::Geometry.
  /// \param ownedRowMap [in] The one-to-one ("owned") row map.
  /// \param ownedPlusSharedRowMap [in] The owned rows, followed by the shared
  ///   rows (rows that locally-owned elements contribute to but this process
  ///   does not own).  Must be locally fitted to ownedRowMap.
  /// \param ownedPlusSharedColMap [in] Optional column map.  If null, a column
  ///   map is built.
  /// \param ownedPlusSharedDomainMap [in] Optional domain map of the
  ///   owned+shared graph.  If null, ownedDomainMap is used.
  /// \param ownedDomainMap [in] Optional domain map of the owned graph.  If
  ///   null, ownedRowMap is used.
  /// \param ownedRangeMap [in] Optional range map of the owned graph.  If null,
  ///   ownedRowMap is used.
  /// \param ownedPlusSharedToOwnedImporter [in] Optional Import from
  ///   ownedRowMap to ownedPlusSharedRowMap.  Reused (reversed) for the
  ///   owned+shared -> owned communication if given.
  GraphAssembly(const geometry_type& geometry,
                const Teuchos::RCP<const map_type>& ownedRowMap,
                const Teuchos::RCP<const map_type>& ownedPlusSharedRowMap,
                const Teuchos::RCP<const map_type>& ownedPlusSharedColMap       = Teuchos::null,
                const Teuchos::RCP<const map_type>& ownedPlusSharedDomainMap    = Teuchos::null,
                const Teuchos::RCP<const map_type>& ownedDomainMap              = Teuchos::null,
                const Teuchos::RCP<const map_type>& ownedRangeMap               = Teuchos::null,
                const Teuchos::RCP<const import_type>& ownedPlusSharedToOwnedImporter = Teuchos::null);

  //! Destructor.
  ~GraphAssembly() = default;

  //@}
  //! @name Assembly
  //@{

  /// \brief Build the graph.  Collective over the maps' communicator.
  void build();

  //@}
  //! @name Accessors
  //@{

  /// \brief The owned+shared graph.  Null before build().
  ///
  /// Fill complete and locally indexed.  Row map: ownedPlusSharedRowMap.  Its
  /// owned rows (the leading rows) hold the complete, merged structure; its
  /// shared rows hold only this process' contributions.
  Teuchos::RCP<crs_graph_type> getOwnedPlusSharedGraph() const;

  /// \brief The owned graph.  Null before build().
  ///
  /// Fill complete; its row map is ownedRowMap.  It is a view of the leading
  /// (owned) rows of the owned+shared graph, sharing its column map.  Built on
  /// the first call.
  Teuchos::RCP<crs_graph_type> getGraph() const;

  //@}

 private:
  geometry_type geometry_;
  Teuchos::RCP<const map_type> ownedRowMap_;
  Teuchos::RCP<const map_type> ownedPlusSharedRowMap_;
  Teuchos::RCP<const map_type> ownedPlusSharedColMap_;
  Teuchos::RCP<const map_type> ownedPlusSharedDomainMap_;
  Teuchos::RCP<const map_type> ownedDomainMap_;
  Teuchos::RCP<const map_type> ownedRangeMap_;
  Teuchos::RCP<const import_type> importer_;

  //! The assembled owned+shared graph (result of build()).
  Teuchos::RCP<crs_graph_type> ownedPlusSharedGraph_;
  //! The owned graph (built lazily by getGraph()).
  mutable Teuchos::RCP<crs_graph_type> graph_;
};

}  // namespace Details
}  // namespace Tpetra

#endif  // TPETRA_DETAILS_GRAPHASSEMBLY_DECL_HPP
