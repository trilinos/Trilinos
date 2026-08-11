// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_GRAPHASSEMBLY_DECL_HPP
#define TPETRA_GRAPHASSEMBLY_DECL_HPP

/// \file Tpetra_GraphAssembly_decl.hpp
///
/// Declaration of Tpetra::Experimental::GraphAssembly.

#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Experimental {

/// \class GraphAssembly
/// \brief Assemble a Tpetra::CrsGraph for a finite-element problem entirely
///   with Kokkos kernels, from a locally-owned element-to-node connectivity.
///
/// \tparam LocalOrdinal  Type of local indices; matches CrsGraph.
/// \tparam GlobalOrdinal Type of global indices; matches CrsGraph.
/// \tparam Node          Kokkos Node type; matches CrsGraph.
///
/// Given the distribution of the mesh nodes (an "owned" row map and an
/// "owned+shared" map) and the connectivity of the locally-owned elements
/// (element -> global node IDs), this class builds the sparse graph of the
/// finite-element operator.  The graph structure is built directly on device
/// (count / prefix-sum / fill kernels using a per-team hash table), wrapped in
/// a packed, locally-indexed CrsGraph over the owned+shared map, and then the
/// owned+shared -> owned communication is performed explicitly with a fused
/// Export + fillComplete.  No host serial insertGlobalIndices loop is used.
///
/// Typical usage:
/// \code
/// using GA = Tpetra::Experimental::GraphAssembly<>;
/// GA assembler(rowMap, ownedPlusSharedMap, ownedElementToNode);
/// assembler.build();
/// auto graph = assembler.getGraph();               // owned CrsGraph
/// auto asmGraph = assembler.getOwnedPlusSharedGraph(); // for matrix/RHS fill
/// auto exporter = assembler.getExporter();          // may be null (serial)
/// \endcode
///
/// The owned+shared assembly graph and the Export object are also exposed so
/// that a caller can assemble the matrix values / right-hand side into
/// owned+shared objects sharing this packed structure and then Export them to
/// the owned objects with the same communication pattern.
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
  using export_type    = ::Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>;

  using device_type     = typename Node::device_type;
  using execution_space = typename device_type::execution_space;

  /// \brief Type of the element-to-node connectivity view.
  ///
  /// A rank-2 Kokkos::View, indexed as (localElementIndex, nodeOfElement),
  /// holding the GLOBAL node IDs of the nodes adjacent to each locally-owned
  /// element.  The first extent is the number of owned elements; the second
  /// extent is the number of nodes per element.  It must be accessible from
  /// the graph's execution space.
  using element_to_node_type =
      Kokkos::View<const global_ordinal_type**, device_type>;
  //@}

  //! @name Constructors/destructor
  //@{

  /// \brief Constructor.
  ///
  /// \param rowMap [in] The one-to-one ("owned") map of the mesh nodes.  This
  ///   becomes the row, domain and range map of the assembled graph.
  /// \param ownedPlusSharedMap [in] The overlapping ("owned+shared") map: every
  ///   node touched by a locally-owned element must be in this map.  Used as
  ///   both the row and column map of the intermediate assembly graph.
  /// \param ownedElementToNode [in] The element-to-node connectivity of the
  ///   locally-owned elements (global node IDs).  See #element_to_node_type.
  GraphAssembly(const Teuchos::RCP<const map_type>& rowMap,
                const Teuchos::RCP<const map_type>& ownedPlusSharedMap,
                const element_to_node_type& ownedElementToNode);

  //! Destructor.
  ~GraphAssembly() = default;

  //@}
  //! @name Assembly
  //@{

  /// \brief Build the graph.
  ///
  /// After this returns, getGraph() gives the assembled, fill-complete owned
  /// CrsGraph.  It is safe to call more than once (the graph is rebuilt).
  void build();

  //@}
  //! @name Accessors
  //@{

  /// \brief The assembled, fill-complete graph over the owned (row) map.
  ///
  /// Returns null if build() has not been called yet.
  Teuchos::RCP<crs_graph_type> getGraph() const;

  /// \brief The intermediate assembly graph over the owned+shared map.
  ///
  /// Packed and locally indexed; row map == column map == owned+shared map,
  /// domain map == owned map.  Useful for assembling matrix values / RHS on the
  /// owned+shared map before Exporting to the owned objects.  Returns null if
  /// build() has not been called yet.
  Teuchos::RCP<crs_graph_type> getOwnedPlusSharedGraph() const;

  /// \brief The Export object (owned+shared -> owned) used to communicate.
  ///
  /// Null in the serial / single-rank case, where no communication is needed.
  Teuchos::RCP<const export_type> getExporter() const;

  //@}

 private:
  //! The owned (one-to-one) row/domain/range map.
  Teuchos::RCP<const map_type> rowMap_;
  //! The owned+shared (overlapping) map.
  Teuchos::RCP<const map_type> ownedPlusSharedMap_;
  //! Element-to-node connectivity (global node IDs) of the owned elements.
  element_to_node_type ownedElementToNode_;

  //! The assembled owned graph (result of build()).
  Teuchos::RCP<crs_graph_type> graph_;
  //! The intermediate owned+shared assembly graph.
  Teuchos::RCP<crs_graph_type> ownedPlusSharedGraph_;
  //! The Export used for the owned+shared -> owned communication (may be null).
  Teuchos::RCP<const export_type> exporter_;
};

}  // namespace Experimental
}  // namespace Tpetra

#endif  // TPETRA_GRAPHASSEMBLY_DECL_HPP
