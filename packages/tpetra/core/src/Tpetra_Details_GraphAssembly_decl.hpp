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
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {

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
/// finite-element operator.  The local owned+shared graph is built entirely on device.
/// The owned+shared -> owned communication is performed explicitly with a fused
/// Export + fillComplete.
///
/// There are two codepaths, selected by whether an owned+shared column map is
/// provided to the constructor:
///   - If a column map is provided, the global node IDs are mapped to local
///     column indices directly inside the Fill functor, using that column map.
///   - If no column map is provided, a globally-indexed local graph is built
///     first, then Tpetra::Details::makeColMap constructs the column map, and a
///     final parallel_for kernel remaps the global indices to local indices.
///
/// Typical usage:
/// \code
/// using Assembly = Tpetra::Details::GraphAssembly<>;
/// Assembly assembly(rowMap, ownedPlusSharedMap, ownedElementToNode);
/// assembly.build();
/// RCP<CrsGraph> graph = assembly.getGraph();                   // owned CrsGraph
/// RCP<CrsGraph> asmGraph = assembly.getOwnedPlusSharedGraph(); // for matrix fill
/// RCP<Export> exporter = assembly.getExporter();               // may be null (serial)
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
  /// \param ownedPlusSharedColMap [in] Optional owned+shared column map.  If
  ///   provided, it is used to map global column IDs to local column indices
  ///   inside the Fill functor.  If null, a column map is built with
  ///   Tpetra::Details::makeColMap after a globally-indexed local graph has been
  ///   assembled.
  GraphAssembly(const Teuchos::RCP<const map_type>& rowMap,
                const Teuchos::RCP<const map_type>& ownedPlusSharedMap,
                const element_to_node_type& ownedElementToNode,
                const Teuchos::RCP<const map_type>& ownedPlusSharedColMap = Teuchos::null);

  //! Destructor.
  ~GraphAssembly() = default;

  //@}
  //! @name Assembly
  //@{

  /// \brief Build the graph.
  ///
  /// This performs the local assembly, export and fillComplete to produce the
  /// owned graph (as well as the owned+shared graph and the exporter).
  void build();

  //@}
  //! @name Accessors
  //@{

  /// \brief Get the owned graph. Returns null if build() has not been called yet.
  Teuchos::RCP<crs_graph_type> getGraph() const;

  /// \brief Get the owned+shared graph.
  ///
  /// It will be packed, locally indexed and fillComplete. row map == column map == owned+shared map,
  /// domain map == owned map. Used for filling matrix values and RHS on the
  /// owned+shared map (using owned elements) before Exporting to the owned matrix/vector.
  /// Returns null if build() has not been called yet.
  Teuchos::RCP<crs_graph_type> getOwnedPlusSharedGraph() const;

  /// \brief The Export object for (owned+shared -> owned) communication.
  ///
  /// Null in the serial / single-rank case, where no communication is needed.
  /// Useful for exporting the owned+shared matrix to owned.
  Teuchos::RCP<const export_type> getExporter() const;

  /// \brief Get the owned+shared column map used by the assembled graph.
  ///
  /// If a column map was provided to the constructor, this returns it.
  /// Otherwise, it returns the column map built by makeColMap during build().
  /// Returns null if build() has not been called yet.
  Teuchos::RCP<const map_type> getOwnedPlusSharedColMap() const;

  //@}

 private:
  //! The owned (one-to-one) row/domain/range map.
  Teuchos::RCP<const map_type> rowMap_;
  //! The owned+shared (overlapping) map.
  Teuchos::RCP<const map_type> ownedPlusSharedMap_;
  //! Element-to-node connectivity (global node IDs) of the owned elements.
  element_to_node_type ownedElementToNode_;
  //! Optional user-provided owned+shared column map (may be null).
  Teuchos::RCP<const map_type> ownedPlusSharedColMap_;

  //! The assembled owned graph (result of build()).
  Teuchos::RCP<crs_graph_type> graph_;
  //! The intermediate owned+shared assembly graph (result of build()).
  Teuchos::RCP<crs_graph_type> ownedPlusSharedGraph_;
  //! The Export used for the owned+shared -> owned communication (may be null).
  Teuchos::RCP<const export_type> exporter_;
};

}  // namespace Details
}  // namespace Tpetra

#endif  // TPETRA_DETAILS_GRAPHASSEMBLY_DECL_HPP
