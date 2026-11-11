// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_GEOMETRY_HPP
#define TPETRA_GEOMETRY_HPP

/// \file Tpetra_Geometry.hpp
/// \brief Declaration and definition of the Tpetra::Geometry class.

#include "Tpetra_Details_DefaultTypes.hpp"
#include "Kokkos_Core.hpp"
#include <type_traits>
#include <utility>
#include <vector>

namespace Tpetra {

/// \class Geometry
/// \brief A collection of one or more finite-element mesh blocks, each
///   described by an element-to-node connectivity view.
///
/// \tparam GlobalOrdinal Type of global node indices; matches CrsGraph.
/// \tparam Node          Kokkos Node type; matches CrsGraph.
///
/// A single \c ownedElementToNode view (a rank-2 Kokkos::View indexed as
/// (localElementIndex, nodeOfElement) of GLOBAL node IDs) can only describe a
/// single element type, because all elements in one non-ragged 2D view have the
/// same number of adjacent nodes.  A finite-element mesh may however contain
/// several element types (mesh "blocks") -- for instance some triangles and
/// some quadrilaterals -- that are assembled into a single graph.
///
/// Geometry bundles any number of these element-to-node views into one object.
/// Each block may have a different number of nodes per element (the second
/// extent of its view).  The blocks are assembled together into a single
/// CrsGraph.
///
/// Geometry is a host-side object: the block views are held in a
/// std::vector, so there is no limit on the number of blocks.  The views
/// themselves live in device memory.  Graph assembly (see
/// Tpetra::Details::GraphAssembly) combines all blocks into a single
/// device-resident element-to-node sparse structure (rows = elements,
/// columns = nodes) before running its device kernels.
///
/// Typical usage:
/// \code
/// using geom_type = Tpetra::Geometry<GO, Node>;
/// // Two blocks: triangles and quads.
/// geom_type geom(triElementToNode, quadElementToNode);
/// // or, from a runtime-sized list of blocks:
/// std::vector<geom_type::element_to_node_type> blocks = ...;
/// geom_type geom2(blocks);
/// // or, equivalently, add blocks one at a time:
/// geom_type geom3;
/// geom3.addBlock(triElementToNode);
/// geom3.addBlock(quadElementToNode);
/// \endcode
template <class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node          = ::Tpetra::Details::DefaultTypes::node_type>
class Geometry {
 public:
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;
  using device_type         = typename Node::device_type;

  /// \brief Type of a single block's element-to-node connectivity view.
  ///
  /// A rank-2 Kokkos::View, indexed as (localElementIndex, nodeOfElement),
  /// holding the GLOBAL node IDs of the nodes adjacent to each locally-owned
  /// element of the block.  The first extent is the number of owned elements in
  /// the block; the second extent is the number of nodes per element of the
  /// block.
  using element_to_node_type =
      Kokkos::View<const global_ordinal_type**, device_type>;

  //! Default constructor: an empty Geometry with no blocks.
  Geometry() = default;

  /// \brief Variadic constructor from one or more element-to-node views.
  ///
  /// Each argument must be convertible to #element_to_node_type.
  template <class... Views,
            class = std::enable_if_t<
                (sizeof...(Views) >= 1) &&
                (std::is_convertible_v<const Views&, element_to_node_type> && ...)>>
  Geometry(const Views&... views) {
    blocks_.reserve(sizeof...(Views));
    (addBlock(element_to_node_type(views)), ...);
  }

  //! Constructor from a runtime-sized list of element-to-node views.
  explicit Geometry(const std::vector<element_to_node_type>& blocks)
    : blocks_(blocks) {}

  //! Constructor from a runtime-sized list of element-to-node views.
  explicit Geometry(std::vector<element_to_node_type>&& blocks)
    : blocks_(std::move(blocks)) {}

  //! Add a mesh block (element-to-node connectivity view) to the Geometry.
  void addBlock(const element_to_node_type& elementToNode) {
    blocks_.push_back(elementToNode);
  }

  //! The number of mesh blocks currently held.
  int getNumBlocks() const { return static_cast<int>(blocks_.size()); }

  //! All mesh blocks' element-to-node connectivity views.
  const std::vector<element_to_node_type>& getBlocks() const { return blocks_; }

  //! The element-to-node connectivity view of block \c b.
  const element_to_node_type& getBlock(const int b) const { return blocks_[b]; }

  //! The number of owned elements in block \c b.
  size_t getNumElements(const int b) const { return blocks_[b].extent(0); }

  //! The number of nodes per element of block \c b.
  int getNodesPerElement(const int b) const {
    return static_cast<int>(blocks_[b].extent(1));
  }

  //! The total number of owned elements across all blocks.
  size_t getTotalNumElements() const {
    size_t total = 0;
    for (const auto& block : blocks_) total += block.extent(0);
    return total;
  }

  //! The maximum number of nodes per element across all blocks.
  int getMaxNodesPerElement() const {
    int m = 0;
    for (const auto& block : blocks_) {
      const int npe = static_cast<int>(block.extent(1));
      if (npe > m) m = npe;
    }
    return m;
  }

 private:
  //! The block views (host-side container of device views).
  std::vector<element_to_node_type> blocks_;
};

}  // namespace Tpetra

#endif  // TPETRA_GEOMETRY_HPP
