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
#include <initializer_list>

namespace Tpetra {

/// \class Geometry
/// \brief A device-compatible collection of one or more finite-element mesh
///   blocks, each described by an element-to-node connectivity view.
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
/// Geometry bundles a variadic number of these element-to-node views into one
/// object.  Each block may have a different number of nodes per element (the
/// second extent of its view).  The blocks are assembled together into a single
/// CrsGraph.
///
/// Geometry is trivially copyable and captures its blocks by (reference-counted
/// Kokkos::View) value, so it can be captured by value in a KOKKOS_LAMBDA and
/// used inside device kernels.  A fixed compile-time capacity
/// (\c maxNumBlocks) is used so that the block views can be stored in a
/// Kokkos::Array (which is device-friendly), avoiding any host-only container.
///
/// Typical usage:
/// \code
/// using geom_type = Tpetra::Geometry<GO, Node>;
/// // Two blocks: triangles and quads.
/// geom_type geom(triElementToNode, quadElementToNode);
/// // or, equivalently, add blocks one at a time:
/// geom_type geom2;
/// geom2.addBlock(triElementToNode);
/// geom2.addBlock(quadElementToNode);
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

  /// \brief Maximum number of mesh blocks a Geometry can hold.
  ///
  /// A compile-time capacity so that the blocks can live in a device-friendly
  /// Kokkos::Array.  Increase this if more element types are needed.
  static constexpr int maxNumBlocks = 8;

  //! Default constructor: an empty Geometry with no blocks.
  KOKKOS_DEFAULTED_FUNCTION
  Geometry() = default;

  /// \brief Variadic constructor from one or more element-to-node views.
  ///
  /// Each argument must be convertible to #element_to_node_type.  There must be
  /// at least one, and no more than #maxNumBlocks arguments.
  template <class... Views,
            class = typename std::enable_if<(sizeof...(Views) >= 1)>::type>
  Geometry(const Views&... views) {
    static_assert(sizeof...(Views) <= static_cast<size_t>(maxNumBlocks),
                  "Tpetra::Geometry: too many mesh blocks (exceeds maxNumBlocks).");
    addBlocks(views...);
  }

  //! Add a mesh block (element-to-node connectivity view) to the Geometry.
  void addBlock(const element_to_node_type& elementToNode) {
    if (numBlocks_ < maxNumBlocks) {
      blocks_[numBlocks_] = elementToNode;
      ++numBlocks_;
    }
  }

  //! The number of mesh blocks currently held.
  KOKKOS_INLINE_FUNCTION
  int getNumBlocks() const { return numBlocks_; }

  //! The element-to-node connectivity view of block \c b.
  KOKKOS_INLINE_FUNCTION
  const element_to_node_type& getBlock(const int b) const { return blocks_[b]; }

  //! The number of owned elements in block \c b.
  KOKKOS_INLINE_FUNCTION
  size_t getNumElements(const int b) const { return blocks_[b].extent(0); }

  //! The number of nodes per element of block \c b.
  KOKKOS_INLINE_FUNCTION
  int getNodesPerElement(const int b) const {
    return static_cast<int>(blocks_[b].extent(1));
  }

  //! The total number of owned elements across all blocks.
  KOKKOS_INLINE_FUNCTION
  size_t getTotalNumElements() const {
    size_t total = 0;
    for (int b = 0; b < numBlocks_; ++b) total += blocks_[b].extent(0);
    return total;
  }

  //! The maximum number of nodes per element across all blocks.
  KOKKOS_INLINE_FUNCTION
  int getMaxNodesPerElement() const {
    int m = 0;
    for (int b = 0; b < numBlocks_; ++b) {
      const int npe = static_cast<int>(blocks_[b].extent(1));
      if (npe > m) m = npe;
    }
    return m;
  }

 private:
  //! Recursively add a parameter pack of views.
  void addBlocks() {}
  template <class First, class... Rest>
  void addBlocks(const First& first, const Rest&... rest) {
    addBlock(element_to_node_type(first));
    addBlocks(rest...);
  }

  //! Fixed-capacity, device-friendly storage of the block views.
  Kokkos::Array<element_to_node_type, maxNumBlocks> blocks_;
  //! The number of blocks actually stored.
  int numBlocks_ = 0;
};

}  // namespace Tpetra

#endif  // TPETRA_GEOMETRY_HPP
