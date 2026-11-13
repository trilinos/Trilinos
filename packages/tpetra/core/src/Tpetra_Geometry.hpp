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
#include "Tpetra_Map_fwd.hpp"
#include "Kokkos_Core.hpp"
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace Tpetra {

/// \class Geometry
/// \brief The element connectivity of a finite-element mesh: one or more
///   mesh blocks, each described by element-to-node views.
///
/// \tparam LocalOrdinal  Type of local indices; matches CrsGraph.
/// \tparam GlobalOrdinal Type of global indices; matches CrsGraph.
/// \tparam Node          Kokkos Node type; matches CrsGraph.
///
/// Geometry is the input to Tpetra::assembleFECrsGraph and
/// Tpetra::assembleFECrsGraphWithMaps, which build the sparse graph of a
/// finite-element operator from it.
///
/// <b>Blocks.</b>  A block is a rank-2 Kokkos::View indexed as
/// (localElementIndex, nodeOfElement) that lists the nodes (or, more generally,
/// degrees of freedom) of each locally-owned element.  All elements of one view
/// have the same number of nodes, so a mesh with several element types (e.g.
/// triangles and quadrilaterals) or several element blocks is described by
/// adding one block per type.  There is no limit on the number of blocks.
///
/// <b>Square and rectangular blocks.</b>  A block added with one view is
/// "square": each element contributes a dense clique, coupling every node of the
/// element (row) with every node of the element (column).  A block added with
/// two views, \c elementToRow and \c elementToCol, is "rectangular": each
/// element couples every node in its \c elementToRow row with every node in its
/// \c elementToCol row.  Both views must have the same number of elements
/// (extent 0), but may have different numbers of nodes per element (extent 1).
/// Use this when the row and column degrees of freedom differ (e.g. an
/// off-diagonal block of a multiphysics operator), or when an element's columns
/// extend beyond its own nodes (e.g. coupling to a neighboring element's
/// nodes).
///
/// <b>Padding.</b>  Entries equal to
/// Teuchos::OrdinalTraits<GlobalOrdinal>::invalid() (for global blocks) or
/// Teuchos::OrdinalTraits<LocalOrdinal>::invalid() (for local blocks) are
/// ignored.  This allows elements with fewer nodes than the view's extent(1).
///
/// <b>Global or local indices.</b>  Blocks hold either global indices (added
/// with addGlobalBlock) or local indices (added with addLocalBlock).  All blocks
/// of one Geometry must use the same kind of index.  Local indices refer to the
/// maps passed to the assembly function:
///   - row indices (the single view of a square block, or \c elementToRow of a
///     rectangular block) are local indices of the owned+shared row map;
///   - column indices of a rectangular block (\c elementToCol) are local indices
///     of the owned+shared column map if one is given, else of the owned+shared
///     domain map if one is given, else of the owned+shared row map.
///
/// The views are held by reference-counted shallow copy, and must remain
/// unmodified until assembly is done.  They are stored with
/// Kokkos::LayoutStride, so a view of any layout (or a strided subview) can be
/// passed in without a copy.
///
/// Typical usage:
/// \code
/// Tpetra::Geometry<LO, GO, Node> geom;
/// geom.addGlobalBlock(triElementToNode);   // triangles
/// geom.addGlobalBlock(quadElementToNode);  // quadrilaterals
/// auto graph = Tpetra::assembleFECrsGraph(geom, ownedPlusSharedMap);
/// \endcode
template <class LocalOrdinal  = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node          = ::Tpetra::Details::DefaultTypes::node_type>
class Geometry {
 public:
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;
  using device_type         = typename Node::device_type;
  using map_type            = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;

  //! Type of a block view holding global indices.
  using global_block_type =
      Kokkos::View<const global_ordinal_type**, Kokkos::LayoutStride, device_type>;
  //! Type of a block view holding local indices.
  using local_block_type =
      Kokkos::View<const local_ordinal_type**, Kokkos::LayoutStride, device_type>;

  //! Default constructor: an empty Geometry with no blocks.
  Geometry() = default;

  /// \brief Construct from one or more square blocks of global indices.
  ///
  /// Equivalent to calling addGlobalBlock(view) for each argument, in order.
  template <class... Views,
            class = std::enable_if_t<
                (sizeof...(Views) >= 1) &&
                (std::is_convertible_v<const Views&, global_block_type> && ...)>>
  explicit Geometry(const Views&... views) {
    blocks_.reserve(sizeof...(Views));
    (addGlobalBlock(global_block_type(views)), ...);
  }

  //! @name Adding blocks
  //@{

  //! Add a square block of global indices (each element is a clique).
  void addGlobalBlock(const global_block_type& elementToNode) {
    setIndexKind(IndexKind::Global);
    Block block;
    block.globalRows = elementToNode;
    blocks_.push_back(block);
  }

  //! Add a rectangular block of global indices.
  void addGlobalBlock(const global_block_type& elementToRow,
                      const global_block_type& elementToCol) {
    setIndexKind(IndexKind::Global);
    checkSameNumElements(elementToRow.extent(0), elementToCol.extent(0));
    Block block;
    block.globalRows = elementToRow;
    block.globalCols = elementToCol;
    block.hasCols    = true;
    blocks_.push_back(block);
  }

  //! Add a square block of local indices (each element is a clique).
  void addLocalBlock(const local_block_type& elementToNode) {
    setIndexKind(IndexKind::Local);
    Block block;
    block.localRows = elementToNode;
    blocks_.push_back(block);
  }

  //! Add a rectangular block of local indices.
  void addLocalBlock(const local_block_type& elementToRow,
                     const local_block_type& elementToCol) {
    setIndexKind(IndexKind::Local);
    checkSameNumElements(elementToRow.extent(0), elementToCol.extent(0));
    Block block;
    block.localRows = elementToRow;
    block.localCols = elementToCol;
    block.hasCols   = true;
    blocks_.push_back(block);
  }

  //@}
  //! @name Queries
  //@{

  //! Whether the blocks hold local indices (false if global, or if there are no blocks).
  bool hasLocalIndices() const { return kind_ == IndexKind::Local; }

  //! Whether the blocks hold global indices (true if there are no blocks).
  bool hasGlobalIndices() const { return kind_ != IndexKind::Local; }

  //! The number of mesh blocks.
  int getNumBlocks() const { return static_cast<int>(blocks_.size()); }

  //! Whether block \c b is rectangular (has separate column connectivity).
  bool isRectangular(const int b) const { return blocks_[b].hasCols; }

  //! The number of owned elements in block \c b.
  size_t getNumElements(const int b) const {
    return hasLocalIndices() ? blocks_[b].localRows.extent(0) : blocks_[b].globalRows.extent(0);
  }

  //! The total number of owned elements across all blocks.
  size_t getTotalNumElements() const {
    size_t total = 0;
    for (int b = 0; b < getNumBlocks(); ++b) total += getNumElements(b);
    return total;
  }

  //! The row view of block \c b (requires global indices).
  const global_block_type& getGlobalRows(const int b) const { return blocks_[b].globalRows; }
  //! The column view of block \c b; the row view if the block is square (requires global indices).
  const global_block_type& getGlobalCols(const int b) const {
    return blocks_[b].hasCols ? blocks_[b].globalCols : blocks_[b].globalRows;
  }
  //! The row view of block \c b (requires local indices).
  const local_block_type& getLocalRows(const int b) const { return blocks_[b].localRows; }
  //! The column view of block \c b; the row view if the block is square (requires local indices).
  const local_block_type& getLocalCols(const int b) const {
    return blocks_[b].hasCols ? blocks_[b].localCols : blocks_[b].localRows;
  }

  //@}

 private:
  enum class IndexKind { None,
                         Global,
                         Local };

  struct Block {
    global_block_type globalRows, globalCols;
    local_block_type localRows, localCols;
    bool hasCols = false;
  };

  void setIndexKind(const IndexKind kind) {
    if (kind_ != IndexKind::None && kind_ != kind)
      throw std::invalid_argument(
          "Tpetra::Geometry: all blocks must use the same kind of index (all "
          "addGlobalBlock, or all addLocalBlock).");
    kind_ = kind;
  }

  static void checkSameNumElements(const size_t numRowElements, const size_t numColElements) {
    if (numRowElements != numColElements)
      throw std::invalid_argument(
          "Tpetra::Geometry: the row and column views of a rectangular block must "
          "have the same number of elements (extent(0)).");
  }

  IndexKind kind_ = IndexKind::None;
  std::vector<Block> blocks_;
};

}  // namespace Tpetra

#endif  // TPETRA_GEOMETRY_HPP
