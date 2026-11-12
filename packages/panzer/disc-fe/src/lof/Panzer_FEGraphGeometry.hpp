// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_FE_GRAPH_GEOMETRY_HPP
#define PANZER_FE_GRAPH_GEOMETRY_HPP

#include "PanzerDiscFE_config.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_ConnManager.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Tpetra_Geometry.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_Core.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace panzer {

/** \brief Describe the element connectivity of a finite-element graph as a
  *        Tpetra::Geometry, for Tpetra::assembleFECrsGraphWithMaps.
  *
  * Each locally owned element (as listed by \c elementIndexer's element blocks)
  * couples every row DOF of the element (from \c rowIndexer) with every column
  * DOF of the element (from \c colIndexer), and, if \c connMgr has associated
  * neighbors, with every column DOF of the element's associated neighbors.  This
  * is exactly the coupling the linear object factories insert into their graphs.
  * There is one Geometry block per element block.
  *
  * If both indexers provide their device local ID table (GlobalIndexer::getLIDs(),
  * built by DOFManager and Filtered_GlobalIndexer), the blocks hold local indices,
  * gathered on device: row indices are local indices of \c rowIndexer's
  * owned-and-ghosted indices (the factory's ghosted row map), and column indices
  * are local indices of \c colIndexer's owned-and-ghosted indices (its ghosted
  * column map).  So the ghosted column map must be passed to
  * assembleFECrsGraphWithMaps as the owned+shared domain map (and no column map),
  * which is what the factories do.  Otherwise (an indexer without a local ID table,
  * or associated neighbors), the blocks hold global indices, gathered on host
  * with GlobalIndexer::getElementGIDs.
  *
  * \param rowIndexer [in] Provides each element's row DOFs.
  * \param colIndexer [in] Provides each element's column DOFs.
  * \param elementIndexer [in] Provides the element blocks and their local elements.
  * \param connMgr [in] If not null and it has associated neighbors, each element's
  *                     columns include its associated neighbors' column DOFs.
  */
template <typename LocalOrdinalT, typename GlobalOrdinalT, typename NodeT>
Tpetra::Geometry<LocalOrdinalT, GlobalOrdinalT, NodeT>
buildFEGraphGeometry(const GlobalIndexer& rowIndexer,
                     const GlobalIndexer& colIndexer,
                     const GlobalIndexer& elementIndexer,
                     const Teuchos::RCP<const ConnManager>& connMgr)
{
  using geometry_type   = Tpetra::Geometry<LocalOrdinalT, GlobalOrdinalT, NodeT>;
  using device_type     = typename NodeT::device_type;
  using execution_space = typename NodeT::execution_space;

  static_assert(Kokkos::SpaceAccessibility<execution_space, PHX::Device::memory_space>::accessible,
                "panzer::buildFEGraphGeometry: the indexers' local ID tables must be accessible "
                "from the Tpetra node's execution space.");

  const bool square            = &rowIndexer == &colIndexer;
  const bool associatedNeighbors = !connMgr.is_null() && connMgr->hasAssociatedNeighbors();

  std::vector<std::string> blockIds;
  elementIndexer.getElementBlockIds(blockIds);

  // Use the local ID tables if both indexers have them, covering every element.
  auto rowLIDs = rowIndexer.getLIDs();
  auto colLIDs = colIndexer.getLIDs();
  bool useLIDs = !associatedNeighbors;
  for (const std::string& blockId : blockIds) {
    const std::vector<panzer::LocalOrdinal>& elements = elementIndexer.getElementBlock(blockId);
    panzer::LocalOrdinal maxElement = -1;
    for (const panzer::LocalOrdinal e : elements) maxElement = std::max(maxElement, e);
    if (maxElement >= 0 && (size_t(maxElement) >= rowLIDs.extent(0) || size_t(maxElement) >= colLIDs.extent(0)))
      useLIDs = false;
    if (rowLIDs.extent(1) < size_t(rowIndexer.getElementBlockGIDCount(blockId)) ||
        colLIDs.extent(1) < size_t(colIndexer.getElementBlockGIDCount(blockId)))
      useLIDs = false;
  }

  geometry_type geometry;
  for (const std::string& blockId : blockIds) {
    const std::vector<panzer::LocalOrdinal>& elements = elementIndexer.getElementBlock(blockId);
    const size_t numElements = elements.size();

    if (useLIDs) {
      // Gather the block's elements' rows of the local ID tables, on device.
      Kokkos::View<panzer::LocalOrdinal*, device_type> elementsDevice(
          Kokkos::ViewAllocateWithoutInitializing("panzer::buildFEGraphGeometry::elements"), numElements);
      Kokkos::deep_copy(elementsDevice,
                        Kokkos::View<const panzer::LocalOrdinal*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                            elements.data(), numElements));
      auto gather = [&](const decltype(rowLIDs)& lids, const size_t numDofs) {
        Kokkos::View<LocalOrdinalT**, Kokkos::LayoutRight, device_type> block(
            Kokkos::ViewAllocateWithoutInitializing("panzer::buildFEGraphGeometry::block"), numElements, numDofs);
        Kokkos::parallel_for(
            "panzer::buildFEGraphGeometry::gather",
            Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<2>>({0, 0}, {numElements, numDofs}),
            KOKKOS_LAMBDA(const size_t e, const size_t k) {
              block(e, k) = static_cast<LocalOrdinalT>(lids(elementsDevice(e), k));
            });
        return block;
      };
      auto rows = gather(rowLIDs, rowIndexer.getElementBlockGIDCount(blockId));
      if (square)
        geometry.addLocalBlock(rows);
      else
        geometry.addLocalBlock(rows, gather(colLIDs, colIndexer.getElementBlockGIDCount(blockId)));
    } else {
      // Gather the block's elements' global indices on host.  Elements may have
      // different numbers of columns (associated neighbors): pad with invalid.
      const GlobalOrdinalT INVALID = Tpetra::Details::OrdinalTraits<GlobalOrdinalT>::invalid();
      std::vector<std::vector<GlobalOrdinalT>> rowGIDs(numElements), colGIDs(numElements);
      size_t maxRows = 0, maxCols = 0;
      std::vector<GlobalOrdinalT> neighborGIDs;
      for (size_t e = 0; e < numElements; ++e) {
        rowIndexer.getElementGIDs(elements[e], rowGIDs[e], blockId);
        colIndexer.getElementGIDs(elements[e], colGIDs[e], blockId);
        if (associatedNeighbors) {
          for (const panzer::LocalOrdinal neighbor : connMgr->getAssociatedNeighbors(elements[e])) {
            colIndexer.getElementGIDs(neighbor, neighborGIDs);
            colGIDs[e].insert(colGIDs[e].end(), neighborGIDs.begin(), neighborGIDs.end());
          }
        }
        maxRows = std::max(maxRows, rowGIDs[e].size());
        maxCols = std::max(maxCols, colGIDs[e].size());
      }
      auto toDevice = [&](const std::vector<std::vector<GlobalOrdinalT>>& gids, const size_t width) {
        Kokkos::View<GlobalOrdinalT**, Kokkos::LayoutRight, device_type> block(
            Kokkos::ViewAllocateWithoutInitializing("panzer::buildFEGraphGeometry::block"), numElements, width);
        auto blockHost = Kokkos::create_mirror_view(block);
        for (size_t e = 0; e < numElements; ++e)
          for (size_t k = 0; k < width; ++k)
            blockHost(e, k) = k < gids[e].size() ? gids[e][k] : INVALID;
        Kokkos::deep_copy(block, blockHost);
        return block;
      };
      auto rows = toDevice(rowGIDs, maxRows);
      if (square && !associatedNeighbors)
        geometry.addGlobalBlock(rows);
      else
        geometry.addGlobalBlock(rows, toDevice(colGIDs, maxCols));
    }
  }
  return geometry;
}

} // namespace panzer

#endif
