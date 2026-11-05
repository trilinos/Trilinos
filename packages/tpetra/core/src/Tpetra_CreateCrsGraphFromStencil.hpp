// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CREATECRSGRAPHFROMSTENCIL_HPP
#define TPETRA_CREATECRSGRAPHFROMSTENCIL_HPP

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_TestForException.hpp"

#include <Kokkos_Sort.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

namespace Tpetra {
namespace Experimental {

/// \brief A relative Cartesian grid position in at most three dimensions.
///
/// Components belonging to dimensions not present in the grid must be zero.
struct StencilOffset {
  int x = 0;
  int y = 0;
  int z = 0;
};

enum class StencilNumbering {
  globalLexicographic,
  rankBlocked
};

namespace Details {

enum class StencilAction : unsigned char {
  invalid = 0,
  local   = 1,
  remote  = 2
};

template <class GlobalOrdinal, class LocalNodes>
KOKKOS_INLINE_FUNCTION void getLocalNodeCoordinates(
    const GlobalOrdinal localNode, const LocalNodes& localNodes,
    GlobalOrdinal& x, GlobalOrdinal& y, GlobalOrdinal& z) {
  x = localNode % localNodes[0];
  y = (localNode / localNodes[0]) % localNodes[1];
  z = localNode / (localNodes[0] * localNodes[1]);
}

template <class GlobalOrdinal>
KOKKOS_INLINE_FUNCTION size_t getBoundaryRegionCoordinate(
    const GlobalOrdinal coordinate, const GlobalOrdinal extent,
    const int radius) {
  if (extent <= 2 * radius)
    return static_cast<size_t>(coordinate);
  if (coordinate < radius)
    return static_cast<size_t>(coordinate);
  if (coordinate >= extent - radius)
    return static_cast<size_t>(radius + 1 + coordinate - (extent - radius));
  return static_cast<size_t>(radius);
}

template <class GlobalOrdinal>
KOKKOS_INLINE_FUNCTION void getBoundaryRegionGeometry(
    const size_t regionCoordinate, const GlobalOrdinal extent,
    const int radius, GlobalOrdinal& start, size_t& count) {
  if (extent <= 2 * radius || regionCoordinate < static_cast<size_t>(radius)) {
    start = static_cast<GlobalOrdinal>(regionCoordinate);
    count = 1;
  } else if (regionCoordinate == static_cast<size_t>(radius)) {
    start = radius;
    count = static_cast<size_t>(extent - 2 * radius);
  } else {
    start = extent - radius +
            static_cast<GlobalOrdinal>(regionCoordinate - radius - 1);
    count = 1;
  }
}

template <class GlobalOrdinal, class ProcGrid, class RankData>
KOKKOS_INLINE_FUNCTION void updateNeighborDimensions(
    const int procCoordinates[3], const ProcGrid& procGrid,
    const RankData& rankData, int& neighborRank,
    GlobalOrdinal dimensions[3]) {
  neighborRank = procCoordinates[2] * procGrid[0] * procGrid[1] +
                 procCoordinates[1] * procGrid[0] + procCoordinates[0];
  dimensions[0] = rankData(4 * neighborRank);
  dimensions[1] = rankData(4 * neighborRank + 1);
  dimensions[2] = rankData(4 * neighborRank + 2);
}

template <class GlobalOrdinal, class LocalNodes, class GlobalNodes,
          class ProcGrid, class RankData>
KOKKOS_INLINE_FUNCTION bool resolveNeighbor(
    const GlobalOrdinal x, const GlobalOrdinal y, const GlobalOrdinal z,
    const StencilOffset& offset, const LocalNodes& localNodes,
    const GlobalNodes& globalNodes, const GlobalNodes& globalNodeStart,
    const ProcGrid& procGrid, const RankData& rankData,
    const int myProcX, const int myProcY, const int myProcZ,
    const StencilNumbering numbering,
    int& neighborRank, GlobalOrdinal& neighborNode) {
  GlobalOrdinal coordinates[3] = {
      x + static_cast<GlobalOrdinal>(offset.x),
      y + static_cast<GlobalOrdinal>(offset.y),
      z + static_cast<GlobalOrdinal>(offset.z)};
  const GlobalOrdinal globalCoordinates[3] = {
      globalNodeStart[0] + coordinates[0],
      globalNodeStart[1] + coordinates[1],
      globalNodeStart[2] + coordinates[2]};
  for (int dim = 0; dim < 3; ++dim)
    if (globalCoordinates[dim] < 0 ||
        globalCoordinates[dim] >= globalNodes[dim])
      return false;
  GlobalOrdinal dimensions[3] = {
      localNodes[0], localNodes[1], localNodes[2]};

  int procCoordinates[3] = {myProcX, myProcY, myProcZ};

  for (int dim = 0; dim < 3; ++dim) {
    if (coordinates[dim] < 0) {
      --procCoordinates[dim];
      if (procCoordinates[dim] < 0)
        return false;
      updateNeighborDimensions(procCoordinates, procGrid, rankData,
                               neighborRank, dimensions);
      coordinates[dim] += dimensions[dim];
    }
    if (coordinates[dim] < 0)
      return false;
    if (coordinates[dim] >= dimensions[dim]) {
      coordinates[dim] -= dimensions[dim];
      ++procCoordinates[dim];
      if (procCoordinates[dim] >= procGrid[dim])
        return false;
      updateNeighborDimensions(procCoordinates, procGrid, rankData,
                               neighborRank, dimensions);
    }
    if (coordinates[dim] >= dimensions[dim])
      return false;
  }

  updateNeighborDimensions(procCoordinates, procGrid, rankData, neighborRank,
                           dimensions);
  if (numbering == StencilNumbering::globalLexicographic) {
    neighborNode = globalCoordinates[2] * globalNodes[0] * globalNodes[1] +
                   globalCoordinates[1] * globalNodes[0] + globalCoordinates[0];
  } else {
    neighborNode = rankData(4 * neighborRank + 3) +
                   coordinates[2] * dimensions[0] * dimensions[1] +
                   coordinates[1] * dimensions[0] + coordinates[0];
  }
  return true;
}

}  // namespace Details

/// \brief Create a fill-complete graph directly in device-resident CSR form.
/// This implementation supports radius-two stencils, at most three
/// dimensions, and at most three degrees of freedom per grid point.
template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>
createCrsGraphFromStencil(
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
    const Teuchos::ArrayView<const LocalOrdinal>& localGridDimensions,
    const Teuchos::ArrayView<const StencilOffset>& stencilOffsets,
    const LocalOrdinal dofsPerNode, const StencilNumbering numbering,
    const bool debug) {
  using graph_type       = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using import_type      = Import<LocalOrdinal, GlobalOrdinal, Node>;
  using local_graph_type = typename graph_type::local_graph_device_type;
  using row_map_type     = typename local_graph_type::row_map_type::non_const_type;
  using entries_type     = typename local_graph_type::entries_type::non_const_type;
  using device_type      = typename Node::device_type;
  using execution_space  = typename device_type::execution_space;
  using range_policy =
      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<size_t>>;

  constexpr size_t maxStencilSize = 125;  // Full radius-two stencil in 3D.
  constexpr size_t maxDofsPerNode = 3;
  TEUCHOS_TEST_FOR_EXCEPTION(rowMap.is_null(), std::invalid_argument,
                             "createCrsGraphFromStencil: rowMap is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      localGridDimensions.size() < 1 || localGridDimensions.size() > 3,
      std::invalid_argument,
      "createCrsGraphFromStencil: the grid must have between one and three "
      "dimensions.");
  TEUCHOS_TEST_FOR_EXCEPTION(stencilOffsets.size() == 0, std::invalid_argument,
                             "createCrsGraphFromStencil: the stencil is empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      Teuchos::as<size_t>(stencilOffsets.size()) > maxStencilSize,
      std::invalid_argument,
      "createCrsGraphFromStencil: at most 125 stencil offsets are supported.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      dofsPerNode <= 0 || Teuchos::as<size_t>(dofsPerNode) > maxDofsPerNode,
      std::invalid_argument,
      "createCrsGraphFromStencil: dofsPerNode must be between one and three.");

  const int numDimensions = localGridDimensions.size();
  Kokkos::Array<GlobalOrdinal, 3> localNodes{{1, 1, 1}};
  GlobalOrdinal localNumNodesGO = 1;
  for (int dim = 0; dim < numDimensions; ++dim) {
    TEUCHOS_TEST_FOR_EXCEPTION(localGridDimensions[dim] <= 0,
                               std::invalid_argument,
                               "createCrsGraphFromStencil: local grid "
                               "dimensions must be positive.");
    localNodes[dim] = Teuchos::as<GlobalOrdinal>(localGridDimensions[dim]);
    localNumNodesGO *= localNodes[dim];
  }

  Kokkos::View<StencilOffset*, device_type> offsets(
      Kokkos::ViewAllocateWithoutInitializing("Tpetra stencil offsets"),
      stencilOffsets.size());
  auto offsetsHost     = Kokkos::create_mirror_view(offsets);
  int stencilRadius[3] = {0, 0, 0};
  for (size_t entry = 0; entry < Teuchos::as<size_t>(stencilOffsets.size());
       ++entry) {
    const auto offset = stencilOffsets[entry];
    TEUCHOS_TEST_FOR_EXCEPTION(
        offset.x < -2 || offset.x > 2 || offset.y < -2 || offset.y > 2 ||
            offset.z < -2 || offset.z > 2,
        std::invalid_argument,
        "createCrsGraphFromStencil: only radius-two offsets are supported.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        (numDimensions < 2 && offset.y != 0) ||
            (numDimensions < 3 && offset.z != 0),
        std::invalid_argument,
        "createCrsGraphFromStencil: an offset uses an inactive dimension.");
    const int magnitudes[3] = {
        offset.x < 0 ? -offset.x : offset.x,
        offset.y < 0 ? -offset.y : offset.y,
        offset.z < 0 ? -offset.z : offset.z};
    for (int dim = 0; dim < 3; ++dim)
      if (magnitudes[dim] > stencilRadius[dim])
        stencilRadius[dim] = magnitudes[dim];
    offsetsHost(entry) = offset;
  }
  Kokkos::deep_copy(offsets, offsetsHost);

  Teuchos::Array<size_t> interiorStencilOrder(stencilOffsets.size());
  for (size_t stencil = 0; stencil < Teuchos::as<size_t>(stencilOffsets.size());
       ++stencil)
    interiorStencilOrder[stencil] = stencil;
  const auto linearNodeOffset = [&](const size_t stencil) {
    const auto offset = stencilOffsets[stencil];
    return Teuchos::as<GlobalOrdinal>(offset.x) +
           localNodes[0] *
               (Teuchos::as<GlobalOrdinal>(offset.y) +
                localNodes[1] * Teuchos::as<GlobalOrdinal>(offset.z));
  };
  std::sort(interiorStencilOrder.begin(), interiorStencilOrder.end(),
            [&](const size_t left, const size_t right) {
              return linearNodeOffset(left) < linearNodeOffset(right);
            });
  Kokkos::View<GlobalOrdinal*, device_type> interiorNodeOffsets(
      Kokkos::ViewAllocateWithoutInitializing(
          "Tpetra interior stencil node offsets"),
      stencilOffsets.size());
  auto interiorNodeOffsetsHost = Kokkos::create_mirror_view(interiorNodeOffsets);
  for (size_t stencil = 0; stencil < Teuchos::as<size_t>(stencilOffsets.size());
       ++stencil)
    interiorNodeOffsetsHost(stencil) =
        linearNodeOffset(interiorStencilOrder[stencil]);
  Kokkos::deep_copy(interiorNodeOffsets, interiorNodeOffsetsHost);

  const GlobalOrdinal dofsPerNodeGO =
      Teuchos::as<GlobalOrdinal>(dofsPerNode);
  const size_t rowsPerNode  = Teuchos::as<size_t>(dofsPerNode);
  const size_t localNumRows = rowMap->getLocalNumElements();
  TEUCHOS_TEST_FOR_EXCEPTION(
      Teuchos::as<GlobalOrdinal>(localNumRows) !=
          localNumNodesGO * dofsPerNodeGO,
      std::invalid_argument,
      "createCrsGraphFromStencil: local grid dimensions and dofsPerNode do "
      "not match the local row count.");

  const GlobalOrdinal numGlobalRows =
      Teuchos::as<GlobalOrdinal>(rowMap->getGlobalNumElements());
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalRows % dofsPerNodeGO != 0,
                             std::invalid_argument,
                             "createCrsGraphFromStencil: global row count is "
                             "not divisible by dofsPerNode.");
  const GlobalOrdinal numGlobalNodes = numGlobalRows / dofsPerNodeGO;
  const GlobalOrdinal globalMinGid   = rowMap->getMinAllGlobalIndex();
  const GlobalOrdinal globalMaxGid   = rowMap->getMaxAllGlobalIndex();
  TEUCHOS_TEST_FOR_EXCEPTION(
      globalMaxGid - globalMinGid + 1 != numGlobalRows,
      std::invalid_argument,
      "createCrsGraphFromStencil: row-map GIDs must form a contiguous range.");

  const auto comm                 = rowMap->getComm();
  const int myRank                = comm->getRank();
  const int numRanks              = comm->getSize();
  const GlobalOrdinal localMinGid = rowMap->getMinGlobalIndex();
  TEUCHOS_TEST_FOR_EXCEPTION(
      (localMinGid - globalMinGid) % dofsPerNodeGO != 0,
      std::invalid_argument,
      "createCrsGraphFromStencil: a local row range does not begin at a "
      "grid-point boundary.");
  const GlobalOrdinal firstLocalNode =
      (localMinGid - globalMinGid) / dofsPerNodeGO;
  // Store dimensions followed by the rank's first node in row-map order.
  Teuchos::Array<GlobalOrdinal> localRankData(4);
  localRankData[0] = localNodes[0];
  localRankData[1] = localNodes[1];
  localRankData[2] = localNodes[2];
  localRankData[3] = firstLocalNode;
  Teuchos::Array<GlobalOrdinal> rankData(4 * numRanks);
  Teuchos::gatherAll(*comm, 4, localRankData.getRawPtr(), 4 * numRanks,
                     rankData.getRawPtr());

  // Infer the same Cartesian process grid used by Galeri's map construction.
  Kokkos::Array<int, 3> procGrid{{1, 1, 1}};
  if (numDimensions == 1) {
    procGrid[0] = numRanks;
  } else if (numDimensions == 2) {
    while ((procGrid[0] + 1) * (procGrid[0] + 1) <= numRanks)
      ++procGrid[0];
    procGrid[1] = numRanks / procGrid[0];
    while (procGrid[0] * procGrid[1] != numRanks)
      procGrid[1] = numRanks / (--procGrid[0]);
  } else {
    int cubeRoot = 1;
    while ((cubeRoot + 1) * (cubeRoot + 1) * (cubeRoot + 1) <= numRanks)
      ++cubeRoot;
    procGrid[0] = procGrid[1] = procGrid[2] = cubeRoot;
    if (procGrid[0] * procGrid[1] * procGrid[2] != numRanks) {
      procGrid[0] = procGrid[1] = procGrid[2] = 1;
      int remainingRanks                      = numRanks;
      constexpr int factorLimit               = 50;
      int factors[factorLimit]                = {};
      for (int factor = 2; factor < factorLimit; ++factor) {
        while (remainingRanks % factor == 0) {
          ++factors[factor];
          remainingRanks /= factor;
        }
      }
      procGrid[0] = remainingRanks;
      for (int factor = factorLimit - 1; factor > 0; --factor) {
        while (factors[factor] != 0) {
          if (procGrid[0] <= procGrid[1] && procGrid[0] <= procGrid[2])
            procGrid[0] *= factor;
          else if (procGrid[1] <= procGrid[0] &&
                   procGrid[1] <= procGrid[2])
            procGrid[1] *= factor;
          else
            procGrid[2] *= factor;
          --factors[factor];
        }
      }
    }
  }

  const int procXY  = procGrid[0] * procGrid[1];
  const int myProcX = myRank % procGrid[0];
  const int myProcY = (myRank % procXY) / procGrid[0];
  const int myProcZ = myRank / procXY;
  TEUCHOS_TEST_FOR_EXCEPTION(
      procGrid[0] <= 0 || procGrid[1] <= 0 || procGrid[2] <= 0 ||
          procXY * procGrid[2] != numRanks,
      std::runtime_error,
      "createCrsGraphFromStencil: failed to infer the process grid.");

  Kokkos::Array<GlobalOrdinal, 3> globalNodes{{0, 0, 0}};
  for (int px = 0; px < procGrid[0]; ++px)
    globalNodes[0] += rankData[4 * px];
  for (int py = 0; py < procGrid[1]; ++py)
    globalNodes[1] += rankData[4 * py * procGrid[0] + 1];
  for (int pz = 0; pz < procGrid[2]; ++pz)
    globalNodes[2] += rankData[4 * pz * procXY + 2];
  Kokkos::Array<GlobalOrdinal, 3> globalNodeStart{{0, 0, 0}};
  for (int px = 0; px < myProcX; ++px)
    globalNodeStart[0] += rankData[4 * px];
  for (int py = 0; py < myProcY; ++py)
    globalNodeStart[1] += rankData[4 * py * procGrid[0] + 1];
  for (int pz = 0; pz < myProcZ; ++pz)
    globalNodeStart[2] += rankData[4 * pz * procXY + 2];
  TEUCHOS_TEST_FOR_EXCEPTION(
      globalNodes[0] * globalNodes[1] * globalNodes[2] != numGlobalNodes,
      std::invalid_argument,
      "createCrsGraphFromStencil: process-grid dimensions do not match the "
      "global row map.");

  Kokkos::Array<size_t, 3> regionDimensions{{1, 1, 1}};
  for (int dim = 0; dim < 3; ++dim)
    regionDimensions[dim] = localNodes[dim] <= 2 * stencilRadius[dim]
                                ? Teuchos::as<size_t>(localNodes[dim])
                                : Teuchos::as<size_t>(2 * stencilRadius[dim] + 1);
  const size_t numRegions = regionDimensions[0] * regionDimensions[1] *
                            regionDimensions[2];
  const size_t stencilSize = stencilOffsets.size();
  Teuchos::Array<Details::StencilAction> planActionsHost(
      numRegions * stencilSize, Details::StencilAction::invalid);
  Teuchos::Array<size_t> planRowLengthsHost(numRegions, 0);
  for (size_t region = 0; region < numRegions; ++region) {
    const size_t regionX = region % regionDimensions[0];
    const size_t regionY =
        (region / regionDimensions[0]) % regionDimensions[1];
    const size_t regionZ =
        region / (regionDimensions[0] * regionDimensions[1]);
    GlobalOrdinal x = 0, y = 0, z = 0;
    size_t unusedCount = 0;
    Details::getBoundaryRegionGeometry(
        regionX, localNodes[0], stencilRadius[0], x, unusedCount);
    Details::getBoundaryRegionGeometry(
        regionY, localNodes[1], stencilRadius[1], y, unusedCount);
    Details::getBoundaryRegionGeometry(
        regionZ, localNodes[2], stencilRadius[2], z, unusedCount);
    for (size_t stencil = 0; stencil < stencilSize; ++stencil) {
      const auto offset                   = stencilOffsets[stencil];
      const GlobalOrdinal neighborGlobalX = globalNodeStart[0] + x + offset.x;
      const GlobalOrdinal neighborGlobalY = globalNodeStart[1] + y + offset.y;
      const GlobalOrdinal neighborGlobalZ = globalNodeStart[2] + z + offset.z;
      if (neighborGlobalX < 0 || neighborGlobalX >= globalNodes[0] ||
          neighborGlobalY < 0 || neighborGlobalY >= globalNodes[1] ||
          neighborGlobalZ < 0 || neighborGlobalZ >= globalNodes[2])
        continue;
      const GlobalOrdinal neighborLocalX = x + offset.x;
      const GlobalOrdinal neighborLocalY = y + offset.y;
      const GlobalOrdinal neighborLocalZ = z + offset.z;
      const bool local                   = neighborLocalX >= 0 && neighborLocalX < localNodes[0] &&
                         neighborLocalY >= 0 && neighborLocalY < localNodes[1] &&
                         neighborLocalZ >= 0 && neighborLocalZ < localNodes[2];
      planActionsHost[region * stencilSize + stencil] =
          local ? Details::StencilAction::local
                : Details::StencilAction::remote;
      planRowLengthsHost[region] += rowsPerNode;
    }
  }

  Kokkos::Array<GlobalOrdinal, 5> regionStartsX{{0, 0, 0, 0, 0}};
  Kokkos::Array<GlobalOrdinal, 5> regionStartsY{{0, 0, 0, 0, 0}};
  Kokkos::Array<GlobalOrdinal, 5> regionStartsZ{{0, 0, 0, 0, 0}};
  Kokkos::Array<size_t, 5> regionCountsX{{0, 0, 0, 0, 0}};
  Kokkos::Array<size_t, 5> regionCountsY{{0, 0, 0, 0, 0}};
  Kokkos::Array<size_t, 5> regionCountsZ{{0, 0, 0, 0, 0}};
  for (size_t coordinate = 0; coordinate < regionDimensions[0]; ++coordinate)
    Details::getBoundaryRegionGeometry(
        coordinate, localNodes[0], stencilRadius[0], regionStartsX[coordinate],
        regionCountsX[coordinate]);
  for (size_t coordinate = 0; coordinate < regionDimensions[1]; ++coordinate)
    Details::getBoundaryRegionGeometry(
        coordinate, localNodes[1], stencilRadius[1], regionStartsY[coordinate],
        regionCountsY[coordinate]);
  for (size_t coordinate = 0; coordinate < regionDimensions[2]; ++coordinate)
    Details::getBoundaryRegionGeometry(
        coordinate, localNodes[2], stencilRadius[2], regionStartsZ[coordinate],
        regionCountsZ[coordinate]);

  Teuchos::Array<size_t> regionXPrefixesHost(numRegions, 0);
  Teuchos::Array<size_t> regionYPrefixesHost(
      regionDimensions[1] * regionDimensions[2], 0);
  Teuchos::Array<size_t> regionZPrefixesHost(regionDimensions[2], 0);
  Teuchos::Array<size_t> lineLengthsHost(
      regionDimensions[1] * regionDimensions[2], 0);
  Teuchos::Array<size_t> planeLengthsHost(regionDimensions[2], 0);
  size_t totalRowLengths = 0;
  for (size_t regionZ = 0; regionZ < regionDimensions[2]; ++regionZ) {
    regionZPrefixesHost[regionZ] = totalRowLengths;
    size_t planeLength           = 0;
    for (size_t regionY = 0; regionY < regionDimensions[1]; ++regionY) {
      const size_t line         = regionY + regionDimensions[1] * regionZ;
      regionYPrefixesHost[line] = planeLength;
      size_t lineLength         = 0;
      for (size_t regionX = 0; regionX < regionDimensions[0]; ++regionX) {
        const size_t region         = regionX + regionDimensions[0] * line;
        regionXPrefixesHost[region] = lineLength;
        lineLength += regionCountsX[regionX] *
                      planRowLengthsHost[region];
      }
      lineLengthsHost[line] = lineLength;
      planeLength += regionCountsY[regionY] * lineLength;
    }
    planeLengthsHost[regionZ] = planeLength;
    totalRowLengths += regionCountsZ[regionZ] * planeLength;
  }

  execution_space executionSpace;
  Kokkos::View<Details::StencilAction*, device_type> planActions(
      "Tpetra stencil plan actions", planActionsHost.size());
  Kokkos::View<size_t*, device_type> planRowLengths(
      "Tpetra stencil plan row lengths", numRegions);
  Kokkos::View<size_t*, device_type> regionXPrefixes(
      "Tpetra stencil region x prefixes", numRegions);
  Kokkos::View<size_t*, device_type> regionYPrefixes(
      "Tpetra stencil region y prefixes", regionYPrefixesHost.size());
  Kokkos::View<size_t*, device_type> regionZPrefixes(
      "Tpetra stencil region z prefixes", regionZPrefixesHost.size());
  Kokkos::View<size_t*, device_type> lineLengths(
      "Tpetra stencil region line lengths", lineLengthsHost.size());
  Kokkos::View<size_t*, device_type> planeLengths(
      "Tpetra stencil region plane lengths", planeLengthsHost.size());
  const auto copyHostArray = [&](auto destination, auto& source) {
    using value_type = typename std::remove_reference_t<decltype(source)>::value_type;
    Kokkos::deep_copy(
        executionSpace, destination,
        Kokkos::View<value_type*, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            source.getRawPtr(), source.size()));
  };
  copyHostArray(planActions, planActionsHost);
  copyHostArray(planRowLengths, planRowLengthsHost);
  copyHostArray(regionXPrefixes, regionXPrefixesHost);
  copyHostArray(regionYPrefixes, regionYPrefixesHost);
  copyHostArray(regionZPrefixes, regionZPrefixesHost);
  copyHostArray(lineLengths, lineLengthsHost);
  copyHostArray(planeLengths, planeLengthsHost);
  Kokkos::View<GlobalOrdinal*, device_type> rankDataDevice(
      Kokkos::ViewAllocateWithoutInitializing("Tpetra stencil rank metadata"),
      rankData.size());
  auto rankDataHost = Kokkos::create_mirror_view(rankDataDevice);
  for (int entry = 0; entry < rankData.size(); ++entry)
    rankDataHost(entry) = rankData[entry];
  Kokkos::deep_copy(executionSpace, rankDataDevice, rankDataHost);

  row_map_type rowptr(
      Kokkos::ViewAllocateWithoutInitializing("Tpetra stencil row pointers"),
      localNumRows + 1);
  const size_t localNnz      = totalRowLengths * rowsPerNode;
  const size_t localNumNodes = localNumRows / rowsPerNode;
  size_t expandedNodes       = 1;
  for (int dim = 0; dim < 3; ++dim)
    expandedNodes *= Teuchos::as<size_t>(
        localNodes[dim] + 2 * stencilRadius[dim]);
  const size_t maxRemoteNodes =
      numRanks == 1 ? 0 : expandedNodes - localNumNodes;
  Kokkos::UnorderedMap<GlobalOrdinal, int, device_type> remoteNodes(
      maxRemoteNodes == 0 ? 1 : maxRemoteNodes);

  Teuchos::Array<GlobalOrdinal> haloCoordinatesHost(3 * maxRemoteNodes);
  size_t haloEntry                = 0;
  const auto appendHaloCoordinate = [&](const GlobalOrdinal x,
                                        const GlobalOrdinal y,
                                        const GlobalOrdinal z) {
    haloCoordinatesHost[3 * haloEntry]     = x;
    haloCoordinatesHost[3 * haloEntry + 1] = y;
    haloCoordinatesHost[3 * haloEntry + 2] = z;
    ++haloEntry;
  };
  const GlobalOrdinal radiusX = stencilRadius[0];
  const GlobalOrdinal radiusY = stencilRadius[1];
  const GlobalOrdinal radiusZ = stencilRadius[2];
  if (numRanks > 1) {
    for (GlobalOrdinal z = -radiusZ; z < localNodes[2] + radiusZ; ++z)
      for (GlobalOrdinal y = -radiusY; y < localNodes[1] + radiusY; ++y) {
        for (GlobalOrdinal x = -radiusX; x < 0; ++x)
          appendHaloCoordinate(x, y, z);
        for (GlobalOrdinal x = localNodes[0];
             x < localNodes[0] + radiusX; ++x)
          appendHaloCoordinate(x, y, z);
      }
    for (GlobalOrdinal z = -radiusZ; z < localNodes[2] + radiusZ; ++z)
      for (GlobalOrdinal x = 0; x < localNodes[0]; ++x) {
        for (GlobalOrdinal y = -radiusY; y < 0; ++y)
          appendHaloCoordinate(x, y, z);
        for (GlobalOrdinal y = localNodes[1];
             y < localNodes[1] + radiusY; ++y)
          appendHaloCoordinate(x, y, z);
      }
    for (GlobalOrdinal y = 0; y < localNodes[1]; ++y)
      for (GlobalOrdinal x = 0; x < localNodes[0]; ++x) {
        for (GlobalOrdinal z = -radiusZ; z < 0; ++z)
          appendHaloCoordinate(x, y, z);
        for (GlobalOrdinal z = localNodes[2];
             z < localNodes[2] + radiusZ; ++z)
          appendHaloCoordinate(x, y, z);
      }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
      haloEntry != maxRemoteNodes, std::logic_error,
      "createCrsGraphFromStencil: internal halo-shell size mismatch.");
  Kokkos::View<GlobalOrdinal*, device_type> haloCoordinates(
      Kokkos::ViewAllocateWithoutInitializing("Tpetra stencil halo coordinates"),
      haloCoordinatesHost.size());
  Kokkos::deep_copy(
      executionSpace, haloCoordinates,
      Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
          haloCoordinatesHost.getRawPtr(), haloCoordinatesHost.size()));
  if (maxRemoteNodes != 0) {
    Kokkos::parallel_for(
        "Tpetra stencil collect halo nodes",
        range_policy(executionSpace, 0, maxRemoteNodes),
        KOKKOS_LAMBDA(const size_t candidate) {
          const GlobalOrdinal neighborX = haloCoordinates(3 * candidate);
          const GlobalOrdinal neighborY = haloCoordinates(3 * candidate + 1);
          const GlobalOrdinal neighborZ = haloCoordinates(3 * candidate + 2);
          const GlobalOrdinal globalX   = globalNodeStart[0] + neighborX;
          const GlobalOrdinal globalY   = globalNodeStart[1] + neighborY;
          const GlobalOrdinal globalZ   = globalNodeStart[2] + neighborZ;
          if (globalX < 0 || globalX >= globalNodes[0] ||
              globalY < 0 || globalY >= globalNodes[1] ||
              globalZ < 0 || globalZ >= globalNodes[2])
            return;
          for (size_t stencil = 0; stencil < stencilSize; ++stencil) {
            const GlobalOrdinal sourceX = neighborX - offsets[stencil].x;
            const GlobalOrdinal sourceY = neighborY - offsets[stencil].y;
            const GlobalOrdinal sourceZ = neighborZ - offsets[stencil].z;
            if (sourceX < 0 || sourceX >= localNodes[0] ||
                sourceY < 0 || sourceY >= localNodes[1] ||
                sourceZ < 0 || sourceZ >= localNodes[2])
              continue;
            int neighborRank           = myRank;
            GlobalOrdinal neighborNode = 0;
            const bool valid           = Details::resolveNeighbor(
                          sourceX, sourceY, sourceZ, offsets[stencil], localNodes,
                          globalNodes, globalNodeStart, procGrid, rankDataDevice,
                          myProcX, myProcY, myProcZ, numbering, neighborRank,
                          neighborNode);
            if (valid) {
              remoteNodes.insert(neighborNode, neighborRank);
              break;
            }
          }
        });
    executionSpace.fence("Tpetra stencil halo nodes ready");
    TEUCHOS_TEST_FOR_EXCEPTION(
        remoteNodes.failed_insert(), std::runtime_error,
        "createCrsGraphFromStencil: remote-node halo capacity was insufficient.");
  }

  const size_t numRemoteNodes =
      maxRemoteNodes == 0 ? 0 : remoteNodes.size();
  Kokkos::View<GlobalOrdinal*, device_type> remoteNodeGids(
      Kokkos::ViewAllocateWithoutInitializing("Tpetra stencil remote nodes"),
      numRemoteNodes);
  Kokkos::View<int*, device_type> remoteNodePids(
      Kokkos::ViewAllocateWithoutInitializing("Tpetra stencil remote node PIDs"),
      numRemoteNodes);
  if (numRemoteNodes != 0) {
    Kokkos::parallel_scan(
        "Tpetra stencil compact remote nodes",
        range_policy(executionSpace, 0, remoteNodes.capacity()),
        KOKKOS_LAMBDA(const size_t slot, size_t& output, const bool final) {
          if (remoteNodes.valid_at(slot)) {
            if (final) {
              remoteNodeGids(output) = remoteNodes.key_at(slot);
              remoteNodePids(output) = remoteNodes.value_at(slot);
            }
            ++output;
          }
        });
  }
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> colMap = rowMap;
  Teuchos::RCP<const import_type> importer;
  if (numRemoteNodes != 0) {
    Kokkos::View<size_t*, device_type> permutation(
        Kokkos::ViewAllocateWithoutInitializing(
            "Tpetra stencil halo permutation"),
        numRemoteNodes);
    Kokkos::parallel_for(
        "Tpetra stencil initialize halo permutation",
        range_policy(executionSpace, 0, numRemoteNodes),
        KOKKOS_LAMBDA(const size_t entry) { permutation(entry) = entry; });
    Kokkos::sort(
        executionSpace, permutation,
        KOKKOS_LAMBDA(const size_t left, const size_t right) {
          return remoteNodePids(left) != remoteNodePids(right)
                     ? remoteNodePids(left) < remoteNodePids(right)
                     : remoteNodeGids(left) < remoteNodeGids(right);
        });

    const size_t numRemoteColumns = numRemoteNodes * rowsPerNode;
    Kokkos::View<GlobalOrdinal*, device_type> columnGids(
        Kokkos::ViewAllocateWithoutInitializing("Tpetra stencil column map"),
        localNumRows + numRemoteColumns);
    Kokkos::View<int*, device_type> orderedRemoteNodePids(
        Kokkos::ViewAllocateWithoutInitializing(
            "Tpetra stencil ordered remote node PIDs"),
        numRemoteNodes);
    const auto localRowGids = rowMap->getMyGlobalIndicesDevice();
    Kokkos::parallel_for(
        "Tpetra stencil fill column map",
        range_policy(executionSpace, 0, localNumRows + numRemoteColumns),
        KOKKOS_LAMBDA(const size_t column) {
          if (column < localNumRows) {
            columnGids(column) = localRowGids(column);
            return;
          }
          const size_t remoteColumn = column - localNumRows;
          const size_t remoteNode   = remoteColumn / rowsPerNode;
          const size_t columnDof    = remoteColumn % rowsPerNode;
          const size_t source       = permutation(remoteNode);
          columnGids(column)        = globalMinGid +
                               remoteNodeGids(source) * dofsPerNodeGO +
                               static_cast<GlobalOrdinal>(columnDof);
          if (columnDof == 0)
            orderedRemoteNodePids(remoteNode) = remoteNodePids(source);
        });
    executionSpace.fence("Tpetra stencil column map ready");

    auto mapParams = Teuchos::rcp(new Teuchos::ParameterList);
    mapParams->set("compute global constants", false);
    const auto invalidGlobalSize =
        Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    colMap = Teuchos::rcp(new Map<LocalOrdinal, GlobalOrdinal, Node>(
        invalidGlobalSize, columnGids, rowMap->getIndexBase(), comm,
        mapParams));

    Teuchos::Array<int> remoteNodePIDs(numRemoteNodes);
    auto remoteNodePIDsHost =
        Kokkos::View<int*, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            remoteNodePIDs.getRawPtr(), numRemoteNodes);
    Kokkos::deep_copy(executionSpace, remoteNodePIDsHost,
                      orderedRemoteNodePids);
    executionSpace.fence("Tpetra stencil remote PIDs ready");
    Teuchos::Array<int> remotePIDs(numRemoteColumns);
    for (size_t node = 0; node < numRemoteNodes; ++node)
      for (size_t dof = 0; dof < rowsPerNode; ++dof)
        remotePIDs[node * rowsPerNode + dof] = remoteNodePIDs[node];
    importer = Teuchos::rcp(new import_type(rowMap, colMap, remotePIDs));
  }

  entries_type localColumns(
      Kokkos::ViewAllocateWithoutInitializing(
          "Tpetra stencil local column indices"),
      localNnz);
  const auto colLocalMap = colMap->getLocalMap();
  Kokkos::parallel_for(
      "Tpetra stencil fill local column indices",
      range_policy(executionSpace, 0, localNumNodes + 1),
      KOKKOS_LAMBDA(const size_t node) {
        if (node == localNumNodes) {
          rowptr(localNumRows) = localNnz;
          return;
        }
        const GlobalOrdinal localNode = static_cast<GlobalOrdinal>(node);
        GlobalOrdinal x = 0, y = 0, z = 0;
        Details::getLocalNodeCoordinates(localNode, localNodes, x, y, z);
        const size_t regionX = Details::getBoundaryRegionCoordinate(
            x, localNodes[0], stencilRadius[0]);
        const size_t regionY = Details::getBoundaryRegionCoordinate(
            y, localNodes[1], stencilRadius[1]);
        const size_t regionZ = Details::getBoundaryRegionCoordinate(
            z, localNodes[2], stencilRadius[2]);
        const size_t line             = regionY + regionDimensions[1] * regionZ;
        const size_t region           = regionX + regionDimensions[0] * line;
        const size_t plannedRowLength = planRowLengths(region);
        const size_t precedingRowLengths =
            regionZPrefixes(regionZ) +
            static_cast<size_t>(z - regionStartsZ[regionZ]) *
                planeLengths(regionZ) +
            regionYPrefixes(line) +
            static_cast<size_t>(y - regionStartsY[regionY]) *
                lineLengths(line) +
            regionXPrefixes(region) +
            static_cast<size_t>(x - regionStartsX[regionX]) *
                plannedRowLength;
        const size_t firstRow = node * rowsPerNode;
        const size_t begin    = rowsPerNode * precedingRowLengths;
        for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof)
          rowptr(firstRow + rowDof) = begin + rowDof * plannedRowLength;

        const bool interior =
            x >= stencilRadius[0] &&
            x + stencilRadius[0] < localNodes[0] &&
            y >= stencilRadius[1] &&
            y + stencilRadius[1] < localNodes[1] &&
            z >= stencilRadius[2] &&
            z + stencilRadius[2] < localNodes[2];
        if (interior) {
          size_t output = begin;
          for (size_t stencil = 0; stencil < stencilSize; ++stencil) {
            const LocalOrdinal firstColumn = static_cast<LocalOrdinal>(
                (localNode + interiorNodeOffsets(stencil)) * dofsPerNodeGO);
            for (size_t columnDof = 0; columnDof < rowsPerNode; ++columnDof)
              localColumns(output++) =
                  firstColumn + static_cast<LocalOrdinal>(columnDof);
          }
          const size_t rowLength = output - begin;
          for (size_t rowDof = 1; rowDof < rowsPerNode; ++rowDof) {
            const size_t rowBegin = rowptr(firstRow + rowDof);
            for (size_t entry = 0; entry < rowLength; ++entry)
              localColumns(rowBegin + entry) = localColumns(begin + entry);
          }
          return;
        }

        size_t output = begin;
        for (size_t stencil = 0; stencil < stencilSize; ++stencil) {
          const Details::StencilAction action =
              planActions(region * stencilSize + stencil);
          if (action == Details::StencilAction::invalid)
            continue;

          const GlobalOrdinal neighborLocalX = x + offsets[stencil].x;
          const GlobalOrdinal neighborLocalY = y + offsets[stencil].y;
          const GlobalOrdinal neighborLocalZ = z + offsets[stencil].z;
          int neighborRank                   = myRank;
          GlobalOrdinal neighborNode         = 0;
          LocalOrdinal firstColumn           = 0;
          if (action == Details::StencilAction::local) {
            const GlobalOrdinal neighborLocalNode =
                neighborLocalZ * localNodes[0] * localNodes[1] +
                neighborLocalY * localNodes[0] + neighborLocalX;
            firstColumn = static_cast<LocalOrdinal>(
                neighborLocalNode * dofsPerNodeGO);
          } else {
            const bool valid = Details::resolveNeighbor(
                x, y, z, offsets[stencil], localNodes,
                globalNodes, globalNodeStart, procGrid, rankDataDevice,
                myProcX, myProcY, myProcZ, numbering, neighborRank,
                neighborNode);
            if (debug)
              KOKKOS_ASSERT(valid);
            (void)valid;
            firstColumn = colLocalMap.getLocalElement(
                globalMinGid + neighborNode * dofsPerNodeGO);
          }
          for (size_t columnDof = 0; columnDof < rowsPerNode; ++columnDof)
            localColumns(output++) =
                firstColumn + static_cast<LocalOrdinal>(columnDof);
        }
        for (size_t entry = begin + 1; entry < output; ++entry) {
          const LocalOrdinal value = localColumns(entry);
          size_t insertion         = entry;
          while (insertion > begin && value < localColumns(insertion - 1)) {
            localColumns(insertion) = localColumns(insertion - 1);
            --insertion;
          }
          localColumns(insertion) = value;
        }
        const size_t rowLength = output - begin;
        for (size_t rowDof = 1; rowDof < rowsPerNode; ++rowDof) {
          const size_t rowBegin = rowptr(firstRow + rowDof);
          for (size_t entry = 0; entry < rowLength; ++entry)
            localColumns(rowBegin + entry) = localColumns(begin + entry);
        }
      });
  executionSpace.fence("Tpetra stencil local columns ready");

  auto params = Teuchos::rcp(new Teuchos::ParameterList);
  params->set("No Nonlocal Changes", true);
  params->set("Optimize Storage", true);
  params->set("compute global constants", true);
  local_graph_type localGraph(localColumns, rowptr);
  return Teuchos::rcp(
      new graph_type(localGraph, rowMap, colMap, rowMap, rowMap, importer,
                     Teuchos::null, params));
}

}  // namespace Experimental
}  // namespace Tpetra

#endif  // TPETRA_CREATECRSGRAPHFROMSTENCIL_HPP
