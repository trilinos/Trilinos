// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_GRAPHASSEMBLY_DEF_HPP
#define TPETRA_DETAILS_GRAPHASSEMBLY_DEF_HPP

/// \file Tpetra_Details_GraphAssembly_def.hpp
///
/// Definition (implementation) of Tpetra::Details::GraphAssembly.

#include "Tpetra_Details_GraphAssembly_decl.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_Profiling.hpp"
// For makeColMapImpl, which builds a column map from the sets of local and
// remote column GIDs exactly as CrsGraph::makeColMap does.
#include "Tpetra_Details_makeColMap_def.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include "Kokkos_Sort.hpp"
#include <memory>
#include <sstream>
#include <stdexcept>

namespace Tpetra {
namespace Details {

namespace Impl {

// ---------------------------------------------------------------------------
// Index converters, used when combining the Geometry's blocks.  Each returns
// the local index in the target map, or INVALID if the input index is not in
// the source/target map.
// ---------------------------------------------------------------------------

//! Global index -> local index of a map.
template <class LocalMap>
struct GlobalToLocal {
  LocalMap map;
  template <class GO>
  KOKKOS_INLINE_FUNCTION typename LocalMap::local_ordinal_type operator()(const GO gid) const {
    return map.getLocalElement(gid);
  }
};

//! Local index of a map -> the same local index, after a range check.
template <class LO>
struct CheckLocal {
  LO numLocal;
  KOKKOS_INLINE_FUNCTION LO operator()(const LO lid) const {
    return (lid >= 0 && lid < numLocal) ? lid : ::Tpetra::Details::OrdinalTraits<LO>::invalid();
  }
};

//! Local index of one map -> local index of another map (via the global index).
template <class LocalMap>
struct LocalToLocal {
  using LO = typename LocalMap::local_ordinal_type;
  LocalMap from;
  LocalMap to;
  LO numFrom;
  KOKKOS_INLINE_FUNCTION LO operator()(const LO lid) const {
    if (lid < 0 || lid >= numFrom) return ::Tpetra::Details::OrdinalTraits<LO>::invalid();
    return to.getLocalElement(from.getGlobalElement(lid));
  }
};

//! Count the non-padding entries of each element of a block view, into
//! counts(elemOffset + e).
template <class ExecSpace, class BlockView, class Counts>
void countBlockEntries(const BlockView& block, const Counts& counts, const size_t elemOffset) {
  using value_type    = typename BlockView::non_const_value_type;
  const value_type INV = ::Tpetra::Details::OrdinalTraits<value_type>::invalid();
  const size_t width   = block.extent(1);
  Kokkos::parallel_for(
      "Tpetra::GraphAssembly::countBlockEntries",
      Kokkos::RangePolicy<ExecSpace>(0, block.extent(0)),
      KOKKOS_LAMBDA(const size_t e) {
        typename Counts::non_const_value_type count = 0;
        for (size_t n = 0; n < width; ++n)
          if (block(e, n) != INV) ++count;
        counts(elemOffset + e) = count;
      });
}

//! Write the non-padding entries of each element of a block view, converted
//! to local indices, into entries(offsets(elemOffset + e) ...).  Sets
//! errorFlag() = 1 if an index could not be converted.
template <class ExecSpace, class BlockView, class Offsets, class Entries, class Converter, class FlagView>
void fillBlockEntries(const BlockView& block, const Offsets& offsets, const Entries& entries,
                      const size_t elemOffset, const Converter& convert, const FlagView& errorFlag) {
  using value_type     = typename BlockView::non_const_value_type;
  using LO             = typename Entries::non_const_value_type;
  const value_type INV = ::Tpetra::Details::OrdinalTraits<value_type>::invalid();
  const LO LO_INV      = ::Tpetra::Details::OrdinalTraits<LO>::invalid();
  const size_t width   = block.extent(1);
  Kokkos::parallel_for(
      "Tpetra::GraphAssembly::fillBlockEntries",
      Kokkos::RangePolicy<ExecSpace>(0, block.extent(0)),
      KOKKOS_LAMBDA(const size_t e) {
        auto k = offsets(elemOffset + e);
        for (size_t n = 0; n < width; ++n) {
          const value_type index = block(e, n);
          if (index == INV) continue;
          const LO lid = convert(index);
          if (lid == LO_INV) errorFlag() = 1;
          entries(k++) = lid;
        }
      });
}

//! Exclusive prefix sum of v(0 .. n-1) into v(0 .. n); v(n) is the total.  The
//! input value v(n) is ignored.  Returns the total.
template <class ExecSpace, class Offsets>
typename Offsets::non_const_value_type exclusiveScan(const Offsets& v, const size_t n) {
  using size_type = typename Offsets::non_const_value_type;
  size_type total = 0;
  Kokkos::parallel_scan(
      "Tpetra::GraphAssembly::exclusiveScan",
      Kokkos::RangePolicy<ExecSpace>(0, n + 1),
      KOKKOS_LAMBDA(const size_t i, size_type& lsum, const bool finalPass) {
        const size_type count = (i < n) ? v(i) : size_type(0);
        if (finalPass) v(i) = lsum;
        lsum += count;
      },
      total);
  return total;
}

// ---------------------------------------------------------------------------
// Local graph assembly.  One team per row: walk the elements adjacent to the
// row, and every column of each such element, inserting each column (a local
// index of the assembly column map) into a per-team scratch hash table to
// deduplicate.  If the scratch table overflows, a global hash set of (row,
// column) pairs is used; if that overflows too, the pass is retried with a
// larger global set.
// ---------------------------------------------------------------------------

//! Insert key into the team's scratch hash table (or, on overflow, insert
//! (key, row) into the global fallback set).  Returns true iff key was not
//! already present.
template <class ScratchHashTable, class FallbackSet, class FlagView, class LO>
KOKKOS_INLINE_FUNCTION bool insertUnique(const ScratchHashTable& hashTable, const LO hashSize,
                                         const FallbackSet& fallbackSet, const FlagView& failFlag,
                                         const LO key, const size_t edge) {
  const LO EMPTY = ::Tpetra::Details::OrdinalTraits<LO>::invalid();
  size_t hash    = key;
  for (unsigned probe = 0; probe < 8; probe++) {
    hash               = KokkosKernels::Impl::xorshiftHash(hash);
    const unsigned pos = hash % hashSize;
    while (true) {
      const LO cellValue = Kokkos::volatile_load(&hashTable(pos));
      if (cellValue == key) return false;
      if (cellValue != EMPTY) break;  // occupied by another key: next probe
      if (EMPTY == Kokkos::atomic_compare_exchange(&hashTable(pos), EMPTY, key)) return true;
      // The compare-exchange lost a race; re-read this slot.
    }
  }
  // The scratch table is full along this key's probe sequence.
  auto result = fallbackSet.insert(edge);
  if (result.failed()) {
    // The global set is full: the whole pass will be retried with a bigger one.
    failFlag() = 1;
    return false;
  }
  return result.success();
}

/// COUNT pass (Fill == false): counts(row) = number of unique columns of row.
/// FILL pass (Fill == true): write the unique columns of row into
///   entries(rowptrs(row) ...), in arbitrary order.
template <bool Fill, class Offsets, class Entries, class ScratchHashTable, class ScratchCounter,
          class FallbackSet, class FlagView, class TeamMember, class LO>
struct FEAssemblyFunctor {
  using size_type = typename Offsets::non_const_value_type;

  Offsets rowOffsets;  // COUNT: output counts.  FILL: input rowptrs.
  Entries entries;     // FILL: output column indices.
  Offsets elemToColPtrs;
  Entries elemToColEntries;
  Offsets rowToElemPtrs;
  Entries rowToElemEntries;
  FallbackSet fallbackSet;
  FlagView failFlag;
  LO hashSize;
  LO numRows;

  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& t) const {
    const LO EMPTY = ::Tpetra::Details::OrdinalTraits<LO>::invalid();
    ScratchHashTable hashTable(t.team_scratch(0), hashSize);
    ScratchCounter insertCounter(t.team_scratch(0));
    Kokkos::single(Kokkos::PerTeam(t), [&]() { insertCounter() = 0; });
    Kokkos::parallel_for(Kokkos::TeamVectorRange(t, hashSize), [&](const LO i) { hashTable(i) = EMPTY; });
    t.team_barrier();

    const LO row               = t.league_rank();
    const size_type elemBegin  = rowToElemPtrs(row);
    const size_type elemEnd    = rowToElemPtrs(row + 1);
    if constexpr (Fill) {
      const size_type rowBegin = rowOffsets(row);
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(t, elemBegin, elemEnd), [&](const size_type i) {
            const LO element = rowToElemEntries(i);
            Kokkos::parallel_for(
                Kokkos::ThreadVectorRange(t, elemToColPtrs(element), elemToColPtrs(element + 1)),
                [&](const size_type j) {
                  const LO col = elemToColEntries(j);
                  if (insertUnique(hashTable, hashSize, fallbackSet, failFlag, col,
                                   size_t(col) * size_t(numRows) + size_t(row))) {
                    const LO pos               = Kokkos::atomic_fetch_add(&insertCounter(), LO(1));
                    entries(rowBegin + pos) = col;
                  }
                });
          });
    } else {
      size_type numEntries = 0;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(t, elemBegin, elemEnd),
          [&](const size_type i, size_type& lTeamCount) {
            const LO element           = rowToElemEntries(i);
            size_type numThreadEntries = 0;
            Kokkos::parallel_reduce(
                Kokkos::ThreadVectorRange(t, elemToColPtrs(element), elemToColPtrs(element + 1)),
                [&](const size_type j, size_type& lThreadCount) {
                  const LO col = elemToColEntries(j);
                  if (insertUnique(hashTable, hashSize, fallbackSet, failFlag, col,
                                   size_t(col) * size_t(numRows) + size_t(row)))
                    lThreadCount++;
                },
                numThreadEntries);
            Kokkos::single(Kokkos::PerThread(t), [&]() { lTeamCount += numThreadEntries; });
          },
          numEntries);
      Kokkos::single(Kokkos::PerTeam(t), [&]() { rowOffsets(row) = numEntries; });
    }
  }
};

/// Build a column map from a list of column GIDs (with repeats), exactly as
/// CrsGraph::makeColMap would for a graph containing those GIDs: the GIDs
/// owned by the domain map come first, in domain map order; the remaining
/// ("remote") GIDs follow, grouped by owning process, ascending within each
/// process.  Collective.
template <class LO, class GO, class NT>
Teuchos::RCP<const Map<LO, GO, NT>>
makeColMapFromGIDs(const Teuchos::RCP<const Map<LO, GO, NT>>& domainMap,
                   const Kokkos::View<const GO*, typename NT::device_type>& gids) {
  using execution_space = typename NT::execution_space;
  using memory_space    = typename NT::memory_space;
  using range_policy    = Kokkos::RangePolicy<execution_space>;

  const LO LO_INV         = ::Tpetra::Details::OrdinalTraits<LO>::invalid();
  const size_t numDomain  = domainMap->getLocalNumElements();
  auto localDomainMap     = domainMap->getLocalMap();
  const size_t numEntries = gids.extent(0);

  // Mark the domain map entries that appear, and count the remote entries.
  Kokkos::View<bool*, memory_space> isLocal("GIDisLocal", numDomain);
  size_t numRemoteEntries = 0;
  Kokkos::parallel_reduce(
      "Tpetra::GraphAssembly::markLocalColumns", range_policy(0, numEntries),
      KOKKOS_LAMBDA(const size_t k, size_t& lcount) {
        const LO lid = localDomainMap.getLocalElement(gids(k));
        if (lid == LO_INV)
          ++lcount;
        else
          isLocal(lid) = true;
      },
      numRemoteEntries);
  size_t numLocalCols = 0;
  Kokkos::parallel_reduce(
      "Tpetra::GraphAssembly::countLocalColumns", range_policy(0, numDomain),
      KOKKOS_LAMBDA(const size_t i, size_t& lcount) { lcount += isLocal(i) ? 1 : 0; },
      numLocalCols);

  // Collect the remote GIDs, sort, and remove duplicates.
  Kokkos::View<GO*, memory_space> remoteEntries(Kokkos::ViewAllocateWithoutInitializing("remoteEntries"), numRemoteEntries);
  Kokkos::parallel_scan(
      "Tpetra::GraphAssembly::collectRemoteColumns", range_policy(0, numEntries),
      KOKKOS_LAMBDA(const size_t k, size_t& lcount, const bool finalPass) {
        const GO gid = gids(k);
        if (localDomainMap.getLocalElement(gid) == LO_INV) {
          if (finalPass) remoteEntries(lcount) = gid;
          ++lcount;
        }
      });
  Kokkos::sort(execution_space(), remoteEntries);
  size_t numRemoteCols = 0;
  Kokkos::parallel_reduce(
      "Tpetra::GraphAssembly::countUniqueRemoteColumns", range_policy(0, numRemoteEntries),
      KOKKOS_LAMBDA(const size_t k, size_t& lcount) {
        if (k == 0 || remoteEntries(k) != remoteEntries(k - 1)) ++lcount;
      },
      numRemoteCols);
  Kokkos::View<GO*, memory_space> remoteCols(Kokkos::ViewAllocateWithoutInitializing("remoteCols"), numRemoteCols);
  Kokkos::parallel_scan(
      "Tpetra::GraphAssembly::uniqueRemoteColumns", range_policy(0, numRemoteEntries),
      KOKKOS_LAMBDA(const size_t k, size_t& lcount, const bool finalPass) {
        if (k == 0 || remoteEntries(k) != remoteEntries(k - 1)) {
          if (finalPass) remoteCols(lcount) = remoteEntries(k);
          ++lcount;
        }
      });

  Teuchos::RCP<const Map<LO, GO, NT>> colMap;
  Teuchos::Array<int> remotePIDs;
  std::ostringstream errStrm;
  const int lclErr = makeColMapImpl<LO, GO, NT>(colMap, remotePIDs, domainMap, numLocalCols, numRemoteCols,
                                                remoteCols, isLocal, &errStrm);
  int gblErr = 0;
  Teuchos::reduceAll<int, int>(*domainMap->getComm(), Teuchos::REDUCE_MAX, lclErr != 0 ? 1 : 0,
                               Teuchos::outArg(gblErr));
  TEUCHOS_TEST_FOR_EXCEPTION(gblErr != 0, std::runtime_error,
                             "Tpetra::Details::GraphAssembly: making the column map failed on at least one "
                             "process (are all column indices in the domain map?)."
                                 << (lclErr != 0 ? " This process: " + errStrm.str() : std::string()));
  return colMap;
}

}  // namespace Impl

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::GraphAssembly(
    const geometry_type& geometry,
    const Teuchos::RCP<const map_type>& ownedRowMap,
    const Teuchos::RCP<const map_type>& ownedPlusSharedRowMap,
    const Teuchos::RCP<const map_type>& ownedPlusSharedColMap,
    const Teuchos::RCP<const map_type>& ownedPlusSharedDomainMap,
    const Teuchos::RCP<const map_type>& ownedDomainMap,
    const Teuchos::RCP<const map_type>& ownedRangeMap,
    const Teuchos::RCP<const import_type>& ownedPlusSharedToOwnedImporter)
  : geometry_(geometry)
  , ownedRowMap_(ownedRowMap)
  , ownedPlusSharedRowMap_(ownedPlusSharedRowMap)
  , ownedPlusSharedColMap_(ownedPlusSharedColMap)
  , ownedPlusSharedDomainMap_(ownedPlusSharedDomainMap)
  , ownedDomainMap_(ownedDomainMap.is_null() ? ownedRowMap : ownedDomainMap)
  , ownedRangeMap_(ownedRangeMap.is_null() ? ownedRowMap : ownedRangeMap)
  , importer_(ownedPlusSharedToOwnedImporter) {
  TEUCHOS_TEST_FOR_EXCEPTION(ownedRowMap_.is_null(), std::invalid_argument,
                             "Tpetra::Details::GraphAssembly: ownedRowMap is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(ownedPlusSharedRowMap_.is_null(), std::invalid_argument,
                             "Tpetra::Details::GraphAssembly: ownedPlusSharedRowMap is null.");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::build() {
  using Teuchos::RCP;
  using Teuchos::rcp;

  using local_graph_type = typename crs_graph_type::local_graph_device_type;
  using rowptrs_t        = typename local_graph_type::row_map_type::non_const_type;
  using entries_t        = typename local_graph_type::entries_type::non_const_type;
  using gids_t           = Kokkos::View<global_ordinal_type*, device_type>;
  using size_type        = typename rowptrs_t::non_const_value_type;
  using LO               = local_ordinal_type;
  using GO               = global_ordinal_type;

  using range_policy       = Kokkos::RangePolicy<execution_space>;
  using team_policy        = Kokkos::TeamPolicy<execution_space>;
  using team_member        = typename team_policy::member_type;
  using scratch_space      = typename execution_space::scratch_memory_space;
  using scratch_hash_table = Kokkos::View<LO*, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using scratch_counter    = Kokkos::View<LO, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using fallback_set_type  = Kokkos::UnorderedMap<size_t, void, device_type>;
  using flag_view          = Kokkos::View<int, typename device_type::memory_space>;

  const LO LO_INV = ::Tpetra::Details::OrdinalTraits<LO>::invalid();
  auto comm       = ownedRowMap_->getComm();

  // Time each phase (enable with the environment variable TPETRA_TIMING).
  ::Tpetra::Details::ProfilingRegion buildRegion("Tpetra::GraphAssembly::build", "GraphAssembly");
  std::unique_ptr<::Tpetra::Details::ProfilingRegion> phaseRegion;
  auto phase = [&phaseRegion](const char name[]) {
    phaseRegion.reset();
    phaseRegion = std::make_unique<::Tpetra::Details::ProfilingRegion>(name, "GraphAssembly");
  };
  phase("Tpetra::GraphAssembly::checkMaps");

  // With a single graph (the owned and owned+shared row maps are the same, e.g.
  // on one process), there is no owned+shared -> owned communication.  This is
  // the same test that FECrsGraph uses.
  const bool twoGraphs = !ownedRowMap_->isSameAs(*ownedPlusSharedRowMap_);

  // The maps that FECrsGraph would use.  The owned+shared graph's domain map
  // defaults to the owned domain map (as in FECrsGraph's constructors).
  RCP<const map_type> ownedPlusSharedDomainMap =
      ownedPlusSharedDomainMap_.is_null() ? ownedDomainMap_ : ownedPlusSharedDomainMap_;
  // The assembly column map: the map whose local indices are the column
  // indices of the locally assembled graph.
  RCP<const map_type> asmColMap = !ownedPlusSharedColMap_.is_null()      ? ownedPlusSharedColMap_
                                  : !ownedPlusSharedDomainMap_.is_null() ? ownedPlusSharedDomainMap_
                                                                         : ownedPlusSharedRowMap_;

  const size_t numRows      = ownedPlusSharedRowMap_->getLocalNumElements();
  const size_t numAsmCols   = asmColMap->getLocalNumElements();
  const size_t numOwnedRows = ownedRowMap_->getLocalNumElements();
  auto localRowMap          = ownedPlusSharedRowMap_->getLocalMap();
  auto localAsmColMap       = asmColMap->getLocalMap();

  flag_view failFlag("failFlag");
  auto failFlagHost = Kokkos::create_mirror_view(failFlag);

  // The owned+shared row map must be locally fitted to the owned row map: its
  // leading entries are the owned rows, in the same order.  (The owned graph is
  // a view of the leading rows of the owned+shared graph.)
  {
    auto localOwnedRowMap = ownedRowMap_->getLocalMap();
    size_t numMismatched  = numOwnedRows > numRows ? 1 : 0;
    if (numOwnedRows <= numRows) {
      Kokkos::parallel_reduce(
          "Tpetra::GraphAssembly::checkLocallyFitted", range_policy(0, numOwnedRows),
          KOKKOS_LAMBDA(const size_t i, size_t& lcount) {
            if (localOwnedRowMap.getGlobalElement(i) != localRowMap.getGlobalElement(i)) ++lcount;
          },
          numMismatched);
    }
    int gblNotFitted = 0;
    Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MAX, numMismatched != 0 ? 1 : 0, Teuchos::outArg(gblNotFitted));
    TEUCHOS_TEST_FOR_EXCEPTION(gblNotFitted != 0, std::invalid_argument,
                               "Tpetra::Details::GraphAssembly: on at least one process, the owned+shared row map "
                               "is not locally fitted to the owned row map (its leading entries must be the owned "
                               "rows, in the same order)."
                                   << (numMismatched != 0 ? " This process' is not." : ""));
  }

  // -- Step 1: combine the blocks into element -> row and element -> column ----
  // structures (rows: elements of every block, in block order), translating
  // every index to a local index once.  Padding entries are dropped.
  phase("Tpetra::GraphAssembly::combineBlocks");
  const int numBlocks      = geometry_.getNumBlocks();
  const size_t numElements = geometry_.getTotalNumElements();
  const bool localIndices  = geometry_.hasLocalIndices();
  size_t maxColsPerElement = 0;
  for (int b = 0; b < numBlocks; ++b)
    maxColsPerElement = std::max<size_t>(maxColsPerElement, localIndices ? geometry_.getLocalCols(b).extent(1)
                                                                         : geometry_.getGlobalCols(b).extent(1));

  rowptrs_t elemToRowPtrs("elemToRowPtrs", numElements + 1);
  rowptrs_t elemToColPtrs("elemToColPtrs", numElements + 1);
  {
    size_t elemOffset = 0;
    for (int b = 0; b < numBlocks; ++b) {
      if (localIndices) {
        Impl::countBlockEntries<execution_space>(geometry_.getLocalRows(b), elemToRowPtrs, elemOffset);
        Impl::countBlockEntries<execution_space>(geometry_.getLocalCols(b), elemToColPtrs, elemOffset);
      } else {
        Impl::countBlockEntries<execution_space>(geometry_.getGlobalRows(b), elemToRowPtrs, elemOffset);
        Impl::countBlockEntries<execution_space>(geometry_.getGlobalCols(b), elemToColPtrs, elemOffset);
      }
      elemOffset += geometry_.getNumElements(b);
    }
  }
  const size_type numRowIncidences = Impl::exclusiveScan<execution_space>(elemToRowPtrs, numElements);
  const size_type numColIncidences = Impl::exclusiveScan<execution_space>(elemToColPtrs, numElements);
  entries_t elemToRowEntries(Kokkos::ViewAllocateWithoutInitializing("elemToRowEntries"), numRowIncidences);
  entries_t elemToColEntries(Kokkos::ViewAllocateWithoutInitializing("elemToColEntries"), numColIncidences);
  {
    using local_map_type = typename map_type::local_map_type;
    const Impl::GlobalToLocal<local_map_type> globalToRow{localRowMap};
    const Impl::GlobalToLocal<local_map_type> globalToAsmCol{localAsmColMap};
    const Impl::CheckLocal<LO> checkRow{static_cast<LO>(numRows)};
    const Impl::CheckLocal<LO> checkAsmCol{static_cast<LO>(numAsmCols)};
    const Impl::LocalToLocal<local_map_type> rowToAsmCol{localRowMap, localAsmColMap, static_cast<LO>(numRows)};
    const bool asmColMapIsRowMap = asmColMap.get() == ownedPlusSharedRowMap_.get();

    size_t elemOffset = 0;
    for (int b = 0; b < numBlocks; ++b) {
      if (localIndices) {
        const auto& rows = geometry_.getLocalRows(b);
        Impl::fillBlockEntries<execution_space>(rows, elemToRowPtrs, elemToRowEntries, elemOffset, checkRow, failFlag);
        if (geometry_.isRectangular(b))
          // Column indices are local indices of the assembly column map.
          Impl::fillBlockEntries<execution_space>(geometry_.getLocalCols(b), elemToColPtrs, elemToColEntries,
                                                  elemOffset, checkAsmCol, failFlag);
        else if (asmColMapIsRowMap)
          Impl::fillBlockEntries<execution_space>(rows, elemToColPtrs, elemToColEntries, elemOffset, checkRow, failFlag);
        else
          // Square block: the columns are the rows, translated to the assembly column map.
          Impl::fillBlockEntries<execution_space>(rows, elemToColPtrs, elemToColEntries, elemOffset, rowToAsmCol, failFlag);
      } else {
        Impl::fillBlockEntries<execution_space>(geometry_.getGlobalRows(b), elemToRowPtrs, elemToRowEntries,
                                                elemOffset, globalToRow, failFlag);
        Impl::fillBlockEntries<execution_space>(geometry_.getGlobalCols(b), elemToColPtrs, elemToColEntries,
                                                elemOffset, globalToAsmCol, failFlag);
      }
      elemOffset += geometry_.getNumElements(b);
    }
    Kokkos::deep_copy(failFlagHost, failFlag);
    int gblFail = 0;
    Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MAX, failFlagHost(), Teuchos::outArg(gblFail));
    TEUCHOS_TEST_FOR_EXCEPTION(gblFail != 0, std::invalid_argument,
                               "Tpetra::Details::GraphAssembly: on at least one process, an element's row index is "
                               "not in the owned+shared row map, or a column index is not in the "
                               "owned+shared column map (or, if none was given, the owned+shared domain map, "
                               "or if none was given, the owned+shared row map)."
                                   << (failFlagHost() ? " This process has such an index." : ""));
    Kokkos::deep_copy(failFlag, 0);
  }

  // -- Step 2: the row -> element structure (transpose of element -> row) -----
  phase("Tpetra::GraphAssembly::transpose");
  rowptrs_t rowToElemPtrs("rowToElemPtrs", numRows + 1);
  Kokkos::parallel_for(
      "Tpetra::GraphAssembly::countRowElements", range_policy(0, numElements),
      KOKKOS_LAMBDA(const size_t element) {
        for (size_type k = elemToRowPtrs(element); k < elemToRowPtrs(element + 1); ++k)
          Kokkos::atomic_inc(&rowToElemPtrs(elemToRowEntries(k)));
      });
  const size_type rowToElemNNZ = Impl::exclusiveScan<execution_space>(rowToElemPtrs, numRows);
  entries_t rowToElemEntries(Kokkos::ViewAllocateWithoutInitializing("rowToElemEntries"), rowToElemNNZ);
  {
    rowptrs_t insertPos(Kokkos::ViewAllocateWithoutInitializing("insertPos"), numRows + 1);
    Kokkos::deep_copy(insertPos, rowToElemPtrs);
    Kokkos::parallel_for(
        "Tpetra::GraphAssembly::fillRowElements", range_policy(0, numElements),
        KOKKOS_LAMBDA(const size_t element) {
          for (size_type k = elemToRowPtrs(element); k < elemToRowPtrs(element + 1); ++k)
            rowToElemEntries(Kokkos::atomic_fetch_add(&insertPos(elemToRowEntries(k)), size_type(1))) =
                static_cast<LO>(element);
        });
  }

  // -- Step 3: tune the team policy and scratch hash table size ---------------
  // Each team handles one row: its threads loop over the row's elements, and
  // each thread's vector lanes over an element's columns.
  const LO avgElementsPerRow = (numRows > 0) ? static_cast<LO>(rowToElemNNZ / numRows) : 1;
  int vectorLength           = 1;
  while (vectorLength < static_cast<int>(maxColsPerElement))
    vectorLength <<= 1;
  if (vectorLength > team_policy::vector_length_max())
    vectorLength = team_policy::vector_length_max();
  int teamSize     = std::max<int>(1, static_cast<int>(avgElementsPerRow));
  LO hashTableSize = 1 + avgElementsPerRow * (static_cast<LO>(maxColsPerElement) - 1);
  if (hashTableSize < 64) hashTableSize = 64;
  if (hashTableSize > 1024) hashTableSize = 1024;
  const int scratchPerTeam = scratch_hash_table::shmem_size(hashTableSize) + scratch_counter::shmem_size();

  // Run the count (Fill = false) or fill (Fill = true) pass, retrying with a
  // bigger fallback set if it overflows.
  size_t fallbackSetSize = 8 * numRows + 1;
  auto runPass           = [&](auto fillTag, const rowptrs_t& rowOffsets, const entries_t& entries) {
    constexpr bool Fill = decltype(fillTag)::value;
    using functor_type  = Impl::FEAssemblyFunctor<Fill, rowptrs_t, entries_t, scratch_hash_table, scratch_counter,
                                                 fallback_set_type, flag_view, team_member, LO>;
    while (true) {
      functor_type functor{rowOffsets, entries, elemToColPtrs, elemToColEntries, rowToElemPtrs, rowToElemEntries,
                           fallback_set_type(fallbackSetSize), failFlag, hashTableSize, static_cast<LO>(numRows)};
      team_policy dummy = team_policy(1, Kokkos::AUTO(), vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam));
      const int teamSizeRecommended = dummy.team_size_recommended(functor, Kokkos::ParallelForTag());
      if (teamSize > teamSizeRecommended) teamSize = teamSizeRecommended;
      Kokkos::parallel_for(Fill ? "Tpetra::GraphAssembly::fill" : "Tpetra::GraphAssembly::count",
                           team_policy(numRows, teamSize, vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam)),
                           functor);
      Kokkos::deep_copy(failFlagHost, failFlag);
      if (!failFlagHost()) break;
      fallbackSetSize *= 2;
      Kokkos::deep_copy(failFlag, 0);
      if (!Fill) Kokkos::deep_copy(rowOffsets, size_type(0));
    }
  };

  // -- Step 4: assemble the local owned+shared graph --------------------------
  phase("Tpetra::GraphAssembly::localAssembly");
  // Rows: ownedPlusSharedRowMap; columns: local indices of asmColMap.
  rowptrs_t asmRowptrs("asmRowptrs", numRows + 1);
  runPass(std::false_type{}, asmRowptrs, entries_t());
  const size_type asmNNZ = Impl::exclusiveScan<execution_space>(asmRowptrs, numRows);
  entries_t asmEntries(Kokkos::ViewAllocateWithoutInitializing("asmEntries"), asmNNZ);
  runPass(std::true_type{}, asmRowptrs, asmEntries);
  local_graph_type asmLocalGraph(asmEntries, asmRowptrs);
  Tpetra::Import_Util::sortCrsGraph(asmLocalGraph);

  // -- Step 5: receive other processes' contributions to the owned rows ------
  // Only the shared rows (contributions to rows owned elsewhere) need to move.
  // They form a small graph, with a compact column map of just their columns,
  // which a fused Export + fillComplete sends to the owning processes.  The
  // result has the owned row map, and holds only the received contributions.
  phase("Tpetra::GraphAssembly::exchangeSharedRows");
  const size_t numSharedRows = twoGraphs ? numRows - numOwnedRows : 0;
  rowptrs_t recvRowptrs;
  entries_t recvEntries;
  typename map_type::local_map_type localRecvColMap;
  if (twoGraphs) {
    // The shared rows' global row and column indices.
    gids_t sharedRowGIDs(Kokkos::ViewAllocateWithoutInitializing("sharedRowGIDs"), numSharedRows);
    rowptrs_t sharedRowptrs(Kokkos::ViewAllocateWithoutInitializing("sharedRowptrs"), numSharedRows + 1);
    Kokkos::parallel_for(
        "Tpetra::GraphAssembly::sharedRows", range_policy(0, numSharedRows + 1), KOKKOS_LAMBDA(const size_t i) {
          if (i < numSharedRows) sharedRowGIDs(i) = localRowMap.getGlobalElement(numOwnedRows + i);
          sharedRowptrs(i) = asmRowptrs(numOwnedRows + i) - asmRowptrs(numOwnedRows);
        });
    size_type sharedNNZ = 0;
    Kokkos::deep_copy(sharedNNZ, Kokkos::subview(sharedRowptrs, numSharedRows));
    gids_t sharedColGIDs(Kokkos::ViewAllocateWithoutInitializing("sharedColGIDs"), sharedNNZ);
    Kokkos::parallel_for(
        "Tpetra::GraphAssembly::sharedColGIDs", range_policy(0, sharedNNZ), KOKKOS_LAMBDA(const size_t k) {
          sharedColGIDs(k) = localAsmColMap.getGlobalElement(asmEntries(asmRowptrs(numOwnedRows) + k));
        });
    // The shared graph: its domain map (the owned domain map) tells the fused
    // Export which process owns each column.
    const global_size_t INV = Teuchos::OrdinalTraits<global_size_t>::invalid();
    RCP<const map_type> sharedRowMap =
        rcp(new map_type(INV, sharedRowGIDs, ownedRowMap_->getIndexBase(), comm));
    RCP<const map_type> sharedColMap = Impl::makeColMapFromGIDs<LO, GO, Node>(ownedDomainMap_, sharedColGIDs);
    auto localSharedColMap           = sharedColMap->getLocalMap();
    entries_t sharedEntries(Kokkos::ViewAllocateWithoutInitializing("sharedEntries"), sharedNNZ);
    Kokkos::parallel_for(
        "Tpetra::GraphAssembly::sharedEntries", range_policy(0, sharedNNZ), KOKKOS_LAMBDA(const size_t k) {
          sharedEntries(k) = localSharedColMap.getLocalElement(sharedColGIDs(k));
        });
    local_graph_type sharedLocalGraph(sharedEntries, sharedRowptrs);
    Tpetra::Import_Util::sortCrsGraph(sharedLocalGraph);
    RCP<const crs_graph_type> sharedGraph =
        rcp(new crs_graph_type(sharedLocalGraph, sharedRowMap, sharedColMap, ownedDomainMap_, sharedRowMap));
    const export_type exporter(sharedRowMap, ownedRowMap_);
    RCP<crs_graph_type> received =
        exportAndFillCompleteCrsGraph<crs_graph_type>(sharedGraph, exporter, ownedDomainMap_, ownedRangeMap_);
    auto recvLocalGraph = received->getLocalGraphDevice();
    recvRowptrs         = rowptrs_t(Kokkos::ViewAllocateWithoutInitializing("recvRowptrs"), recvLocalGraph.row_map.extent(0));
    Kokkos::deep_copy(recvRowptrs, recvLocalGraph.row_map);
    recvEntries     = recvLocalGraph.entries;
    localRecvColMap = received->getColMap()->getLocalMap();
  }

  // -- Step 6: the final owned+shared structure, in global indices ------------
  // Each owned row: its local contributions, followed by the received ones
  // (with possible duplicates, removed below).  Each shared row (or every row,
  // with a single graph): its local contributions.
  phase("Tpetra::GraphAssembly::finalGlobalIndices");
  const size_t numMergedRows = twoGraphs ? numOwnedRows : 0;
  rowptrs_t finalRowptrs("finalRowptrs", numRows + 1);
  Kokkos::parallel_for(
      "Tpetra::GraphAssembly::countFinalRows", range_policy(0, numRows), KOKKOS_LAMBDA(const size_t row) {
        finalRowptrs(row) = (asmRowptrs(row + 1) - asmRowptrs(row)) +
                            ((row < numMergedRows) ? recvRowptrs(row + 1) - recvRowptrs(row) : size_type(0));
      });
  const size_type finalGIDCount = Impl::exclusiveScan<execution_space>(finalRowptrs, numRows);
  gids_t finalGIDs(Kokkos::ViewAllocateWithoutInitializing("finalGIDs"), finalGIDCount);
  Kokkos::parallel_for(
      "Tpetra::GraphAssembly::finalGIDs", range_policy(0, numRows), KOKKOS_LAMBDA(const size_t row) {
        size_type out = finalRowptrs(row);
        for (size_type k = asmRowptrs(row); k < asmRowptrs(row + 1); ++k)
          finalGIDs(out++) = localAsmColMap.getGlobalElement(asmEntries(k));
        if (row < numMergedRows)
          for (size_type k = recvRowptrs(row); k < recvRowptrs(row + 1); ++k)
            finalGIDs(out++) = localRecvColMap.getGlobalElement(recvEntries(k));
      });

  // -- Step 7: the final column map ------------------------------------------
  phase("Tpetra::GraphAssembly::makeColMap");
  // As in FECrsGraph: the given column map, or else one made from all the
  // owned+shared graph's columns (owned rows merged, shared rows local), with
  // the owned+shared domain map (or, with a single graph, the owned domain map).
  RCP<const map_type> finalColMap = ownedPlusSharedColMap_;
  if (finalColMap.is_null())
    finalColMap = Impl::makeColMapFromGIDs<LO, GO, Node>(twoGraphs ? ownedPlusSharedDomainMap : ownedDomainMap_,
                                                         finalGIDs);

  // -- Step 8: translate to local indices of the final column map, and sort ----
  phase("Tpetra::GraphAssembly::finalLocalIndices");
  // With a given column map, columns not in it are dropped (as FECrsGraph does
  // for a locally indexed graph).
  auto localFinalColMap = finalColMap->getLocalMap();
  rowptrs_t lclRowptrs("lclRowptrs", numRows + 1);
  Kokkos::parallel_for(
      "Tpetra::GraphAssembly::countFinalEntries", range_policy(0, numRows), KOKKOS_LAMBDA(const size_t row) {
        size_type count = 0;
        for (size_type k = finalRowptrs(row); k < finalRowptrs(row + 1); ++k)
          if (localFinalColMap.getLocalElement(finalGIDs(k)) != LO_INV) ++count;
        lclRowptrs(row) = count;
      });
  const size_type lclNNZ = Impl::exclusiveScan<execution_space>(lclRowptrs, numRows);
  entries_t lclEntries(Kokkos::ViewAllocateWithoutInitializing("lclEntries"), lclNNZ);
  Kokkos::parallel_for(
      "Tpetra::GraphAssembly::fillFinalEntries", range_policy(0, numRows), KOKKOS_LAMBDA(const size_t row) {
        size_type out = lclRowptrs(row);
        for (size_type k = finalRowptrs(row); k < finalRowptrs(row + 1); ++k) {
          const LO lid = localFinalColMap.getLocalElement(finalGIDs(k));
          if (lid != LO_INV) lclEntries(out++) = lid;
        }
      });
  // Sort each row, and remove the duplicates between local and received
  // contributions to the owned rows.
  rowptrs_t mergedRowptrs;
  entries_t mergedEntries;
  KokkosSparse::sort_and_merge_graph(execution_space(), lclRowptrs, lclEntries, mergedRowptrs, mergedEntries,
                                     static_cast<LO>(finalColMap->getLocalNumElements()));
  local_graph_type finalLocalGraph(mergedEntries, mergedRowptrs);

  // The owned+shared graph, with the domain and range maps FECrsGraph uses.
  phase("Tpetra::GraphAssembly::makeGraph");
  if (twoGraphs)
    ownedPlusSharedGraph_ = rcp(new crs_graph_type(finalLocalGraph, ownedPlusSharedRowMap_, finalColMap,
                                                   ownedPlusSharedDomainMap, ownedPlusSharedRowMap_));
  else
    ownedPlusSharedGraph_ = rcp(new crs_graph_type(finalLocalGraph, ownedPlusSharedRowMap_, finalColMap,
                                                   ownedDomainMap_, ownedRangeMap_));
  graph_ = Teuchos::null;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::crs_graph_type>
GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::getOwnedPlusSharedGraph() const {
  return ownedPlusSharedGraph_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::crs_graph_type>
GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::getGraph() const {
  if (graph_.is_null() && !ownedPlusSharedGraph_.is_null()) {
    if (ownedRowMap_->isSameAs(*ownedPlusSharedRowMap_)) {
      graph_ = ownedPlusSharedGraph_;
    } else {
      // A view of the leading (owned) rows, as FECrsGraph makes it.
      graph_ = Teuchos::rcp(new crs_graph_type(*ownedPlusSharedGraph_, ownedRowMap_));
      graph_->fillComplete(ownedDomainMap_, ownedRangeMap_);
    }
  }
  return graph_;
}

}  // namespace Details
}  // namespace Tpetra

//
// Explicit template instantiation macro for
// Tpetra::Details::GraphAssembly.  NOT FOR USERS!!!  Must be used
// inside the Tpetra namespace.
//
#define TPETRA_DETAILS_GRAPHASSEMBLY_INSTANT(LO, GO, NODE) \
  template class Details::GraphAssembly<LO, GO, NODE>;

#endif  // TPETRA_DETAILS_GRAPHASSEMBLY_DEF_HPP
