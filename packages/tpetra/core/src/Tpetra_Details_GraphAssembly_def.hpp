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
#include "Tpetra_Details_makeColMap.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include <type_traits>
#include <set>

namespace Tpetra {
namespace Details {

namespace Impl {

// ---------------------------------------------------------------------------
// Owned-only graph-assembly functors.  These build the local sparsity pattern
// of the FE operator for the OWNED+SHARED rows.  An owned element only touches
// owned or shared nodes, so the owned+shared map serves as both the row map and
// the column map of the locally-assembled graph, and its LIDs are the column
// indices written out.
//
// These are adapted from the prototype in the FE-assembly example
// (fem_assembly_InsertGlobalIndices_FE.hpp), generalized over the CrsGraph
// template parameters.
// ---------------------------------------------------------------------------

/// COUNT pass: one team per row (owned+shared node).  Walk the adjacent owned
/// elements and count the unique neighbor node LIDs, using a per-team scratch
/// linear-probing hash table with a global fallback set for overflow.  The hash
/// table stores LOCAL column IDs (LIDs) rather than global IDs, since LIDs are
/// typically a narrower type -- this conserves scratch (shared) memory.
template <class ElementToNode, class RowptrsType, class EntriesType,
          class LocalMapType, class ScratchHashTable, class GlobalEdgeSet,
          class FlagView, class TeamMember, class LocalOrdinal, class GlobalOrdinal>
struct FECountEntriesFunctor {
  using size_type = typename RowptrsType::value_type;

  FECountEntriesFunctor(
      const RowptrsType& counts_,
      const ElementToNode& ownedElementToNode_,
      const RowptrsType& nodeToElementRowptrs_,
      const EntriesType& nodeToElementEntries_,
      const LocalMapType& columnMap_,
      GlobalEdgeSet fallbackEdgeSet_,
      const FlagView& globalFailFlag_,
      LocalOrdinal hashSize_)
    : counts(counts_)
    , ownedElementToNode(ownedElementToNode_)
    , nodeToElementRowptrs(nodeToElementRowptrs_)
    , nodeToElementEntries(nodeToElementEntries_)
    , columnMap(columnMap_)
    , fallbackEdgeSet(fallbackEdgeSet_)
    , globalFailFlag(globalFailFlag_)
    , hashSize(hashSize_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& t) const {
    const LocalOrdinal LO_INVALID = ~LocalOrdinal(0);
    ScratchHashTable hashTable(t.team_scratch(0), hashSize);
    // Initially, fill the table with INVALID to mean 'empty'.
    Kokkos::parallel_for(Kokkos::TeamVectorRange(t, hashSize),
                         [&](int i) { hashTable(i) = LO_INVALID; });
    t.team_barrier();
    LocalOrdinal localRow            = t.league_rank();
    const LocalOrdinal numLocalNodes = counts.extent(0) - 1;
    size_type elementBegin           = nodeToElementRowptrs(localRow);
    size_type elementEnd             = nodeToElementRowptrs(localRow + 1);
    size_type numEntries;
    // Iterate over the elements adjacent to this node.
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(t, elementBegin, elementEnd),
        [&](size_type i, size_type& lTeamCount) {
          LocalOrdinal localElement = nodeToElementEntries(i);
          size_type numThreadEntries;
          // Iterate over the nodes adjacent to localElement.
          Kokkos::parallel_reduce(
              Kokkos::ThreadVectorRange(t, ownedElementToNode.extent(1)),
              [&](int j, size_type& lThreadCount) {
                GlobalOrdinal nei = ownedElementToNode(localElement, j);
                // Work with the LOCAL column ID; the hash table stores LIDs.
                LocalOrdinal neiLid = columnMap.getLocalElement(nei);
                // Try to insert neiLid into the scratch hash table, if it's not
                // already there.
                size_t hash          = neiLid;
                bool foundOrInserted = false;
                for (unsigned probe = 0; probe < 8; probe++) {
                  hash         = KokkosKernels::Impl::xorshiftHash(hash);
                  unsigned pos = hash % hashSize;
                  while (true) {
                    LocalOrdinal cellValue = Kokkos::volatile_load(&hashTable(pos));
                    if (cellValue == neiLid) {
                      foundOrInserted = true;
                      break;
                    } else if (cellValue == LO_INVALID) {
                      if (LO_INVALID == Kokkos::atomic_compare_exchange(&hashTable(pos), LO_INVALID, neiLid)) {
                        foundOrInserted = true;
                        lThreadCount++;
                        break;
                      }
                      // The cmp-xchg failed, so we don't know yet what another
                      // thread put at pos.  Retry this slot.
                    } else {
                      // This cell holds an entry other than neiLid; next probe.
                      break;
                    }
                  }
                  if (foundOrInserted)
                    break;
                }
                if (!foundOrInserted) {
                  // Scratch hash table ran out of space: use the global set.
                  size_t edge       = size_t(neiLid) * size_t(numLocalNodes) + localRow;
                  auto insertResult = fallbackEdgeSet.insert(edge);
                  if (insertResult.failed()) {
                    // The global table filled up: grow it and start over.
                    globalFailFlag() = 1;
                  }
                  if (!insertResult.existing()) {
                    lThreadCount++;
                  }
                }
              },
              numThreadEntries);
          Kokkos::single(Kokkos::PerThread(t), [&]() { lTeamCount += numThreadEntries; });
        },
        numEntries);
    Kokkos::single(Kokkos::PerTeam(t), [&]() { counts(localRow) = numEntries; });
  }

  RowptrsType counts;
  ElementToNode ownedElementToNode;
  RowptrsType nodeToElementRowptrs;
  EntriesType nodeToElementEntries;
  LocalMapType columnMap;
  GlobalEdgeSet fallbackEdgeSet;
  FlagView globalFailFlag;
  LocalOrdinal hashSize;
};

/// FILL pass: identical traversal to the count pass, but each time a new unique
/// neighbor is discovered a column index is written into the entries array.
///
/// The scratch hash table always keys on the owned+shared LID (obtained from
/// keyMap), which is a compact, dedup-correct identifier because every neighbor
/// node of an owned element is in the owned+shared map.  What gets WRITTEN into
/// the entries array depends on StoreGlobal:
///   - StoreGlobal == false: the local column index obtained from the (provided)
///     output column map outColMap is written (path 1: colMap supplied).
///   - StoreGlobal == true:  the global node ID is written directly, producing a
///     globally-indexed local graph (path 2: no colMap; makeColMap later).
template <bool StoreGlobal, class ElementToNode, class RowptrsType,
          class EntriesType, class NodeToElemEntriesType, class KeyMapType,
          class OutColMapType, class ScratchHashTable, class ScratchCounter,
          class GlobalEdgeSet, class FlagView, class TeamMember,
          class LocalOrdinal, class GlobalOrdinal>
struct FEFillEntriesFunctor {
  using size_type    = typename RowptrsType::value_type;
  using out_value_type = typename EntriesType::non_const_value_type;

  FEFillEntriesFunctor(
      const RowptrsType& rowptrs_,
      const EntriesType& entries_,
      const ElementToNode& ownedElementToNode_,
      const RowptrsType& nodeToElementRowptrs_,
      const NodeToElemEntriesType& nodeToElementEntries_,
      const KeyMapType& keyMap_,
      const OutColMapType& outColMap_,
      GlobalEdgeSet fallbackEdgeSet_,
      const FlagView& globalFailFlag_,
      LocalOrdinal hashSize_)
    : rowptrs(rowptrs_)
    , entries(entries_)
    , ownedElementToNode(ownedElementToNode_)
    , nodeToElementRowptrs(nodeToElementRowptrs_)
    , nodeToElementEntries(nodeToElementEntries_)
    , keyMap(keyMap_)
    , outColMap(outColMap_)
    , fallbackEdgeSet(fallbackEdgeSet_)
    , globalFailFlag(globalFailFlag_)
    , hashSize(hashSize_) {}

  //! Compute the value written into the entries array for global node ID nei.
  KOKKOS_INLINE_FUNCTION out_value_type
  outputValue(const GlobalOrdinal nei, const LocalOrdinal neiKey) const {
    if (StoreGlobal)
      return static_cast<out_value_type>(nei);
    else
      return static_cast<out_value_type>(outColMap.getLocalElement(nei));
  }

  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& t) const {
    const LocalOrdinal LO_INVALID = ~LocalOrdinal(0);
    ScratchHashTable hashTable(t.team_scratch(0), hashSize);
    // Use this rank-0 view to atomically count up positions in the row as
    // entries are inserted.
    ScratchCounter insertCounter(t.team_scratch(0));
    Kokkos::single(Kokkos::PerTeam(t), [&]() { insertCounter() = 0; });
    Kokkos::parallel_for(Kokkos::TeamVectorRange(t, hashSize),
                         [&](int i) { hashTable(i) = LO_INVALID; });
    t.team_barrier();
    const LocalOrdinal numLocalNodes = rowptrs.extent(0) - 1;
    LocalOrdinal localRow            = t.league_rank();
    size_type elementBegin           = nodeToElementRowptrs(localRow);
    size_type elementEnd             = nodeToElementRowptrs(localRow + 1);
    // Iterate over the elements adjacent to this node.
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(t, elementBegin, elementEnd),
        [&](size_type i) {
          LocalOrdinal localElement = nodeToElementEntries(i);
          // Iterate over the nodes adjacent to localElement.
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(t, ownedElementToNode.extent(1)),
              [&](int j) {
                GlobalOrdinal nei = ownedElementToNode(localElement, j);
                // The hash table keys on the owned+shared LID (compact and
                // dedup-correct); the value written to entries is computed by
                // outputValue (either an output-colMap LID or the GID).
                LocalOrdinal neiKey  = keyMap.getLocalElement(nei);
                size_t hash          = neiKey;
                bool foundOrInserted = false;
                for (unsigned probe = 0; probe < 8; probe++) {
                  hash         = KokkosKernels::Impl::xorshiftHash(hash);
                  unsigned pos = (hash + probe) % hashSize;
                  while (true) {
                    LocalOrdinal cellValue = Kokkos::volatile_load(&hashTable(pos));
                    if (cellValue == neiKey) {
                      foundOrInserted = true;
                      break;
                    } else if (cellValue == LO_INVALID) {
                      if (LO_INVALID == Kokkos::atomic_compare_exchange(&hashTable(pos), LO_INVALID, neiKey)) {
                        foundOrInserted                        = true;
                        LocalOrdinal insertPos                 = Kokkos::atomic_fetch_add(&insertCounter(), 1);
                        entries(rowptrs(localRow) + insertPos) = outputValue(nei, neiKey);
                        break;
                      }
                      // cmp-xchg failed; retry this slot.
                    } else {
                      break;
                    }
                  }
                  if (foundOrInserted)
                    break;
                }
                if (!foundOrInserted) {
                  size_t edge       = size_t(neiKey) * size_t(numLocalNodes) + localRow;
                  auto insertResult = fallbackEdgeSet.insert(edge);
                  if (insertResult.failed()) {
                    globalFailFlag() = 1;
                  }
                  if (!insertResult.existing()) {
                    LocalOrdinal insertPos                 = Kokkos::atomic_fetch_add(&insertCounter(), 1);
                    entries(rowptrs(localRow) + insertPos) = outputValue(nei, neiKey);
                  }
                }
              });
        });
  }

  RowptrsType rowptrs;
  EntriesType entries;
  ElementToNode ownedElementToNode;
  RowptrsType nodeToElementRowptrs;
  NodeToElemEntriesType nodeToElementEntries;
  KeyMapType keyMap;
  OutColMapType outColMap;
  GlobalEdgeSet fallbackEdgeSet;
  FlagView globalFailFlag;
  LocalOrdinal hashSize;
};

}  // namespace Impl

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::GraphAssembly(
    const Teuchos::RCP<const map_type>& rowMap,
    const Teuchos::RCP<const map_type>& ownedPlusSharedMap,
    const element_to_node_type& ownedElementToNode,
    const Teuchos::RCP<const map_type>& ownedPlusSharedColMap)
  : rowMap_(rowMap)
  , ownedPlusSharedMap_(ownedPlusSharedMap)
  , ownedElementToNode_(ownedElementToNode)
  , ownedPlusSharedColMap_(ownedPlusSharedColMap) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::build() {
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Local types used only by the kernels below.
  using local_graph_type = typename crs_graph_type::local_graph_device_type;
  using rowptrs_t        = typename local_graph_type::row_map_type::non_const_type;
  using entries_t        = typename local_graph_type::entries_type;
  using size_type        = typename rowptrs_t::value_type;
  using local_map_type   = typename map_type::local_map_type;

  using range_policy  = Kokkos::RangePolicy<execution_space>;
  using team_policy   = Kokkos::TeamPolicy<execution_space>;
  using team_member   = typename team_policy::member_type;
  using scratch_space = typename execution_space::scratch_memory_space;
  using scratch_hash_table =
      Kokkos::View<local_ordinal_type*, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using scratch_counter =
      Kokkos::View<local_ordinal_type, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using global_edge_set = Kokkos::UnorderedMap<size_t, void, device_type>;
  using flag_view       = Kokkos::View<int, typename device_type::memory_space>;

  const local_ordinal_type LO_INVALID = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();

  // Whether the user provided an owned+shared column map.  If so, we take the
  // "path 1" codepath and map global column IDs to local column indices inside
  // the fill functor.  Otherwise we take "path 2": build a globally-indexed
  // local graph, call makeColMap, and remap to local indices.
  const bool haveColMap = !ownedPlusSharedColMap_.is_null();

  auto owned_element_to_node_ids = ownedElementToNode_;

  // The hash table always keys on the owned+shared map's LIDs (compact and
  // dedup-correct): an owned element's nodes are, by construction, all either
  // owned or shared on this process.
  auto localOwnedPlusSharedMap = ownedPlusSharedMap_->getLocalMap();
  auto numLocalNodes           = ownedPlusSharedMap_->getLocalNumElements();

  flag_view failFlag("failFlag");
  auto failFlagHost = Kokkos::create_mirror_view(failFlag);

  // -- Step 1: build the node -> element graph (transpose of element -> node) --
  // Rows are the local (owned+shared) nodes; entries are owned element indices.
  rowptrs_t nodeToElementRowptrs("nodeToElementRowptrs", numLocalNodes + 1);
  Kokkos::parallel_for(
      range_policy(0, owned_element_to_node_ids.size()),
      KOKKOS_LAMBDA(size_t i) {
        const local_ordinal_type nodesPerElement = owned_element_to_node_ids.extent(1);
        local_ordinal_type ownedElementIndex     = i / nodesPerElement;
        local_ordinal_type nodeOfElem            = i % nodesPerElement;
        global_ordinal_type globalNode           = owned_element_to_node_ids(ownedElementIndex, nodeOfElem);
        local_ordinal_type localNode             = localOwnedPlusSharedMap.getLocalElement(globalNode);
        if (localNode != LO_INVALID)
          Kokkos::atomic_inc(&nodeToElementRowptrs(localNode));
      });
  typename rowptrs_t::value_type nodeToElementNNZ;
  Kokkos::parallel_scan(
      range_policy(0, numLocalNodes + 1),
      KOKKOS_LAMBDA(local_ordinal_type i, size_type & lsum, bool finalPass) {
        size_type count = nodeToElementRowptrs(i);
        if (finalPass)
          nodeToElementRowptrs(i) = lsum;
        lsum += count;
      },
      nodeToElementNNZ);
  entries_t nodeToElementEntries(Kokkos::ViewAllocateWithoutInitializing("nodeToElementEntries"), nodeToElementNNZ);
  {
    rowptrs_t insertPos(Kokkos::ViewAllocateWithoutInitializing("insertPos"), numLocalNodes + 1);
    Kokkos::deep_copy(insertPos, nodeToElementRowptrs);
    Kokkos::parallel_for(
        range_policy(0, owned_element_to_node_ids.size()),
        KOKKOS_LAMBDA(size_t i) {
          const local_ordinal_type nodesPerElement = owned_element_to_node_ids.extent(1);
          local_ordinal_type ownedElementIndex     = i / nodesPerElement;
          local_ordinal_type nodeOfElem            = i % nodesPerElement;
          global_ordinal_type globalNode           = owned_element_to_node_ids(ownedElementIndex, nodeOfElem);
          local_ordinal_type localNode             = localOwnedPlusSharedMap.getLocalElement(globalNode);
          if (localNode != LO_INVALID)
            nodeToElementEntries(Kokkos::atomic_fetch_add(&insertPos(localNode), size_type(1))) = ownedElementIndex;
        });
  }

  // -- Step 2: tune the team policy and scratch hash table size --
  local_ordinal_type avgElementsPerNode = (numLocalNodes > 0) ? (nodeToElementNNZ / numLocalNodes) : 1;
  int vectorLength                      = 1;
  while (vectorLength < avgElementsPerNode)
    vectorLength <<= 1;
  if (vectorLength > team_policy::vector_length_max())
    vectorLength = team_policy::vector_length_max();
  const local_ordinal_type nodesPerElement = owned_element_to_node_ids.extent(1);
  int teamSize                             = nodesPerElement;
  local_ordinal_type hashTableSize         = 1 + avgElementsPerNode * (nodesPerElement - 1);
  if (hashTableSize < 64)
    hashTableSize = 64;
  if (hashTableSize > 1024)
    hashTableSize = 1024;

  // -- Step 3: COUNT pass: unique entries per row --
  rowptrs_t localRowptrs("localRowptrs", numLocalNodes + 1);
  size_t fallbackSetSize = 8 * numLocalNodes + 1;
  {
    using count_functor_type =
        Impl::FECountEntriesFunctor<element_to_node_type, rowptrs_t, entries_t,
                                    local_map_type, scratch_hash_table, global_edge_set,
                                    flag_view, team_member, local_ordinal_type, global_ordinal_type>;
    while (true) {
      global_edge_set fallbackSet(fallbackSetSize);
      count_functor_type functor(
          localRowptrs, owned_element_to_node_ids,
          nodeToElementRowptrs, nodeToElementEntries,
          localOwnedPlusSharedMap, fallbackSet, failFlag, hashTableSize);
      int scratchPerTeam      = scratch_hash_table::shmem_size(hashTableSize);
      team_policy dummy       = team_policy(1, Kokkos::AUTO(), vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam));
      int teamSizeRecommended = dummy.team_size_recommended(functor, Kokkos::ParallelForTag());
      if (teamSize > teamSizeRecommended)
        teamSize = teamSizeRecommended;
      Kokkos::parallel_for(
          team_policy(numLocalNodes, teamSize, vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam)),
          functor);
      Kokkos::deep_copy(failFlagHost, failFlag);
      if (failFlagHost()) {
        // Need to retry the whole pass with a bigger fallback set.
        fallbackSetSize *= 2;
        Kokkos::deep_copy(failFlag, 0);
        Kokkos::deep_copy(localRowptrs, 0);
      } else
        break;
    }
  }

  // -- Step 4: SCAN: counts -> row offsets --
  typename rowptrs_t::value_type nnz;
  Kokkos::parallel_scan(
      range_policy(0, numLocalNodes + 1),
      KOKKOS_LAMBDA(local_ordinal_type i, size_type & lsum, bool finalPass) {
        size_type count = localRowptrs(i);
        if (finalPass)
          localRowptrs(i) = lsum;
        lsum += count;
      },
      nnz);

  // The (possibly reused) counter type for the fill functor.
  auto runFill = [&](auto entriesView, auto outColLocalMap, auto storeGlobalTag) {
    using entries_view_type = decltype(entriesView);
    using out_col_map_type  = decltype(outColLocalMap);
    constexpr bool StoreGlobal = decltype(storeGlobalTag)::value;
    using fill_functor_type =
        Impl::FEFillEntriesFunctor<StoreGlobal, element_to_node_type, rowptrs_t,
                                   entries_view_type, entries_t, local_map_type, out_col_map_type,
                                   scratch_hash_table, scratch_counter,
                                   global_edge_set, flag_view, team_member,
                                   local_ordinal_type, global_ordinal_type>;
    size_t localFallbackSetSize = fallbackSetSize;
    while (true) {
      global_edge_set fallbackSet(localFallbackSetSize);
      fill_functor_type functor(
          localRowptrs, entriesView, owned_element_to_node_ids,
          nodeToElementRowptrs, nodeToElementEntries,
          localOwnedPlusSharedMap, outColLocalMap, fallbackSet, failFlag, hashTableSize);
      int scratchPerTeam      = scratch_hash_table::shmem_size(hashTableSize) + scratch_counter::shmem_size();
      team_policy dummy       = team_policy(1, Kokkos::AUTO(), vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam));
      int teamSizeRecommended = dummy.team_size_recommended(functor, Kokkos::ParallelForTag());
      if (teamSize > teamSizeRecommended)
        teamSize = teamSizeRecommended;
      Kokkos::parallel_for(
          team_policy(numLocalNodes, teamSize, vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam)),
          functor);
      Kokkos::deep_copy(failFlagHost, failFlag);
      if (failFlagHost()) {
        localFallbackSetSize *= 2;
        Kokkos::deep_copy(failFlag, 0);
      } else
        break;
    }
  };

  // The owned+shared column map used by the assembled owned+shared graph, and
  // the local (colMap-indexed) entries produced by the fill pass.  These are
  // populated by one of the two paths below.
  RCP<const map_type> ownedPlusSharedColMap;
  entries_t localEntries(Kokkos::ViewAllocateWithoutInitializing("localEntries"), nnz);

  if (haveColMap) {
    // -- Path 1: the user provided a column map. --------------------------
    // FILL directly into a locally-indexed entries array, mapping global node
    // IDs to local column indices via the provided column map inside the
    // functor.
    ownedPlusSharedColMap = ownedPlusSharedColMap_;
    auto localColMap      = ownedPlusSharedColMap->getLocalMap();
    runFill(localEntries, localColMap, std::false_type{});
  } else {
    // -- Path 2: no column map provided. ----------------------------------
    // The assembled owned+shared graph's rows and columns are both indexed by
    // the owned+shared map: an owned element's nodes are, by construction, all
    // either owned or shared on this process, so the owned+shared map serves as
    // the column map too and no makeColMap call is needed at this stage.  The
    // fill pass therefore writes owned+shared LIDs directly.
    ownedPlusSharedColMap = ownedPlusSharedMap_;
    runFill(localEntries, localOwnedPlusSharedMap, std::false_type{});
  }

  ownedPlusSharedColMap_ = ownedPlusSharedColMap;

  // Wrap the assembled (local, sorted) structure in a fill-complete owned+shared
  // graph.  The 5-argument (lclGraph, row, col, domain, range) constructor sets
  // domain = owned and range = owned+shared without any column reindexing, so
  // the local column indices we just computed stay valid.  The owned domain map
  // is what makes CrsGraph build the owned -> owned+shared Importer that the
  // fused Export below relies on.  Its rows for shared nodes hold this process'
  // partial contributions to rows owned elsewhere.
  local_graph_type localGraph(localEntries, localRowptrs);
  Tpetra::Import_Util::sortCrsGraph(localGraph);
  ownedPlusSharedGraph_ =
      rcp(new crs_graph_type(localGraph,
                             ownedPlusSharedMap_,   // row map
                             ownedPlusSharedColMap, // column map
                             rowMap_,               // domain map (owned)
                             ownedPlusSharedMap_)); // range map

  // -- Final step: build the owned graph and finalize the owned+shared graph. --
  //
  // For compatibility with FECrsMatrix (whose owned matrix aliases the first
  // chunk of the owned+shared matrix's values), the owned graph must be a
  // first-chunk subview of the owned+shared graph, sharing its column map.  For
  // that to be correct, the owned rows of the owned+shared graph must contain
  // the fully MERGED structure (the union of every process' contributions),
  // while the shared rows keep only this process' partial contributions (so the
  // matrix Export with ADD does not double count).  This mirrors exactly what
  // FECrsGraph::doOwnedPlusSharedToOwned does: a restricted-mode self-export
  // (which merges the owned rows in place but leaves the shared rows untouched),
  // followed by makeColMap over all rows and an owned first-chunk subview.
  if (ownedPlusSharedMap_->isSameAs(*rowMap_)) {
    // Serial / single-rank case: nothing to communicate.  The assembly graph
    // already has exactly the owned structure; just re-map domain/range.
    ownedPlusSharedGraph_->resumeFill();
    ownedPlusSharedGraph_->fillComplete(rowMap_, rowMap_);
    graph_    = ownedPlusSharedGraph_;
    exporter_ = Teuchos::null;
  } else {
    RCP<const export_type> exporter = rcp(new export_type(ownedPlusSharedMap_, rowMap_));
    exporter_                       = exporter;

    // (a) Build the merged owned graph via a fused Export + fillComplete.  This
    //     gives the complete owned rows (union of all processes' contributions).
    //     Its column map is built by the fused path
    //     (lowCommunicationMakeColMapAndReindex) with the OWNED domain map
    //     (rowMap_), so the owned GIDs form the leading chunk followed by remote
    //     GIDs grouped by owning process -- identical to what the serial
    //     FECrsGraph::doOwnedPlusSharedToOwned path produces via
    //     CrsGraph::makeColMap(getDomainMap()).  This is the graph returned by
    //     getGraph(), and it defines the owned matrix's column map, keeping the
    //     assembled owned matrix bit-for-bit identical to the serial path.
    graph_ = Tpetra::exportAndFillCompleteCrsGraph<crs_graph_type>(
        ownedPlusSharedGraph_, *exporter, rowMap_, rowMap_);

    // (b) Finalize the owned+shared graph so that it shares graph_'s column map
    //     and so that its owned rows carry the merged structure, exactly as the
    //     serial owned+shared graph does after the restricted-mode self-export.
    //     We assemble its local structure directly, without any further
    //     communication:
    //       - owned rows [0, numOwnedRows)          <- graph_'s merged rows,
    //       - shared rows [numOwnedRows, numOARows)  <- this process' partial
    //                                                   contributions (the
    //                                                   original assembled graph,
    //                                                   reindexed to graph_'s
    //                                                   column map).
    //     The owned map is the leading chunk of the owned+shared map (locally
    //     fitted), so the owned rows are precisely the first numOwnedRows rows.
    //     graph_'s column map is a superset of the columns used by every
    //     owned+shared row (the shared-row columns are all owned+shared nodes,
    //     which the fused Export's colMap retains), so the reindex is exact.
    {
      RCP<const map_type> finalColMap = graph_->getColMap();
      auto finalLocalColMap           = finalColMap->getLocalMap();

      const local_ordinal_type numOwnedRows =
          static_cast<local_ordinal_type>(rowMap_->getLocalNumElements());
      const local_ordinal_type numOARows =
          static_cast<local_ordinal_type>(ownedPlusSharedMap_->getLocalNumElements());

      // The merged owned rows (already indexed against finalColMap).
      auto ownedLclGraph          = graph_->getLocalGraphDevice();
      auto ownedRowptrs           = ownedLclGraph.row_map;
      auto ownedEntries           = ownedLclGraph.entries;

      // The original assembled owned+shared structure (indexed against the
      // assembly column map, i.e. the owned+shared map for path 2, or the
      // user-provided column map for path 1).  Reindex its shared rows to
      // finalColMap via global IDs.
      auto asmLclGraph            = ownedPlusSharedGraph_->getLocalGraphDevice();
      auto asmRowptrs             = asmLclGraph.row_map;
      auto asmEntries             = asmLclGraph.entries;
      auto asmLocalColMap         = ownedPlusSharedColMap->getLocalMap();

      // Build the final rowptrs: owned rows take their length from graph_, the
      // shared rows keep their assembled length.
      rowptrs_t finalRowptrs("finalRowptrs", numOARows + 1);
      Kokkos::parallel_for(
          "finalOwnedPlusSharedRowCounts",
          range_policy(0, numOARows),
          KOKKOS_LAMBDA(local_ordinal_type row) {
            if (row < numOwnedRows)
              finalRowptrs(row) = ownedRowptrs(row + 1) - ownedRowptrs(row);
            else
              finalRowptrs(row) = asmRowptrs(row + 1) - asmRowptrs(row);
          });
      size_type finalNNZ;
      Kokkos::parallel_scan(
          range_policy(0, numOARows + 1),
          KOKKOS_LAMBDA(local_ordinal_type i, size_type & lsum, bool finalPass) {
            size_type count = finalRowptrs(i);
            if (finalPass)
              finalRowptrs(i) = lsum;
            lsum += count;
          },
          finalNNZ);

      entries_t finalEntries(Kokkos::ViewAllocateWithoutInitializing("finalEntries"), finalNNZ);
      Kokkos::parallel_for(
          "finalOwnedPlusSharedFill",
          range_policy(0, numOARows),
          KOKKOS_LAMBDA(local_ordinal_type row) {
            const size_type out = finalRowptrs(row);
            if (row < numOwnedRows) {
              // Copy the merged owned row verbatim (already finalColMap LIDs).
              const size_type beg = ownedRowptrs(row);
              const size_type end = ownedRowptrs(row + 1);
              for (size_type k = beg; k < end; ++k)
                finalEntries(out + (k - beg)) = ownedEntries(k);
            } else {
              // Reindex the assembled shared row into finalColMap via GID.
              const size_type beg = asmRowptrs(row);
              const size_type end = asmRowptrs(row + 1);
              for (size_type k = beg; k < end; ++k) {
                const global_ordinal_type gid =
                    asmLocalColMap.getGlobalElement(asmEntries(k));
                finalEntries(out + (k - beg)) =
                    finalLocalColMap.getLocalElement(gid);
              }
            }
          });

      local_graph_type finalLocalGraph(finalEntries, finalRowptrs);
      // Owned rows are already sorted (copied from graph_); the shared rows were
      // sorted in the owned+shared column map, but finalColMap may order columns
      // differently, so sort again to be safe.
      Tpetra::Import_Util::sortCrsGraph(finalLocalGraph);

      ownedPlusSharedGraph_ =
          rcp(new crs_graph_type(finalLocalGraph,
                                 ownedPlusSharedMap_,  // row map
                                 finalColMap,          // column map (== graph_'s)
                                 rowMap_,              // domain map (owned)
                                 ownedPlusSharedMap_));// range map
    }
    ownedPlusSharedColMap_ = ownedPlusSharedGraph_->getColMap();
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::crs_graph_type>
GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::getGraph() const {
  return graph_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::crs_graph_type>
GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::getOwnedPlusSharedGraph() const {
  return ownedPlusSharedGraph_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::export_type>
GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::getExporter() const {
  return exporter_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::map_type>
GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::getOwnedPlusSharedColMap() const {
  return ownedPlusSharedColMap_;
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
