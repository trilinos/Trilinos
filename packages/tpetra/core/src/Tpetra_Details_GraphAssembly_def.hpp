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
  // Globally-indexed entries array (path 2): holds GlobalOrdinal column IDs.
  using gbl_entries_t =
      Kokkos::View<global_ordinal_type*, typename node_type::memory_space>;

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

  // The owned+shared column map used by the assembled owned+shared graph.
  RCP<const map_type> ownedPlusSharedColMap;

  if (haveColMap) {
    // -- Path 1: the user provided a column map. --------------------------
    // FILL directly into a locally-indexed entries array, mapping global node
    // IDs to local column indices via the provided column map inside the
    // functor.
    ownedPlusSharedColMap = ownedPlusSharedColMap_;
    auto localColMap      = ownedPlusSharedColMap->getLocalMap();

    entries_t localEntries(Kokkos::ViewAllocateWithoutInitializing("localEntries"), nnz);
    runFill(localEntries, localColMap, std::false_type{});

    local_graph_type localGraph(localEntries, localRowptrs);
    Tpetra::Import_Util::sortCrsGraph(localGraph);

    // Build a fill-complete owned+shared graph directly from the local graph.
    // Using the 5-argument (lclGraph, row, col, domain, range) constructor sets
    // domain = owned and range = owned+shared without any column reindexing, so
    // the local column indices we just computed stay valid.  The owned domain
    // map is what makes CrsGraph build the owned -> owned+shared Importer that
    // the fused Export below relies on.
    ownedPlusSharedGraph_ =
        rcp(new crs_graph_type(localGraph,
                               ownedPlusSharedMap_,   // row map
                               ownedPlusSharedColMap, // column map
                               rowMap_,               // domain map (owned)
                               ownedPlusSharedMap_)); // range map
  } else {
    // -- Path 2: no column map provided. ----------------------------------
    // FILL a globally-indexed entries array, then build the column map with
    // Tpetra::Details::makeColMap and remap to local indices.
    gbl_entries_t globalEntries(Kokkos::ViewAllocateWithoutInitializing("globalEntries"), nnz);
    // For the global path, outColMap is unused; pass the owned+shared local map
    // as a placeholder of the right type.
    runFill(globalEntries, localOwnedPlusSharedMap, std::true_type{});

    // Build the column map from the (global) column indices.  The domain map is
    // the owned row map: makeColMap will locally-fit it into the column map.
    {
      RCP<const map_type> colMap;
      const int err = Tpetra::Details::makeColMap<local_ordinal_type, global_ordinal_type, node_type>(
          colMap, rowMap_, globalEntries);
      TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error,
                                 "Tpetra::Details::GraphAssembly::build: makeColMap failed.");
      ownedPlusSharedColMap = colMap;
    }

    // Remap the global column indices to local indices w.r.t. the new colMap.
    entries_t localEntries(Kokkos::ViewAllocateWithoutInitializing("localEntries"), nnz);
    {
      auto localColMap = ownedPlusSharedColMap->getLocalMap();
      Kokkos::parallel_for(
          "remapGlobalToLocalColumnIndices",
          range_policy(0, nnz),
          KOKKOS_LAMBDA(size_t i) {
            localEntries(i) = localColMap.getLocalElement(globalEntries(i));
          });
    }

    local_graph_type localGraph(localEntries, localRowptrs);
    Tpetra::Import_Util::sortCrsGraph(localGraph);

    ownedPlusSharedGraph_ =
        rcp(new crs_graph_type(localGraph,
                               ownedPlusSharedMap_,   // row map
                               ownedPlusSharedColMap, // column map
                               rowMap_,               // domain map (owned)
                               ownedPlusSharedMap_)); // range map
  }

  ownedPlusSharedColMap_ = ownedPlusSharedColMap;

  // -- Final step: build the owned graph. --
  //
  // For compatibility with FECrsMatrix (whose owned matrix aliases the first
  // chunk of the owned+shared matrix's values), the owned graph must be a
  // first-chunk subview of the owned+shared graph, sharing its column map.  For
  // that to be correct, the owned rows of the owned+shared graph must contain
  // the fully MERGED structure (the union of every process' contributions),
  // while the shared rows keep only this process' partial contributions (so the
  // matrix Export with ADD does not double count).  This mirrors what
  // FECrsGraph::doOwnedPlusSharedToOwned does.
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
    //     gives the complete owned rows (union of all processes' contributions),
    //     with a column map built by the fused path.
    graph_ = Tpetra::exportAndFillCompleteCrsGraph<crs_graph_type>(
        ownedPlusSharedGraph_, *exporter, rowMap_, rowMap_);

    // (b) Build the owned+shared graph by Importing the merged owned graph onto
    //     the owned+shared row map.  Every owned+shared row then holds the
    //     complete (merged) structure of the corresponding owned row, and both
    //     graphs share the same column map.  This makes the owned graph a proper
    //     first-chunk subview of the owned+shared graph, which is exactly what
    //     FECrsMatrix's value aliasing requires.  Extra structural positions in
    //     the shared rows simply receive zero contributions during assembly, so
    //     the Export-with-ADD of the matrix values is still correct.
    using import_t = ::Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type>;
    RCP<const import_t> ownedToOwnedPlusShared =
        rcp(new import_t(rowMap_, ownedPlusSharedMap_));
    Teuchos::RCP<const crs_graph_type> mergedOwnedConst = graph_;
    RCP<crs_graph_type> mergedOwnedPlusShared =
        Tpetra::importAndFillCompleteCrsGraph<crs_graph_type>(
            mergedOwnedConst, *ownedToOwnedPlusShared,
            ownedPlusSharedMap_,  // domainMap of the owned+shared graph
            ownedPlusSharedMap_); // rangeMap of the owned+shared graph

    // (c) Rebuild the owned+shared graph with a column map whose first chunk is
    //     exactly the owned+shared map (same GIDs, same local ordering), followed
    //     by any extra columns that the merge pulled in from across process
    //     boundaries.  Keeping the owned+shared map as the column-map prefix
    //     means matrix values and the right-hand side -- which the caller
    //     assembles using owned+shared local IDs -- land in the right place,
    //     while the appended columns capture the full merged owned-row structure.
    {
      auto mergedLocal   = mergedOwnedPlusShared->getLocalGraphDevice();
      auto mergedRowptrs = mergedLocal.row_map;
      auto mergedEntries = mergedLocal.entries;
      auto mergedColMap  = mergedOwnedPlusShared->getColMap()->getLocalMap();
      const size_type osNNZ = static_cast<size_type>(mergedEntries.extent(0));
      auto localOwnedPlusSharedMap2 = ownedPlusSharedMap_->getLocalMap();
      const local_ordinal_type LO_INVALID2 =
          Teuchos::OrdinalTraits<local_ordinal_type>::invalid();

      // Global column IDs of the merged owned+shared structure.
      gbl_entries_t osGlobal(Kokkos::ViewAllocateWithoutInitializing("osGlobal"), osNNZ);
      Kokkos::parallel_for(
          "owned+sharedLocalToGlobal", range_policy(0, osNNZ),
          KOKKOS_LAMBDA(size_type i) {
            osGlobal(i) = mergedColMap.getGlobalElement(mergedEntries(i));
          });

      // Collect (on host) the set of extra column GIDs that are NOT already in
      // the owned+shared map, preserving a deterministic (sorted) order.
      auto osGlobalHost = Kokkos::create_mirror_view(osGlobal);
      Kokkos::deep_copy(osGlobalHost, osGlobal);
      std::set<global_ordinal_type> extraSet;
      for (size_type i = 0; i < osNNZ; ++i) {
        const global_ordinal_type gid = osGlobalHost(i);
        if (localOwnedPlusSharedMap2.getLocalElement(gid) == LO_INVALID2)
          extraSet.insert(gid);
      }

      // Build the column-map GID list: owned+shared GIDs first (in their local
      // order), then the extra GIDs.
      const local_ordinal_type numOS =
          static_cast<local_ordinal_type>(ownedPlusSharedMap_->getLocalNumElements());
      Teuchos::Array<global_ordinal_type> colGIDs(numOS + extraSet.size());
      {
        auto osGidsHost = ownedPlusSharedMap_->getMyGlobalIndices();
        for (local_ordinal_type i = 0; i < numOS; ++i)
          colGIDs[i] = osGidsHost(i);
        local_ordinal_type k = numOS;
        for (const global_ordinal_type gid : extraSet)
          colGIDs[k++] = gid;
      }
      const Tpetra::global_size_t GST_INVALID =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      RCP<const map_type> osColMap =
          rcp(new map_type(GST_INVALID, colGIDs(), ownedPlusSharedMap_->getIndexBase(),
                           ownedPlusSharedMap_->getComm()));

      // Remap the structure to local indices w.r.t. this column map.
      rowptrs_t osRowptrs("osRowptrs", mergedRowptrs.extent(0));
      Kokkos::deep_copy(osRowptrs, mergedRowptrs);
      entries_t osLocal(Kokkos::ViewAllocateWithoutInitializing("osLocal"), osNNZ);
      {
        auto localOsColMap = osColMap->getLocalMap();
        Kokkos::parallel_for(
            "owned+sharedGlobalToLocal", range_policy(0, osNNZ),
            KOKKOS_LAMBDA(size_type i) {
              osLocal(i) = localOsColMap.getLocalElement(osGlobal(i));
            });
      }

      local_graph_type osGraph(osLocal, osRowptrs);
      Tpetra::Import_Util::sortCrsGraph(osGraph);

      ownedPlusSharedGraph_ =
          rcp(new crs_graph_type(osGraph,
                                 ownedPlusSharedMap_,  // row map
                                 osColMap,             // column map
                                 rowMap_,              // domain map (owned)
                                 ownedPlusSharedMap_));// range map
      ownedPlusSharedColMap_ = osColMap;
    }

    // (d) The clean merged owned graph (from the fused Export) is what getGraph()
    //     returns.  It is already stored in graph_.
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
