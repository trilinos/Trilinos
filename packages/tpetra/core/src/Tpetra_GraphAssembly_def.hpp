// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_GRAPHASSEMBLY_DEF_HPP
#define TPETRA_GRAPHASSEMBLY_DEF_HPP

/// \file Tpetra_GraphAssembly_def.hpp
///
/// Definition (implementation) of Tpetra::Experimental::GraphAssembly.

#include "Tpetra_GraphAssembly_decl.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include "Kokkos_UnorderedMap.hpp"

namespace Tpetra {
namespace Experimental {

namespace Details {

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
                GlobalOrdinal nei    = ownedElementToNode(localElement, j);
                // Work with the LOCAL column ID; the hash table stores LIDs.
                LocalOrdinal neiLid  = columnMap.getLocalElement(nei);
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
/// neighbor is discovered its LID is written into the entries array.
template <class ElementToNode, class RowptrsType, class EntriesType,
          class LocalMapType, class ScratchHashTable, class ScratchCounter,
          class GlobalEdgeSet, class FlagView, class TeamMember,
          class LocalOrdinal, class GlobalOrdinal>
struct FEFillEntriesFunctor {
  using size_type = typename RowptrsType::value_type;

  FEFillEntriesFunctor(
      const RowptrsType& rowptrs_,
      const EntriesType& entries_,
      const ElementToNode& ownedElementToNode_,
      const RowptrsType& nodeToElementRowptrs_,
      const EntriesType& nodeToElementEntries_,
      const LocalMapType& columnMap_,
      GlobalEdgeSet fallbackEdgeSet_,
      const FlagView& globalFailFlag_,
      LocalOrdinal hashSize_)
    : rowptrs(rowptrs_)
    , entries(entries_)
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
                GlobalOrdinal nei    = ownedElementToNode(localElement, j);
                // Work with the LOCAL column ID; the hash table stores LIDs,
                // which are exactly what gets written into the entries array.
                LocalOrdinal neiLid  = columnMap.getLocalElement(nei);
                size_t hash          = neiLid;
                bool foundOrInserted = false;
                for (unsigned probe = 0; probe < 8; probe++) {
                  hash         = KokkosKernels::Impl::xorshiftHash(hash);
                  unsigned pos = (hash + probe) % hashSize;
                  while (true) {
                    LocalOrdinal cellValue = Kokkos::volatile_load(&hashTable(pos));
                    if (cellValue == neiLid) {
                      foundOrInserted = true;
                      break;
                    } else if (cellValue == LO_INVALID) {
                      if (LO_INVALID == Kokkos::atomic_compare_exchange(&hashTable(pos), LO_INVALID, neiLid)) {
                        foundOrInserted              = true;
                        LocalOrdinal insertPos       = Kokkos::atomic_fetch_add(&insertCounter(), 1);
                        entries(rowptrs(localRow) + insertPos) = neiLid;
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
                  size_t edge       = size_t(neiLid) * size_t(numLocalNodes) + localRow;
                  auto insertResult = fallbackEdgeSet.insert(edge);
                  if (insertResult.failed()) {
                    globalFailFlag() = 1;
                  }
                  if (!insertResult.existing()) {
                    LocalOrdinal insertPos       = Kokkos::atomic_fetch_add(&insertCounter(), 1);
                    entries(rowptrs(localRow) + insertPos) = neiLid;
                  }
                }
              });
        });
  }

  RowptrsType rowptrs;
  EntriesType entries;
  ElementToNode ownedElementToNode;
  RowptrsType nodeToElementRowptrs;
  EntriesType nodeToElementEntries;
  LocalMapType columnMap;
  GlobalEdgeSet fallbackEdgeSet;
  FlagView globalFailFlag;
  LocalOrdinal hashSize;
};

}  // namespace Details

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::GraphAssembly(
    const Teuchos::RCP<const map_type>& rowMap,
    const Teuchos::RCP<const map_type>& ownedPlusSharedMap,
    const element_to_node_type& ownedElementToNode)
  : rowMap_(rowMap)
  , ownedPlusSharedMap_(ownedPlusSharedMap)
  , ownedElementToNode_(ownedElementToNode) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GraphAssembly<LocalOrdinal, GlobalOrdinal, Node>::build() {
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Local types used only by the kernels below.
  using local_graph_type    = typename crs_graph_type::local_graph_device_type;
  using rowptrs_t           = typename local_graph_type::row_map_type::non_const_type;
  using entries_t           = typename local_graph_type::entries_type;
  using size_type           = typename rowptrs_t::value_type;
  using local_map_type      = typename map_type::local_map_type;

  using range_policy    = Kokkos::RangePolicy<execution_space>;
  using team_policy     = Kokkos::TeamPolicy<execution_space>;
  using team_member     = typename team_policy::member_type;
  using scratch_space   = typename execution_space::scratch_memory_space;
  using scratch_hash_table =
      Kokkos::View<local_ordinal_type*, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using scratch_counter =
      Kokkos::View<local_ordinal_type, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using global_edge_set = Kokkos::UnorderedMap<size_t, void, device_type>;
  using flag_view       = Kokkos::View<int, typename device_type::memory_space>;

  const local_ordinal_type LO_INVALID = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();

  auto owned_element_to_node_ids = ownedElementToNode_;

  // The assembled graph's rows AND columns are both indexed by the owned+shared
  // map: an owned element's nodes are, by construction, all either owned or
  // shared on this process, so there are no remote columns to worry about at
  // this stage.  That means no makeColMap call is needed.
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
        local_ordinal_type nodeOfElem             = i % nodesPerElement;
        global_ordinal_type globalNode            = owned_element_to_node_ids(ownedElementIndex, nodeOfElem);
        local_ordinal_type localNode              = localOwnedPlusSharedMap.getLocalElement(globalNode);
        if (localNode != LO_INVALID)
          Kokkos::atomic_inc(&nodeToElementRowptrs(localNode));
      });
  typename rowptrs_t::value_type nodeToElementNNZ;
  Kokkos::parallel_scan(
      range_policy(0, numLocalNodes + 1),
      KOKKOS_LAMBDA(local_ordinal_type i, size_type& lsum, bool finalPass) {
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
          local_ordinal_type nodeOfElem             = i % nodesPerElement;
          global_ordinal_type globalNode            = owned_element_to_node_ids(ownedElementIndex, nodeOfElem);
          local_ordinal_type localNode              = localOwnedPlusSharedMap.getLocalElement(globalNode);
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
        Details::FECountEntriesFunctor<element_to_node_type, rowptrs_t, entries_t,
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
      KOKKOS_LAMBDA(local_ordinal_type i, size_type& lsum, bool finalPass) {
        size_type count = localRowptrs(i);
        if (finalPass)
          localRowptrs(i) = lsum;
        lsum += count;
      },
      nnz);

  // -- Step 5: FILL pass: write the local column indices --
  entries_t localEntries(Kokkos::ViewAllocateWithoutInitializing("localEntries"), nnz);
  {
    using fill_functor_type =
        Details::FEFillEntriesFunctor<element_to_node_type, rowptrs_t, entries_t,
                                       local_map_type, scratch_hash_table, scratch_counter,
                                       global_edge_set, flag_view, team_member,
                                       local_ordinal_type, global_ordinal_type>;
    while (true) {
      global_edge_set fallbackSet(fallbackSetSize);
      fill_functor_type functor(
          localRowptrs, localEntries, owned_element_to_node_ids,
          nodeToElementRowptrs, nodeToElementEntries,
          localOwnedPlusSharedMap, fallbackSet, failFlag, hashTableSize);
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
        fallbackSetSize *= 2;
        Kokkos::deep_copy(failFlag, 0);
      } else
        break;
    }
  }

  // -- Step 6: wrap in a CrsGraph and communicate explicitly --
  local_graph_type localGraph(localEntries, localRowptrs);
  Tpetra::Import_Util::sortCrsGraph(localGraph);

  // The locally-assembled graph.  Row map == col map == owned+shared map, but
  // the DOMAIN map is the owned map: that is what makes this a meaningful
  // "assembly" graph (overlapping rows, one-to-one domain), and it makes
  // CrsGraph build the owned -> owned+shared Importer that both the fused
  // Export below and the RHS assembly rely on.  This graph is fill-complete
  // immediately (packed, locally indexed); no insertGlobalIndices allocation
  // happens anywhere.  Its rows for shared nodes hold this process' partial
  // contributions to rows owned elsewhere.
  ownedPlusSharedGraph_ =
      rcp(new crs_graph_type(localGraph,
                             ownedPlusSharedMap_,  // row map
                             ownedPlusSharedMap_,  // column map
                             rowMap_,              // domain map (owned)
                             ownedPlusSharedMap_)  // range map
      );

  if (ownedPlusSharedMap_->isSameAs(*rowMap_)) {
    // Serial / single-rank case: nothing to communicate.  The assembly graph
    // already has exactly the owned structure; just re-map domain/range.
    ownedPlusSharedGraph_->resumeFill();
    ownedPlusSharedGraph_->fillComplete(rowMap_, rowMap_);
    graph_    = ownedPlusSharedGraph_;
    exporter_ = Teuchos::null;
  } else {
    // Parallel case: Export from the overlapping owned+shared row map to the
    // one-to-one owned row map.  This is exactly the communication that
    // FECrsGraph::endAssembly() would have performed, but here it is a fused
    // Export + fillComplete on a packed graph:
    //   - rows that are shared get sent to their owning process,
    //   - the destination graph's column map is built by the fused path
    //     (lowCommunicationMakeColMapAndReindex), so it correctly picks up the
    //     cross-boundary columns that arrive from other processes,
    //   - duplicate column entries are sorted and merged.
    RCP<const export_type> exporter = rcp(new export_type(ownedPlusSharedMap_, rowMap_));
    exporter_                       = exporter;
    graph_                          = Tpetra::exportAndFillCompleteCrsGraph<crs_graph_type>(
        ownedPlusSharedGraph_, *exporter, rowMap_, rowMap_);
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

}  // namespace Experimental
}  // namespace Tpetra

//
// Explicit template instantiation macro for
// Tpetra::Experimental::GraphAssembly.  NOT FOR USERS!!!  Must be used
// inside the Tpetra namespace.
//
#define TPETRA_GRAPHASSEMBLY_INSTANT(LO, GO, NODE) \
  template class Experimental::GraphAssembly<LO, GO, NODE>;

#endif  // TPETRA_GRAPHASSEMBLY_DEF_HPP
