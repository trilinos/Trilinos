// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_TOTALELEMENTLOOP_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_TOTALELEMENTLOOP_HPP

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "Tpetra_Core.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "KokkosSparse_SortCrs.hpp"

#include "fem_assembly_typedefs.hpp"
#include "fem_assembly_MeshDatabase.hpp"
#include "fem_assembly_Element.hpp"
#include "fem_assembly_utility.hpp"
#include "fem_assembly_commandLineOpts.hpp"

namespace TpetraExamples {

int executeTotalElementLoopSP_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                               const struct CmdLineOpts& opts);

int executeTotalElementLoopSPKokkos_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                     const struct CmdLineOpts& opts);

int executeTotalElementLoopSP(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                              const struct CmdLineOpts& opts) {
  using Teuchos::RCP;

  // The output stream 'out' will ignore any output not from Process 0.
  RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
  Teuchos::FancyOStream& out      = *pOut;

  std::string useKokkos = opts.useKokkosAssembly ? "Kokkos Assembly" : "Serial Assembly";
  out << "================================================================================" << std::endl
      << "=  Total Element Loop (Static Profile; " << useKokkos << ")" << std::endl
      << "================================================================================" << std::endl
      << std::endl;

  int status = 0;
  for (size_t i = 0; i < opts.repetitions; ++i) {
    if (opts.useKokkosAssembly)
      status += executeTotalElementLoopSPKokkos_(comm, opts);
    else
      status += executeTotalElementLoopSP_(comm, opts);
  }
  return status;
}

int executeTotalElementLoopSP_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                               const struct CmdLineOpts& opts) {
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  const global_ordinal_type GO_INVALID = Teuchos::OrdinalTraits<global_ordinal_type>::invalid();

  // The output stream 'out' will ignore any output not from Process 0.
  RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
  Teuchos::FancyOStream& out      = *pOut;

  // Processor decomp (only works on perfect squares)
  int numProcs  = comm->getSize();
  int sqrtProcs = sqrt(numProcs);

  if (sqrtProcs * sqrtProcs != numProcs) {
    if (0 == comm->getRank())
      std::cerr << "Error: Invalid number of processors provided, num processors must be a perfect square." << std::endl;
    return -1;
  }
  int procx = sqrtProcs;
  int procy = sqrtProcs;

  // Generate a simple 3x3 mesh
  int nex = opts.numElementsX;
  int ney = opts.numElementsY;

  MeshDatabase mesh(comm, nex, ney, procx, procy);

  if (opts.verbose) mesh.print(std::cout);

  int maxEntriesPerRow = 16;

  // Build Tpetra Maps
  // -----------------
  // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
  RCP<const map_type> row_map            = rcp(new map_type(GO_INVALID, mesh.getOwnedNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
  RCP<const map_type> owned_element_map  = rcp(new map_type(GO_INVALID, mesh.getOwnedElementGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
  RCP<const map_type> ghost_element_map  = rcp(new map_type(GO_INVALID, mesh.getGhostElementGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
  RCP<const import_type> elementImporter = rcp(new import_type(owned_element_map, ghost_element_map));

  if (opts.verbose) row_map->describe(out);

  // Graph Construction
  // ------------------
  // - Loop over every element in the mesh.
  //   - Get list of nodes associated with each element.
  //   - Insert node contributions if the node (row) is owned locally.
  //
  auto domain_map = row_map;
  auto range_map  = row_map;

  auto owned_element_to_node_ids = mesh.getOwnedElementToNode().getHostView(Tpetra::Access::ReadOnly);
  auto ghost_element_to_node_ids = mesh.getGhostElementToNode().getHostView(Tpetra::Access::ReadOnly);

  Teuchos::TimeMonitor::getStackedTimer()->startBaseTimer();
  RCP<TimeMonitor> timerElementLoopGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("1) ElementLoop  (Graph)")));

  RCP<crs_graph_type> crs_graph = rcp(new crs_graph_type(row_map, maxEntriesPerRow));

  // Using 4 because we're using quads for this example, so there will be 4 nodes associated with each element.
  Teuchos::Array<global_ordinal_type> global_ids_in_row(4);

  // Insert node contributions for every OWNED element:
  for (size_t element_gidx = 0; element_gidx < mesh.getNumOwnedElements(); element_gidx++) {
    // Populate global_ids_in_row:
    // - Copy the global node ids for current owned element into an array.
    // - Since each element's contribution is a clique, we can re-use this for
    //   each row associated with this element's contribution.
    for (size_t element_node_idx = 0; element_node_idx < owned_element_to_node_ids.extent(1); element_node_idx++) {
      global_ids_in_row[element_node_idx] = owned_element_to_node_ids(element_gidx, element_node_idx);
    }

    // Add the contributions from the current row into the graph if the node is owned.
    // - For example, if Element 0 contains nodes [0,1,4,5] and nodes 0 and 4 are owned
    //   by the current processor, then:
    //   - node 0 inserts [0, 1, 4, 5]
    //   - node 1 <skip>
    //   - node 4 inserts [0, 1, 4, 5]
    //   - node 5 <skip>
    for (size_t element_node_idx = 0; element_node_idx < owned_element_to_node_ids.extent(1); element_node_idx++) {
      if (mesh.nodeIsOwned(global_ids_in_row[element_node_idx])) {
        crs_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
      }
    }
  }

  // Insert the node contributions for every GHOST element:
  for (size_t element_gidx = 0; element_gidx < mesh.getNumGhostElements(); element_gidx++) {
    for (size_t element_node_idx = 0; element_node_idx < ghost_element_to_node_ids.extent(1); element_node_idx++) {
      global_ids_in_row[element_node_idx] = ghost_element_to_node_ids(element_gidx, element_node_idx);
    }
    for (size_t element_node_idx = 0; element_node_idx < ghost_element_to_node_ids.extent(1); element_node_idx++) {
      if (mesh.nodeIsOwned(global_ids_in_row[element_node_idx])) {
        crs_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
      }
    }
  }

  timerElementLoopGraph = Teuchos::null;

  // 'finalize' the crs_graph by calling fillComplete().
  {
    TimeMonitor timer(*TimeMonitor::getNewTimer("2) FillComplete (Graph)"));
    crs_graph->fillComplete();
  }

  // Print out the crs_graph in detail...
  if (opts.verbose) crs_graph->describe(out, Teuchos::VERB_EXTREME);

  // Simulated Ghosting of Material State
  // -------------------
  {
    GhostState state(elementImporter, opts.numStateDoublesPerElement);
    TimeMonitor timer(*TimeMonitor::getNewTimer("3.1) Ghosting Material State (Matrix)"));
    state.doGhost();
  }

  // Matrix Fill
  // -------------------
  // In this example, we're using a simple stencil of values for the stiffness matrix:
  //
  //    +-----+-----+-----+-----+
  //    |  2  | -1  |     | -1  |
  //    +-----+-----+-----+-----+
  //    | -1  |  2  | -1  |     |
  //    +-----+-----+-----+-----+
  //    |     | -1  |  2  | -1  |
  //    +-----+-----+-----+-----+
  //    | -1  |     | -1  |  2  |
  //    +-----+-----+-----+-----+
  //
  // For matrix fill, we create the crs_matrix object and will fill it
  // in the same manner as we filled in the graph but in this case, nodes
  // associated with each element will receive contributions according to
  // the row in this stencil.
  //
  // In this example, the calls to sumIntoGlobalValues() on 1 core will look like:
  //   Element 0
  // - sumIntoGlobalValues( 0,  [  0  1  5  4  ],  [  2  -1  0  -1  ])
  // - sumIntoGlobalValues( 1,  [  0  1  5  4  ],  [  -1  2  -1  0  ])
  // - sumIntoGlobalValues( 5,  [  0  1  5  4  ],  [  0  -1  2  -1  ])
  // - sumIntoGlobalValues( 4,  [  0  1  5  4  ],  [  -1  0  -1  2  ])
  // Element 1
  // - sumIntoGlobalValues( 1,  [  1  2  6  5  ],  [  2  -1  0  -1  ])
  // - sumIntoGlobalValues( 2,  [  1  2  6  5  ],  [  -1  2  -1  0  ])
  // - sumIntoGlobalValues( 6,  [  1  2  6  5  ],  [  0  -1  2  -1  ])
  // - sumIntoGlobalValues( 5,  [  1  2  6  5  ],  [  -1  0  -1  2  ])
  // Element 2
  // - sumIntoGlobalValues( 2,  [  2  3  7  6  ],  [  2  -1  0  -1  ])
  // - sumIntoGlobalValues( 3,  [  2  3  7  6  ],  [  -1  2  -1  0  ])
  // - sumIntoGlobalValues( 7,  [  2  3  7  6  ],  [  0  -1  2  -1  ])
  // - sumIntoGlobalValues( 6,  [  2  3  7  6  ],  [  -1  0  -1  2  ])
  //
  // Similarly to the Graph construction above, we loop over both local and global
  // elements and insert rows for only the locally owned rows.
  //
  RCP<TimeMonitor> timerElementLoopMatrix = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("3) ElementLoop  (Matrix)")));

  crs_matrix_type crs_matrix(crs_graph);
  multivector_type rhs(crs_graph->getRowMap(), 1);

  Scalar element_matrix[4][4];
  Scalar element_rhs[4];

  Teuchos::Array<global_ordinal_type> column_global_ids(4);  // global column ids list
  Teuchos::Array<Scalar> column_scalar_values(4);            // scalar values for each column

  // Loop over owned elements:
  for (size_t element_gidx = 0;
       element_gidx < mesh.getNumOwnedElements(); ++element_gidx) {
    // Get the stiffness matrix for this element
    ReferenceQuad4(element_matrix, element_rhs);

    // Fill the global column ids array for this element
    for (size_t element_node_idx = 0;
         element_node_idx < owned_element_to_node_ids.extent(1);
         ++element_node_idx) {
      column_global_ids[element_node_idx] =
          owned_element_to_node_ids(element_gidx, element_node_idx);
    }

    // For each node (row) on the current element:
    // - populate the values array
    // - add values to crs_matrix if the row is owned.
    //   Note: hardcoded 4 here because we're using quads.
    for (size_t element_node_idx = 0; element_node_idx < 4;
         ++element_node_idx) {
      const global_ordinal_type global_row_id =
          owned_element_to_node_ids(element_gidx, element_node_idx);
      if (mesh.nodeIsOwned(global_row_id)) {
        for (size_t col_idx = 0; col_idx < 4; ++col_idx) {
          column_scalar_values[col_idx] =
              element_matrix[element_node_idx][col_idx];
        }
        crs_matrix.sumIntoGlobalValues(global_row_id,
                                       column_global_ids,
                                       column_scalar_values);
        rhs.sumIntoGlobalValue(global_row_id, 0,
                               element_rhs[element_node_idx]);
      }
    }
  }

  // Loop over ghost elements:
  // - This loop is the same as the element loop for owned elements, but this one
  //   is for ghost elements.
  for (size_t element_gidx = 0; element_gidx < mesh.getNumGhostElements(); element_gidx++) {
    ReferenceQuad4(element_matrix, element_rhs);

    for (size_t element_node_idx = 0; element_node_idx < ghost_element_to_node_ids.extent(1); element_node_idx++) {
      column_global_ids[element_node_idx] = ghost_element_to_node_ids(element_gidx, element_node_idx);
    }

    for (size_t element_node_idx = 0; element_node_idx < 4; element_node_idx++) {
      global_ordinal_type global_row_id = ghost_element_to_node_ids(element_gidx, element_node_idx);
      if (mesh.nodeIsOwned(global_row_id)) {
        for (size_t col_idx = 0; col_idx < 4; col_idx++) {
          column_scalar_values[col_idx] = element_matrix[element_node_idx][col_idx];
        }
        crs_matrix.sumIntoGlobalValues(global_row_id, column_global_ids, column_scalar_values);
        rhs.sumIntoGlobalValue(global_row_id, 0, element_rhs[element_node_idx]);
      }
    }
  }
  timerElementLoopMatrix = Teuchos::null;

  // After the contributions are added, 'finalize' the matrix using fillComplete()
  {
    TimeMonitor timer(*TimeMonitor::getNewTimer("4) FillComplete (Matrix)"));
    crs_matrix.fillComplete();
  }

  // Print out crs_matrix details.
  if (opts.verbose) {
    crs_matrix.describe(out, Teuchos::VERB_EXTREME);
  }

  Teuchos::TimeMonitor::getStackedTimer()->stopBaseTimer();

  // Save crs_matrix as a MatrixMarket file.
  if (opts.saveMM) {
    using writer_type = Tpetra::MatrixMarket::Writer<crs_matrix_type>;
    std::ofstream ofs("crsMatrix_TotalElementLoop_SP.out", std::ofstream::out);
    writer_type::writeSparse(ofs, crs_matrix);
    std::ofstream ofs2("rhs_TotalElementLoop_SP.out", std::ofstream::out);
    writer_type::writeDense(ofs2, rhs);
  }

  return 0;
}

struct PrefixSumCountsFunctor {
  PrefixSumCountsFunctor(const rowptrs_t& localRowptrs_)
    : localRowptrs(localRowptrs_) {}

  KOKKOS_INLINE_FUNCTION void operator()(local_ordinal_type i, size_type& lsum, bool finalPass) const {
    size_type count = localRowptrs(i);
    if (finalPass)
      localRowptrs(i) = lsum;
    lsum += count;
  }

  rowptrs_t localRowptrs;
};

template <typename ElementToNode>
struct CountEntriesFunctor {
  CountEntriesFunctor(
      const rowptrs_t& counts_,
      const ElementToNode& ownedElementToNode_,
      const ElementToNode& ghostElementToNode_,
      const rowptrs_t& nodeToElementRowptrs_,
      const entries_t& nodeToElementEntries_,
      global_edge_set fallbackEdgeSet_,
      const flag_view& globalFailFlag_,
      local_ordinal_type hashSize_)
    : counts(counts_)
    , ownedElementToNode(ownedElementToNode_)
    , ghostElementToNode(ghostElementToNode_)
    , nodeToElementRowptrs(nodeToElementRowptrs_)
    , nodeToElementEntries(nodeToElementEntries_)
    , fallbackEdgeSet(fallbackEdgeSet_)
    , globalFailFlag(globalFailFlag_)
    , hashSize(hashSize_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& t) const {
    const global_ordinal_type GO_INVALID = ~global_ordinal_type(0);
    scratch_hash_table hashTable(t.team_scratch(0), hashSize);
    // Initially, fill the table with INVALID to mean 'empty'.
    Kokkos::parallel_for(Kokkos::TeamVectorRange(t, hashSize),
                         [&](int i) {
                           hashTable(i) = GO_INVALID;
                         });
    t.team_barrier();
    local_ordinal_type localRow            = t.league_rank();
    const local_ordinal_type numLocalNodes = counts.extent(0) - 1;
    size_type elementBegin                 = nodeToElementRowptrs(localRow);
    size_type elementEnd                   = nodeToElementRowptrs(localRow + 1);
    // Each thread will count entries, and then reduce at the end.
    size_type numEntries;
    // Iterate over adjacent elements
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(t, elementBegin, elementEnd),
        [&](size_type i, size_type& lTeamCount) {
          local_ordinal_type localElement = nodeToElementEntries(i);
          size_type numThreadEntries;
          // Iterate over nodes adjacent to localElement
          Kokkos::parallel_reduce(
              Kokkos::ThreadVectorRange(t, ownedElementToNode.extent(1)),
              [&](int j, size_type& lThreadCount) {
                global_ordinal_type nei;
                if (localElement < ownedElementToNode.extent_int(0))
                  nei = ownedElementToNode(localElement, j);
                else
                  nei = ghostElementToNode(localElement - ownedElementToNode.extent_int(0), j);
                // Try to insert nei into the scratch hash table, if it's not already there.
                // TODO: xorshift* is overkill for this, plain xorshift (without the multiplication) would be enough.
                size_t hash          = nei;
                bool foundOrInserted = false;
                for (unsigned probe = 0; probe < 8; probe++) {
                  hash         = KokkosKernels::Impl::xorshiftHash(hash);
                  unsigned pos = hash % hashSize;
                  while (true) {
                    global_ordinal_type cellValue = Kokkos::volatile_load(&hashTable(pos));
                    if (cellValue == nei) {
                      foundOrInserted = true;
                      break;
                    } else if (cellValue == GO_INVALID) {
                      if (GO_INVALID == Kokkos::atomic_compare_exchange(&hashTable(pos), GO_INVALID, nei)) {
                        foundOrInserted = true;
                        lThreadCount++;
                        break;
                      }
                      // If the cmp-xchg failed, we don't know yet what was inserted at pos by some other thread.
                      // So we have to go back to the top of the while and try again.
                    } else {
                      // this cell is populated with an entry other than nei
                      break;
                    }
                  }
                  if (foundOrInserted)
                    break;
                }
                if (!foundOrInserted) {
                  // Scratch hash table ran out of space - must put edge in the global fallback set.
                  size_t edge       = size_t(nei) * size_t(numLocalNodes) + localRow;
                  auto insertResult = fallbackEdgeSet.insert(edge);
                  if (insertResult.failed()) {
                    // the global table filled up, so we'll have to make it bigger and start over
                    globalFailFlag() = 1;
                  }
                  if (!insertResult.existing()) {
                    // This was a new edge
                    lThreadCount++;
                  }
                }
              },
              numThreadEntries);
          Kokkos::single(Kokkos::PerThread(t),
                         [&]() {
                           lTeamCount += numThreadEntries;
                         });
        },
        numEntries);
    // Write out the final count
    Kokkos::single(Kokkos::PerTeam(t),
                   [&]() {
                     counts(localRow) = numEntries;
                   });
  }

  rowptrs_t counts;
  ElementToNode ownedElementToNode;
  ElementToNode ghostElementToNode;
  rowptrs_t nodeToElementRowptrs;
  entries_t nodeToElementEntries;
  global_edge_set fallbackEdgeSet;
  flag_view globalFailFlag;
  local_ordinal_type hashSize;
};

template <typename ElementToNode>
struct FillEntriesFunctor {
  FillEntriesFunctor(
      const rowptrs_t& rowptrs_,
      const entries_t& entries_,
      const ElementToNode& ownedElementToNode_,
      const ElementToNode& ghostElementToNode_,
      const rowptrs_t& nodeToElementRowptrs_,
      const entries_t& nodeToElementEntries_,
      const local_map_type& columnMap_,
      global_edge_set fallbackEdgeSet_,
      const flag_view& globalFailFlag_,
      local_ordinal_type hashSize_)
    : rowptrs(rowptrs_)
    , entries(entries_)
    , ownedElementToNode(ownedElementToNode_)
    , ghostElementToNode(ghostElementToNode_)
    , nodeToElementRowptrs(nodeToElementRowptrs_)
    , nodeToElementEntries(nodeToElementEntries_)
    , columnMap(columnMap_)
    , fallbackEdgeSet(fallbackEdgeSet_)
    , globalFailFlag(globalFailFlag_)
    , hashSize(hashSize_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const team_member& t) const {
    const global_ordinal_type GO_INVALID = ~global_ordinal_type(0);
    scratch_hash_table hashTable(t.team_scratch(0), hashSize);
    // Use this rank-0 size_type view to atomically count up
    // positions in the row as entries are inserted.
    scratch_counter insertCounter(t.team_scratch(0));
    Kokkos::single(Kokkos::PerTeam(t),
                   [&]() {
                     insertCounter() = 0;
                   });
    // Initially, fill the table with INVALID to mean 'empty'.
    Kokkos::parallel_for(Kokkos::TeamVectorRange(t, hashSize),
                         [&](int i) {
                           hashTable(i) = GO_INVALID;
                         });
    t.team_barrier();
    const local_ordinal_type numLocalNodes = rowptrs.extent(0) - 1;
    local_ordinal_type localRow            = t.league_rank();
    size_type elementBegin                 = nodeToElementRowptrs(localRow);
    size_type elementEnd                   = nodeToElementRowptrs(localRow + 1);
    // Iterate over adjacent elements
    Kokkos::parallel_for(Kokkos::TeamThreadRange(t, elementBegin, elementEnd),
                         [&](size_type i) {
                           local_ordinal_type localElement = nodeToElementEntries(i);
                           // Iterate over nodes adjacent to localElement
                           Kokkos::parallel_for(Kokkos::ThreadVectorRange(t, ownedElementToNode.extent(1)),
                                                [&](int j) {
                                                  global_ordinal_type nei;
                                                  if (localElement < ownedElementToNode.extent_int(0))
                                                    nei = ownedElementToNode(localElement, j);
                                                  else
                                                    nei = ghostElementToNode(localElement - ownedElementToNode.extent(0), j);
                                                  // Try to insert nei into the scratch hash table, if it's not already there.
                                                  // TODO: xorshift* is overkill for this, plain xorshift (without the multiplication) would be enough.
                                                  size_t hash          = nei;
                                                  bool foundOrInserted = false;
                                                  for (unsigned probe = 0; probe < 8; probe++) {
                                                    hash         = KokkosKernels::Impl::xorshiftHash(hash);
                                                    unsigned pos = (hash + probe) % hashSize;
                                                    while (true) {
                                                      global_ordinal_type cellValue = Kokkos::volatile_load(&hashTable(pos));
                                                      if (cellValue == nei) {
                                                        foundOrInserted = true;
                                                        break;
                                                      } else if (cellValue == GO_INVALID) {
                                                        if (GO_INVALID == Kokkos::atomic_compare_exchange(&hashTable(pos), GO_INVALID, nei)) {
                                                          foundOrInserted = true;
                                                          // Insert LID for nei
                                                          local_ordinal_type insertPos           = Kokkos::atomic_fetch_add(&insertCounter(), 1);
                                                          entries(rowptrs(localRow) + insertPos) = columnMap.getLocalElement(nei);
                                                          break;
                                                        }
                                                        // If the cmp-xchg failed, we don't know yet what was inserted at pos by some other thread.
                                                        // So we have to go back to the top of the while and try again.
                                                      } else {
                                                        // this cell is populated with an entry other than nei
                                                        break;
                                                      }
                                                    }
                                                    if (foundOrInserted)
                                                      break;
                                                  }
                                                  if (!foundOrInserted) {
                                                    // Scratch hash table ran out of space - must put edge in the global fallback set.
                                                    size_t edge       = size_t(nei) * size_t(numLocalNodes) + localRow;
                                                    auto insertResult = fallbackEdgeSet.insert(edge);
                                                    if (insertResult.failed()) {
                                                      // the global table filled up, so we'll have to make it bigger and start over
                                                      globalFailFlag() = 1;
                                                    }
                                                    if (!insertResult.existing()) {
                                                      // This was a new edge, so insert LID into row
                                                      local_ordinal_type insertPos           = Kokkos::atomic_fetch_add(&insertCounter(), 1);
                                                      entries(rowptrs(localRow) + insertPos) = columnMap.getLocalElement(nei);
                                                    }
                                                  }
                                                });
                         });
  }
  rowptrs_t rowptrs;
  entries_t entries;
  ElementToNode ownedElementToNode;
  ElementToNode ghostElementToNode;
  rowptrs_t nodeToElementRowptrs;
  entries_t nodeToElementEntries;
  local_map_type columnMap;
  global_edge_set fallbackEdgeSet;
  flag_view globalFailFlag;
  local_ordinal_type hashSize;
};

int executeTotalElementLoopSPKokkos_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                     const struct CmdLineOpts& opts) {
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  const global_ordinal_type GO_INVALID = Teuchos::OrdinalTraits<global_ordinal_type>::invalid();
  const local_ordinal_type LO_INVALID  = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();
  using pair_type                      = Kokkos::pair<int, int>;

  // The output stream 'out' will ignore any output not from Process 0.
  RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
  Teuchos::FancyOStream& out      = *pOut;

  // Processor decomp (only works on perfect squares)
  int numProcs  = comm->getSize();
  int sqrtProcs = sqrt(numProcs);

  if (sqrtProcs * sqrtProcs != numProcs) {
    if (0 == comm->getRank()) {
      std::cerr << "Error: Invalid number " << numProcs << " of MPI "
                                                           "processes provided.  This number must be a perfect square."
                << std::endl;
    }
    return -1;
  }
  int procx = sqrtProcs;
  int procy = sqrtProcs;

  // Generate a simple 3x3 mesh
  int nex = opts.numElementsX;
  int ney = opts.numElementsY;

  MeshDatabase mesh(comm, nex, ney, procx, procy);

  if (opts.verbose) mesh.print(std::cout);

  // Build Tpetra Maps
  // -----------------
  // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
  RCP<const map_type> row_map =
      rcp(new map_type(GO_INVALID, mesh.getOwnedNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
  RCP<const map_type> owned_element_map;
  RCP<const map_type> ghost_element_map;
  RCP<const import_type> elementImporter;

  if (opts.verbose) row_map->describe(out);

  // Graph Construction
  // ------------------
  // - Loop over every element in the mesh.
  //   - Get list of nodes associated with each element.
  //   - Insert node contributions if the node (row) is owned locally.
  //
  auto domain_map = row_map;
  auto range_map  = row_map;

  // Build the column map first. This is exactly the

  Teuchos::TimeMonitor::getStackedTimer()->startBaseTimer();

  RCP<crs_graph_type> crs_graph;
  {
    RCP<TimeMonitor> timerElementLoopGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("1) ElementLoop  (Graph)")));
    auto owned_element_to_node_ids         = mesh.getOwnedElementToNode().getDeviceView(Tpetra::Access::ReadOnly);
    auto ghost_element_to_node_ids         = mesh.getGhostElementToNode().getDeviceView(Tpetra::Access::ReadOnly);

    // First, build the column map. Each process has the owned GIDs of the domain map, plus all non-owned node GIDs adjacent to ghost elements.
    // Use makeColMap utility function to do this, given a list of GIDs (duplicates allowed)
    // TODO: add a version of makeColMap to expect non-duplicated remote GIDs, and which takes the local GIDs from domain map.
    // For now, we are just passing in the domain map's owned nodes plus all the nodes adjacent to ghost elements (so there will be duplicates)
    RCP<TimeMonitor> colMapTimer = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("1a) Graph Column Map")));
    auto localDomMap             = domain_map->getLocalMap();
    auto numLocalNodes           = domain_map->getLocalNumElements();
    flag_view failFlag("failFlag");
    auto failFlagHost = Kokkos::create_mirror_view(failFlag);
    Kokkos::View<global_ordinal_type*, memory_space> colMapGIDs(Kokkos::ViewAllocateWithoutInitializing("colMapGIDs"), owned_element_to_node_ids.size() + ghost_element_to_node_ids.size());
    Kokkos::parallel_for(
        range_policy(0, owned_element_to_node_ids.size() + ghost_element_to_node_ids.size()),
        KOKKOS_LAMBDA(size_t i) {
          const local_ordinal_type nodesPerElement = owned_element_to_node_ids.extent(1);
          if (i < owned_element_to_node_ids.size()) {
            local_ordinal_type ownedElementIndex = i / nodesPerElement;
            local_ordinal_type nodeOfElem        = i % nodesPerElement;
            colMapGIDs(i)                        = owned_element_to_node_ids(ownedElementIndex, nodeOfElem);
          } else {
            size_t ghostI                        = i - owned_element_to_node_ids.size();
            local_ordinal_type ghostElementIndex = ghostI / nodesPerElement;
            local_ordinal_type nodeOfElem        = ghostI % nodesPerElement;
            colMapGIDs(i)                        = ghost_element_to_node_ids(ghostElementIndex, nodeOfElem);
          }
        });
    execution_space().fence();
    RCP<const map_type> column_map;
    Tpetra::Details::makeColMap(column_map, domain_map, colMapGIDs);
    auto localColMap = column_map->getLocalMap();
    colMapTimer      = Teuchos::null;
    // Now build the transpose of the element to node graph (node -> element)
    // The rows are locally indexed, owned nodes (according to row/domain map)
    // The columns are indices into concatenated owned and ghost element lists. So 0..mesh.getNumOwnedElements() are owned, and the rest are ghost.
    // Unlike the element_to_node graphs, this one has variable entries per row.
    rowptrs_t nodeToElementRowptrs("nodeToElementRowptrs", numLocalNodes + 1);
    Kokkos::parallel_for(
        range_policy(0, owned_element_to_node_ids.size() + ghost_element_to_node_ids.size()),
        KOKKOS_LAMBDA(size_t i) {
          const local_ordinal_type nodesPerElement = owned_element_to_node_ids.extent(1);
          if (i < owned_element_to_node_ids.size()) {
            local_ordinal_type ownedElementIndex = i / nodesPerElement;
            local_ordinal_type nodeOfElem        = i % nodesPerElement;
            global_ordinal_type globalNode       = owned_element_to_node_ids(ownedElementIndex, nodeOfElem);
            local_ordinal_type localNode         = localDomMap.getLocalElement(globalNode);
            if (localNode != LO_INVALID)
              Kokkos::atomic_inc(&nodeToElementRowptrs(localNode));
          } else {
            i -= owned_element_to_node_ids.size();
            local_ordinal_type ghostElementIndex = i / nodesPerElement;
            local_ordinal_type nodeOfElem        = i % nodesPerElement;
            global_ordinal_type globalNode       = ghost_element_to_node_ids(ghostElementIndex, nodeOfElem);
            local_ordinal_type localNode         = localDomMap.getLocalElement(globalNode);
            if (localNode != LO_INVALID)
              Kokkos::atomic_inc(&nodeToElementRowptrs(localNode));
          }
        });
    // Prefix-sum to go from counts to offsts
    typename rowptrs_t::value_type nodeToElementNNZ;
    Kokkos::parallel_scan(range_policy(0, numLocalNodes + 1),
                          PrefixSumCountsFunctor(nodeToElementRowptrs), nodeToElementNNZ);
    entries_t nodeToElementEntries(Kokkos::ViewAllocateWithoutInitializing("nodeToElementEntries"), nodeToElementNNZ);
    {
      rowptrs_t insertPos(Kokkos::ViewAllocateWithoutInitializing("insertPos"), numLocalNodes + 1);
      Kokkos::deep_copy(insertPos, nodeToElementRowptrs);
      Kokkos::parallel_for(
          range_policy(0, owned_element_to_node_ids.size() + ghost_element_to_node_ids.size()),
          KOKKOS_LAMBDA(size_t i) {
            const local_ordinal_type nodesPerElement = owned_element_to_node_ids.extent(1);
            if (i < owned_element_to_node_ids.size()) {
              local_ordinal_type ownedElementIndex = i / nodesPerElement;
              local_ordinal_type nodeOfElem        = i % nodesPerElement;
              global_ordinal_type globalNode       = owned_element_to_node_ids(ownedElementIndex, nodeOfElem);
              local_ordinal_type localNode         = localDomMap.getLocalElement(globalNode);
              // Know that localNode is valid since it's adjacent to an owned element.
              if (localNode != LO_INVALID)
                nodeToElementEntries(Kokkos::atomic_fetch_add(&insertPos(localNode), size_type(1))) = ownedElementIndex;
            } else {
              i -= owned_element_to_node_ids.size();
              local_ordinal_type ghostElementIndex = i / nodesPerElement;
              local_ordinal_type nodeOfElem        = i % nodesPerElement;
              global_ordinal_type globalNode       = ghost_element_to_node_ids(ghostElementIndex, nodeOfElem);
              local_ordinal_type localNode         = localDomMap.getLocalElement(globalNode);
              if (localNode != LO_INVALID)
                nodeToElementEntries(Kokkos::atomic_fetch_add(&insertPos(localNode), size_type(1))) = owned_element_to_node_ids.extent(0) + ghostElementIndex;
            }
          });
    }
    // Using this transpose graph, can efficiently fill the desired node-node graph (packed and locally indexed).
    // For each row, use the nodeToElement graph to iterate over the adjacent elements.
    //   For each of those elements, insert the adjacent nodes' LIDs.
    //
    // Use a simple shared-memory linear probe hash table to store neighboring nodes, and a global UnorderedMap to store remaining edges if that fills up.
    // For nicely structured meshes (quad/cube, like the example) this should always fit in shared. For unstructured (e.g. tet),
    // where a node can be adjacent to many elements, we can only estimate an upper bound but shared may not be large enough.
    // Try to get power-of-2 vector length for the TeamPolicy that is at least as big as the avg #elements adjacent to each node.
    local_ordinal_type avgElementsPerNode = nodeToElementNNZ / numLocalNodes;
    int vectorLength                      = 1;
    while (vectorLength < avgElementsPerNode)
      vectorLength <<= 1;
    if (vectorLength > team_policy::vector_length_max())
      vectorLength = team_policy::vector_length_max();
    const local_ordinal_type nodesPerElement = owned_element_to_node_ids.extent(1);
    int teamSize                             = nodesPerElement;
    // Hash table size: ideally big enough to fit the median number of entries in a row (plus a few extra), but not
    // so big that occupancy suffers. For a reasonable guess, just pick 1 (for the self-loop) plus the expected average
    // entries per row including all the duplicates, but excluding self loops.
    local_ordinal_type hashTableSize = 1 + avgElementsPerNode * (nodesPerElement - 1);
    if (hashTableSize < 64)
      hashTableSize = 64;
    if (hashTableSize > 1024)
      hashTableSize = 1024;
    // First, count entries per row.
    rowptrs_t localRowptrs("localRowptrs", numLocalNodes + 1);
    // Start with a pretty small guess for the fallback table, since most rows should fit entirely in scratch
    size_t fallbackSetSize = 8 * numLocalNodes;

    while (true) {
      global_edge_set fallbackSet(fallbackSetSize);
      CountEntriesFunctor<decltype(owned_element_to_node_ids)> functor(
          localRowptrs, owned_element_to_node_ids, ghost_element_to_node_ids,
          nodeToElementRowptrs, nodeToElementEntries, fallbackSet, failFlag, hashTableSize);
      int scratchPerTeam      = scratch_hash_table::shmem_size(hashTableSize);
      team_policy dummy       = team_policy(1, Kokkos::AUTO(), vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam));
      int teamSizeRecommended = dummy.team_size_recommended(functor, Kokkos::ParallelForTag());
      if (teamSize > teamSizeRecommended)
        teamSize = teamSizeRecommended;
      Kokkos::parallel_for(team_policy(numLocalNodes, teamSize, vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam)),
                           functor);
      Kokkos::deep_copy(failFlagHost, failFlag);
      if (failFlagHost()) {
        // Need to retry whole operation with a bigger fallback set
        fallbackSetSize *= 2;
        Kokkos::deep_copy(failFlag, 0);
        Kokkos::deep_copy(localRowptrs, 0);
      } else
        break;
    }
    // Prefix-sum
    typename rowptrs_t::value_type nnz;
    Kokkos::parallel_scan(range_policy(0, numLocalNodes + 1),
                          PrefixSumCountsFunctor(localRowptrs), nnz);
    // Allocate entries
    entries_t localEntries(Kokkos::ViewAllocateWithoutInitializing("localEntries"), nnz);
    while (true) {
      global_edge_set fallbackSet(fallbackSetSize);
      FillEntriesFunctor<decltype(owned_element_to_node_ids)> functor(
          localRowptrs, localEntries, owned_element_to_node_ids, ghost_element_to_node_ids, nodeToElementRowptrs, nodeToElementEntries, column_map->getLocalMap(), fallbackSet, failFlag, hashTableSize);
      int scratchPerTeam      = scratch_hash_table::shmem_size(hashTableSize) + scratch_counter::shmem_size();
      team_policy dummy       = team_policy(1, Kokkos::AUTO(), vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam));
      int teamSizeRecommended = dummy.team_size_recommended(functor, Kokkos::ParallelForTag());
      if (teamSize > teamSizeRecommended)
        teamSize = teamSizeRecommended;
      Kokkos::parallel_for(team_policy(numLocalNodes, teamSize, vectorLength).set_scratch_size(0, Kokkos::PerTeam(scratchPerTeam)),
                           functor);
      // KokkosKernels::Impl::print_1Dview(localEntries, true);
      Kokkos::deep_copy(failFlagHost, failFlag);
      if (failFlagHost()) {
        // Need to retry whole operation with a bigger fallback set
        fallbackSetSize *= 2;
        Kokkos::deep_copy(failFlag, 0);
      } else
        break;
    }
    // The local graph is now fully constructed. Sort it and build the fill-complete CrsGraph.
    local_graph_type localGraph(localEntries, localRowptrs);
    KokkosSparse::sort_crs_graph(localGraph);
    timerElementLoopGraph = Teuchos::null;
    {
      TimeMonitor timer(*TimeMonitor::getNewTimer("2) FillComplete (Graph)"));
      crs_graph = Teuchos::rcp(new crs_graph_type(localGraph, row_map, column_map, domain_map, range_map));
    }
  }

  // Print out the crs_graph in detail...
  if (opts.verbose) crs_graph->describe(out, Teuchos::VERB_EXTREME);

  // Simulated Ghosting of Material State
  // -------------------
  {
    GhostState state(elementImporter, opts.numStateDoublesPerElement);
    TimeMonitor timer(*TimeMonitor::getNewTimer("3.1) Ghosting Material State (Matrix)"));
    state.doGhost();
  }

  // Matrix Fill
  // -------------------
  // In this example, we're using a simple stencil of values for the stiffness matrix:
  //
  //    +-----+-----+-----+-----+
  //    |  2  | -1  |     | -1  |
  //    +-----+-----+-----+-----+
  //    | -1  |  2  | -1  |     |
  //    +-----+-----+-----+-----+
  //    |     | -1  |  2  | -1  |
  //    +-----+-----+-----+-----+
  //    | -1  |     | -1  |  2  |
  //    +-----+-----+-----+-----+
  //
  // For matrix fill, we create the crs_matrix object and will fill it
  // in the same manner as we filled in the graph but in this case, nodes
  // associated with each element will receive contributions according to
  // the row in this stencil.
  //
  // In this example, the calls to sumIntoGlobalValues() on 1 core will look like:
  //   Element 0
  // - sumIntoGlobalValues( 0,  [  0  1  5  4  ],  [  2  -1  0  -1  ])
  // - sumIntoGlobalValues( 1,  [  0  1  5  4  ],  [  -1  2  -1  0  ])
  // - sumIntoGlobalValues( 5,  [  0  1  5  4  ],  [  0  -1  2  -1  ])
  // - sumIntoGlobalValues( 4,  [  0  1  5  4  ],  [  -1  0  -1  2  ])
  // Element 1
  // - sumIntoGlobalValues( 1,  [  1  2  6  5  ],  [  2  -1  0  -1  ])
  // - sumIntoGlobalValues( 2,  [  1  2  6  5  ],  [  -1  2  -1  0  ])
  // - sumIntoGlobalValues( 6,  [  1  2  6  5  ],  [  0  -1  2  -1  ])
  // - sumIntoGlobalValues( 5,  [  1  2  6  5  ],  [  -1  0  -1  2  ])
  // Element 2
  // - sumIntoGlobalValues( 2,  [  2  3  7  6  ],  [  2  -1  0  -1  ])
  // - sumIntoGlobalValues( 3,  [  2  3  7  6  ],  [  -1  2  -1  0  ])
  // - sumIntoGlobalValues( 7,  [  2  3  7  6  ],  [  0  -1  2  -1  ])
  // - sumIntoGlobalValues( 6,  [  2  3  7  6  ],  [  -1  0  -1  2  ])
  //
  // Similarly to the Graph construction above, we loop over both local and global
  // elements and insert rows for only the locally owned rows.
  //
  RCP<TimeMonitor> timerElementLoopMemory = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("3.2) ElementLoop  (Memory)")));

  RCP<crs_matrix_type> crs_matrix = rcp(new crs_matrix_type(crs_graph));
  RCP<multivector_type> rhs       = rcp(new multivector_type(crs_graph->getRowMap(), 1));
  {
    auto owned_element_to_node_ids = mesh.getOwnedElementToNode().getDeviceView(Tpetra::Access::ReadOnly);
    auto ghost_element_to_node_ids = mesh.getGhostElementToNode().getDeviceView(Tpetra::Access::ReadOnly);
    auto localMatrix               = crs_matrix->getLocalMatrixDevice();
    auto localRHS                  = rhs->getLocalViewDevice(Tpetra::Access::OverwriteAll);
    auto localRowMap               = crs_matrix->getRowMap()->getLocalMap();
    auto localColMap               = crs_matrix->getColMap()->getLocalMap();

    // Because we're processing elements in parallel, we need storage for all of them
    int numOwnedElements = mesh.getNumOwnedElements();
    int numGhostElements = mesh.getNumGhostElements();
    int nperel           = owned_element_to_node_ids.extent(1);
    local_ordinal_single_view_type all_lcids("all_lids", nperel * std::max(numOwnedElements, numGhostElements));

    timerElementLoopMemory                  = Teuchos::null;
    RCP<TimeMonitor> timerElementLoopMatrix = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("3.3) ElementLoop  (Matrix)")));

    // Loop over owned elements:
    Kokkos::parallel_for(
        Kokkos::RangePolicy<execution_space>(0, numOwnedElements), KOKKOS_LAMBDA(const size_t& element_gidx) {
          // Get subviews
          pair_type location_pair = pair_type(nperel * element_gidx, nperel * (element_gidx + 1));
          auto element_lcids      = Kokkos::subview(all_lcids, location_pair);

          // Get the contributions for the current element
          Scalar element_matrix[4][4];
          Scalar element_rhs[4];
          ReferenceQuad4(element_matrix, element_rhs);

          // Get the local column ids array for this element
          for (int element_node_idx = 0; element_node_idx < nperel; element_node_idx++) {
            element_lcids(element_node_idx) = localColMap.getLocalElement(owned_element_to_node_ids(element_gidx, element_node_idx));
          }

          // For each node (row) on the current element:
          // - populate the values array
          // - add the values to the fe_matrix.
          for (int element_node_idx = 0; element_node_idx < nperel; element_node_idx++) {
            global_ordinal_type global_row_id = owned_element_to_node_ids(element_gidx, element_node_idx);
            local_ordinal_type local_row_id   = localRowMap.getLocalElement(global_row_id);
            if (local_row_id != LO_INVALID) {
              // Force atomics on sums
              for (int col_idx = 0; col_idx < nperel; col_idx++)
                localMatrix.sumIntoValues(local_row_id, &element_lcids(col_idx), 1, &(element_matrix[element_node_idx][col_idx]), true, true);
              Kokkos::atomic_add(&(localRHS(local_row_id, 0)), element_rhs[element_node_idx]);
            }
          }
        });

    // Loop over ghost elements:
    // - This loop is the same as the element loop for owned elements, but this one
    //   is for ghost elements.
    Kokkos::parallel_for(
        Kokkos::RangePolicy<execution_space>(0, numGhostElements), KOKKOS_LAMBDA(const size_t& element_gidx) {
          // Get subviews
          pair_type location_pair = pair_type(nperel * element_gidx, nperel * (element_gidx + 1));
          auto element_lcids      = Kokkos::subview(all_lcids, location_pair);

          // Get the contributions for the current element
          Scalar element_matrix[4][4];
          Scalar element_rhs[4];
          ReferenceQuad4(element_matrix, element_rhs);

          // Get the local column ids array for this element
          for (int element_node_idx = 0; element_node_idx < nperel; element_node_idx++) {
            element_lcids(element_node_idx) = localColMap.getLocalElement(ghost_element_to_node_ids(element_gidx, element_node_idx));
          }

          for (int element_node_idx = 0; element_node_idx < nperel; element_node_idx++) {
            global_ordinal_type global_row_id = ghost_element_to_node_ids(element_gidx, element_node_idx);
            local_ordinal_type local_row_id   = localRowMap.getLocalElement(global_row_id);
            if (local_row_id != LO_INVALID) {
              // Force atomics on sums
              for (int col_idx = 0; col_idx < nperel; col_idx++)
                localMatrix.sumIntoValues(local_row_id, &element_lcids(col_idx), 1, &(element_matrix[element_node_idx][col_idx]), true, true);
              Kokkos::atomic_add(&(localRHS(local_row_id, 0)), element_rhs[element_node_idx]);
            }
          }
        });
    execution_space().fence();
    timerElementLoopMatrix = Teuchos::null;
  }

  // After the contributions are added, 'finalize' the matrix using fillComplete()
  {
    TimeMonitor timer(*TimeMonitor::getNewTimer("4) FillComplete (Matrix)"));
    crs_matrix->fillComplete();
  }

  // Print out crs_matrix details.
  if (opts.verbose) crs_matrix->describe(out, Teuchos::VERB_EXTREME);

  Teuchos::TimeMonitor::getStackedTimer()->stopBaseTimer();

  // Save crs_matrix as a MatrixMarket file.
  if (opts.saveMM) {
    std::ofstream ofs("crsMatrix_TotalElementLoop_SPKokkos.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparse(ofs, crs_matrix);
    std::ofstream ofs2("rhs_TotalElementLoop_SPKokkos.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<multivector_type>::writeDense(ofs2, rhs);
  }

  return 0;
}

}  // namespace TpetraExamples

#endif  // TPETRAEXAMPLES_FEM_ASSEMBLY_TOTALELEMENTLOOP_HPP
