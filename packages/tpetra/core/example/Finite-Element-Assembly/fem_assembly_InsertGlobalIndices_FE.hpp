// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_INSERTGLOBALINDICES_FE_SP_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_INSERTGLOBALINDICES_FE_SP_HPP

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "Tpetra_Core.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Tpetra_Assembly_Helpers.hpp"
#include "Tpetra_GraphAssembly.hpp"

#include "fem_assembly_typedefs.hpp"
#include "fem_assembly_MeshDatabase.hpp"
#include "fem_assembly_Element.hpp"
#include "fem_assembly_utility.hpp"
#include "fem_assembly_commandLineOpts.hpp"

namespace TpetraExamples {

int executeInsertGlobalIndicesFESP_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const struct CmdLineOpts& opts);

int executeInsertGlobalIndicesFESPKokkos_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                          const struct CmdLineOpts& opts);

int executeInsertGlobalIndicesFESP(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                   const struct CmdLineOpts& opts) {
  using Teuchos::RCP;

  // The output stream 'out' will ignore any output not from Process 0.
  RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
  Teuchos::FancyOStream& out      = *pOut;

  std::string useKokkos = opts.useKokkosAssembly ? "Kokkos Assembly" : "Serial Assembly";
  out << "================================================================================" << std::endl
      << "=  FECrsMatrix Insert Global Indices (Static Profile; " << useKokkos << ")" << std::endl
      << "================================================================================" << std::endl
      << std::endl;

  int status = 0;

  for (size_t i = 0; i < opts.repetitions; ++i) {
    if (opts.useKokkosAssembly)
      status += executeInsertGlobalIndicesFESPKokkos_(comm, opts);
    else
      status += executeInsertGlobalIndicesFESP_(comm, opts);
  }
  return status;
}

int executeInsertGlobalIndicesFESP_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
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
  constexpr int MAX_NODES_PER_ELEM = 4;
  const int nodesPerElem           = (int)mesh.getNodesPerElem();
  assert(nodesPerElem <= MAX_NODES_PER_ELEM);
  if (opts.verbose) mesh.print(std::cout);

  // Build Tpetra Maps
  // -----------------
  // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
  RCP<const map_type> row_map =
      rcp(new map_type(GO_INVALID, mesh.getOwnedNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly),
                       0, comm));
  RCP<const map_type> owned_plus_shared_map =
      rcp(new map_type(GO_INVALID, mesh.getOwnedAndGhostNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly),
                       0, comm));

  if (opts.verbose) row_map->describe(out);

  // Graph Construction
  // ------------------
  // - Loop over every element in the mesh.
  //   - Get list of nodes associated with each element.
  //   - Insert the clique of nodes associated with each element into the graph.
  //
  auto domain_map = row_map;
  auto range_map  = row_map;

  auto owned_element_to_node_gids = mesh.getOwnedElementToNode().getHostView(Tpetra::Access::ReadOnly);

  Teuchos::TimeMonitor::getStackedTimer()->startBaseTimer();
  RCP<TimeMonitor> timerElementLoopGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("1) ElementLoop  (Graph)")));

  RCP<fe_graph_type> fe_graph = rcp(new fe_graph_type(row_map, owned_plus_shared_map, nodesPerElem * nodesPerElem));

  // Because we're using quads for this example, there will be 4 nodes associated with each element.
  Teuchos::Array<global_ordinal_type> global_ids_in_row(nodesPerElem);

  // for each element in the mesh...
  Tpetra::beginAssembly(*fe_graph);
  for (size_t element_gidx = 0; element_gidx < mesh.getNumOwnedElements(); element_gidx++) {
    // Populate global_ids_in_row:
    // - Copy the global node ids for current element into an array.
    // - Since each element's contribution is a clique, we can re-use this for
    //   each row associated with this element's contribution.
    for (int element_node_idx = 0; element_node_idx < nodesPerElem; element_node_idx++) {
      global_ids_in_row[element_node_idx] =
          owned_element_to_node_gids(element_gidx, element_node_idx);
    }

    // Add the contributions from the current row into the graph.
    // - For example, if Element 0 contains nodes [0,1,4,5] then we insert the nodes:
    //   - node 0 inserts [0, 1, 4, 5]
    //   - node 1 inserts [0, 1, 4, 5]
    //   - node 4 inserts [0, 1, 4, 5]
    //   - node 5 inserts [0, 1, 4, 5]
    for (int element_node_idx = 0; element_node_idx < nodesPerElem; element_node_idx++) {
      fe_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
    }
  }
  timerElementLoopGraph = Teuchos::null;

  // Call fillComplete on the fe_graph to 'finalize' it.
  {
    TimeMonitor timer(*TimeMonitor::getNewTimer("2) FillComplete (Graph)"));
    Tpetra::endAssembly(*fe_graph);
  }

  // Print out verbose information about the fe_graph.
  if (opts.verbose) fe_graph->describe(out, Teuchos::VERB_EXTREME);

  // Matrix Fill
  // -------------------
  // In this example, we're using a simple stencil of values for the matrix fill:
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
  // For Type 1 matrix fill, we create the fe_matrix object and will fill it
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

  const int numOwnedElements = mesh.getNumOwnedElements();
  RCP<fe_matrix_type> fe_matrix;
  RCP<fe_multivector_type> rhs;

  {
    TimeMonitor timerElementLoopMatrix(*TimeMonitor::getNewTimer("3) ElementLoop  (Matrix)"));

    fe_matrix = rcp(new fe_matrix_type(fe_graph));
    rhs       = rcp(new fe_multivector_type(domain_map, fe_graph->getImporter(), 1));

    Scalar element_matrix[MAX_NODES_PER_ELEM][MAX_NODES_PER_ELEM];
    Scalar element_rhs[MAX_NODES_PER_ELEM];

    Teuchos::Array<global_ordinal_type> column_global_ids(nodesPerElem);  // global column ids list
    Teuchos::Array<Scalar> column_scalar_values(nodesPerElem);            // scalar values for each column

    // Loop over elements
    Tpetra::beginAssembly(*fe_matrix, *rhs);
    for (int element_gidx = 0; element_gidx < numOwnedElements;
         ++element_gidx) {
      // Get the contributions for the current element
      // shape info injected here in real life
      // "GetElementMatrix"
      ReferenceQuad4(element_matrix, element_rhs);

      for (int element_node_idx = 0;
           element_node_idx < nodesPerElem;
           ++element_node_idx) {
        column_global_ids[element_node_idx] =
            owned_element_to_node_gids(element_gidx, element_node_idx);
      }

      // For each node (row) on the current element:
      // - populate the values array
      // - add the values to the fe_matrix.
      for (int element_node_idx = 0; element_node_idx < nodesPerElem; ++element_node_idx) {
        global_ordinal_type global_row_id =
            owned_element_to_node_gids(element_gidx, element_node_idx);

        for (int col_idx = 0; col_idx < nodesPerElem; col_idx++) {
          column_scalar_values[col_idx] =
              element_matrix[element_node_idx][col_idx];
        }

        fe_matrix->sumIntoGlobalValues(
            global_row_id, column_global_ids, column_scalar_values);
        rhs->sumIntoGlobalValue(
            global_row_id, 0, element_rhs[element_node_idx]);
      }
    }
  }  // timerElementLoopMatrix

  // After the contributions are added, 'finalize' the matrix using fillComplete()
  {
    TimeMonitor timer(*TimeMonitor::getNewTimer("4) FillComplete (Matrix)"));
    Tpetra::endAssembly(*fe_matrix);
  }

  {
    // Global assemble the RHS
    TimeMonitor timer(*TimeMonitor::getNewTimer("5) GlobalAssemble (RHS)"));
    Tpetra::endAssembly(*rhs);
  }

  Teuchos::TimeMonitor::getStackedTimer()->stopBaseTimer();

  // Print out fe_matrix details.
  if (opts.verbose) fe_matrix->describe(out, Teuchos::VERB_EXTREME);

  // Save crs_matrix as a MatrixMarket file.
  if (opts.saveMM) {
    std::ofstream ofs("crsMatrix_InsertGlobalIndices_FESP.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparse(ofs, fe_matrix);
    std::ofstream ofs2("rhs_InsertGlobalIndices_FESP.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDense(ofs2, rhs);
  }

  return 0;
}

int executeInsertGlobalIndicesFESPKokkos_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
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

  // We're going to allocate device register storage for the matrix/rhs
  // at compile time, so we need a upper bound on the # nodes per element
  constexpr int MAX_NODES_PER_ELEM = 4;
  const int nodesPerElem           = (int)mesh.getNodesPerElem();
  assert(nodesPerElem <= MAX_NODES_PER_ELEM);

  if (opts.verbose) mesh.print(std::cout);

  // Build Tpetra Maps
  // -----------------
  // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
  RCP<const map_type> row_map =
      rcp(new map_type(GO_INVALID, mesh.getOwnedNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly),
                       0, comm));
  RCP<const map_type> owned_plus_shared_map =
      rcp(new map_type(GO_INVALID, mesh.getOwnedAndGhostNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly),
                       0, comm));

  if (opts.verbose) {
    row_map->describe(out);
  }

  // Graph Construction
  // ------------------
  // Delegate the entire graph assembly to Tpetra::Experimental::GraphAssembly.
  // Given the owned (row) map, the owned+shared map, and the owned
  // element-to-node connectivity, it builds the local sparsity pattern on
  // device with Kokkos kernels, wraps it in a packed CrsGraph over the
  // owned+shared map, and does the owned+shared -> owned communication with a
  // fused Export + fillComplete.
  //
  // Besides the final owned graph, we also grab the intermediate owned+shared
  // "assembly" graph and the Export object it built.  We reuse both below to
  // assemble the matrix values / RHS on the owned+shared map and then Export
  // them to the owned objects with the same communication pattern.
  //
  auto domain_map = row_map;
  auto range_map  = row_map;

  Teuchos::TimeMonitor::getStackedTimer()->startBaseTimer();

  using graph_assembly_type =
      Tpetra::Experimental::GraphAssembly<local_ordinal_type, global_ordinal_type,
                                          map_type::node_type>;

  RCP<crs_graph_type> crs_graph;
  // The locally-assembled (owned+shared) graph and the Export used to move data
  // from it to the owned graph.  Both are reused for the matrix/RHS fill below.
  RCP<crs_graph_type> owned_plus_shared_graph;
  RCP<const export_type> graph_exporter;
  {
    RCP<TimeMonitor> timerElementLoopGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("1) ElementLoop  (Graph)")));

    // The element-to-node connectivity (global node IDs) of the owned elements,
    // as a dynamic-rank-2 device view (what GraphAssembly expects).
    graph_assembly_type::element_to_node_type owned_element_to_node_ids =
        mesh.getOwnedElementToNode().getDeviceView(Tpetra::Access::ReadOnly);

    graph_assembly_type assembler(row_map, owned_plus_shared_map, owned_element_to_node_ids);

    {
      TimeMonitor timer(*TimeMonitor::getNewTimer("2) FillComplete (Graph)"));
      assembler.build();
    }
    timerElementLoopGraph = Teuchos::null;

    crs_graph               = assembler.getGraph();
    owned_plus_shared_graph = assembler.getOwnedPlusSharedGraph();
    graph_exporter          = assembler.getExporter();
  }

  // Print out verbose information about the graph.
  if (opts.verbose) crs_graph->describe(out, Teuchos::VERB_EXTREME);

  // Matrix Fill
  // -------------------
  // In this example, we're using a simple stencil of values for the matrix fill:
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
  // For Type 1 matrix fill, we create the fe_matrix object and will fill it
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
  RCP<TimeMonitor> timerElementLoopMemory = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("3.1) ElementLoop  (Memory)")));

  // the number of finite elements this process will handle, and the number of
  // nodes associated with each element.
  // This information will be used to pre-allocate storage to process elements
  // in parallel.
  const int numOwnedElements = mesh.getNumOwnedElements();

  // The final (owned) matrix and RHS.
  RCP<crs_matrix_type> crs_matrix = rcp(new crs_matrix_type(crs_graph));
  RCP<multivector_type> rhs       = rcp(new multivector_type(domain_map, 1));

  // The assembly matrix/RHS live on the owned+shared map and share the packed
  // sparsity pattern that the Kokkos kernels above produced.  Assembling into
  // these means every element contribution is a purely local, atomic
  // sumIntoValues -- no "nonowned entry" bookkeeping at all.
  RCP<crs_matrix_type> owned_plus_shared_matrix =
      rcp(new crs_matrix_type(owned_plus_shared_graph));
  RCP<multivector_type> owned_plus_shared_rhs =
      rcp(new multivector_type(owned_plus_shared_map, 1));

  // Row map and column map of the assembly objects are both the owned+shared
  // map, so a single local map handles both row and column lookups.
  auto localMap = owned_plus_shared_map->getLocalMap();

  timerElementLoopMemory = Teuchos::null;
  {
    TimeMonitor timerElementLoopMatrix(*TimeMonitor::getNewTimer("3.2) ElementLoop  (Matrix)"));

    auto owned_element_to_node_gids =
        mesh.getOwnedElementToNode().getDeviceView(
            Tpetra::Access::ReadOnly);

    // Loop over elements
    auto localRHS =
        owned_plus_shared_rhs->getLocalViewDevice(Tpetra::Access::OverwriteAll);
    auto localMatrix = owned_plus_shared_matrix->getLocalMatrixDevice();

    // We're not doing explicit worksetting, but if we did, we'd do it here
    Kokkos::parallel_for(
        "Assemble FE matrix and right-hand side",
        Kokkos::RangePolicy<execution_space, int>(0, numOwnedElements),
        KOKKOS_LAMBDA(const size_t element_idx) {
          Scalar element_matrix[MAX_NODES_PER_ELEM][MAX_NODES_PER_ELEM];
          local_ordinal_type element_lcids[MAX_NODES_PER_ELEM];
          Scalar element_rhs[MAX_NODES_PER_ELEM];

          // Get the contributions for the current element
          ReferenceQuad4(element_matrix, element_rhs);

          // Get the local column ids array for this element
          for (int element_node_idx = 0; element_node_idx < nodesPerElem; ++element_node_idx) {
            element_lcids[element_node_idx] =
                localMap.getLocalElement(owned_element_to_node_gids(element_idx, element_node_idx));
          }

          // For each node (row) on the current element:
          // - populate the values array
          // - add the values to the assembly matrix.
          for (int element_node_idx = 0; element_node_idx < nodesPerElem; ++element_node_idx) {
            const local_ordinal_type local_row_id = element_lcids[element_node_idx];

            // Atomically contribute for sums: parallel elements may be contributing
            // to the same node at the same time
            for (int col_idx = 0; col_idx < nodesPerElem; ++col_idx) {
              localMatrix.sumIntoValues(local_row_id, &element_lcids[col_idx], 1,
                                        &(element_matrix[element_node_idx][col_idx]),
                                        true, true);
            }
            Kokkos::atomic_add(&(localRHS(local_row_id, 0)), element_rhs[element_node_idx]);
          }
        });
    execution_space().fence();
  }

  // Communicate the matrix values explicitly: Export owned+shared -> owned with
  // ADD.  The target matrix has a static (already fill-complete) graph, so every
  // incoming entry lands via sumIntoGlobalValues on an existing position -- no
  // structure changes, no reallocation.
  {
    TimeMonitor timer(*TimeMonitor::getNewTimer("4) FillComplete (Matrix)"));
    if (graph_exporter.is_null()) {
      // Serial / single-rank: the assembly matrix already IS the answer, so just
      // copy the values across (same graph, same ordering).
      Kokkos::deep_copy(crs_matrix->getLocalMatrixDevice().values,
                        owned_plus_shared_matrix->getLocalMatrixDevice().values);
    } else {
      crs_matrix->doExport(*owned_plus_shared_matrix, *graph_exporter, Tpetra::ADD);
    }
    crs_matrix->fillComplete(domain_map, range_map);
  }

  {
    // Communicate the RHS the same way: Export owned+shared -> owned with ADD.
    TimeMonitor timer(*TimeMonitor::getNewTimer("5) GlobalAssemble (RHS)"));
    if (graph_exporter.is_null()) {
      Kokkos::deep_copy(rhs->getLocalViewDevice(Tpetra::Access::OverwriteAll),
                        owned_plus_shared_rhs->getLocalViewDevice(Tpetra::Access::ReadOnly));
    } else {
      rhs->doExport(*owned_plus_shared_rhs, *graph_exporter, Tpetra::ADD);
    }
  }

  Teuchos::TimeMonitor::getStackedTimer()->stopBaseTimer();

  // Print out crs_matrix details.
  if (opts.verbose) crs_matrix->describe(out, Teuchos::VERB_EXTREME);

  // Save crs_matrix as a MatrixMarket file.
  if (opts.saveMM) {
    std::ofstream ofs("crsMatrix_InsertGlobalIndices_FESPKokkos.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparse(ofs, crs_matrix);
    std::ofstream ofs2("rhs_InsertGlobalIndices_FESPKokkos.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDense(ofs2, rhs);
  }

  return 0;
}

}  // namespace TpetraExamples

#endif  // TPETRAEXAMPLES_FEM_ASSEMBLY_INSERTGLOBALINDICES_FE_SP_HPP
