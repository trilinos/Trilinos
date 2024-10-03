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

#include "fem_assembly_typedefs.hpp"
#include "fem_assembly_MeshDatabase.hpp"
#include "fem_assembly_Element.hpp"
#include "fem_assembly_utility.hpp"
#include "fem_assembly_commandLineOpts.hpp"


namespace TpetraExamples
{

int executeInsertGlobalIndicesFESP_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const struct CmdLineOpts& opts);

int executeInsertGlobalIndicesFESPKokkos_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const struct CmdLineOpts& opts);



int executeInsertGlobalIndicesFESP(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                   const struct CmdLineOpts & opts)
{
  using Teuchos::RCP;

  // The output stream 'out' will ignore any output not from Process 0.
  RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
  Teuchos::FancyOStream& out = *pOut;

  std::string useKokkos = opts.useKokkosAssembly ? "Kokkos Assembly" : "Serial Assembly";
  out << "================================================================================" << std::endl
      << "=  FECrsMatrix Insert Global Indices (Static Profile; "<<useKokkos <<")"    << std::endl
      << "================================================================================" << std::endl
      << std::endl;

  int status = 0;

  for(size_t i=0; i<opts.repetitions; ++i)
  {

    if(opts.useKokkosAssembly)
      status += executeInsertGlobalIndicesFESPKokkos_(comm, opts);
    else
      status += executeInsertGlobalIndicesFESP_(comm, opts);

  }
  return status;
}

int executeInsertGlobalIndicesFESP_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const struct CmdLineOpts& opts)
{
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  const global_ordinal_type GO_INVALID = Teuchos::OrdinalTraits<global_ordinal_type>::invalid();


  // The output stream 'out' will ignore any output not from Process 0.
  RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
  Teuchos::FancyOStream& out = *pOut;

  // Processor decomp (only works on perfect squares)
  int numProcs  = comm->getSize();
  int sqrtProcs = sqrt(numProcs);

  if(sqrtProcs*sqrtProcs != numProcs)
  {
    if(0 == comm->getRank())
      std::cerr << "Error: Invalid number of processors provided, num processors must be a perfect square." << std::endl;
    return -1;
  }
  int procx = sqrtProcs;
  int procy = sqrtProcs;

  // Generate a simple 3x3 mesh
  int nex = opts.numElementsX;
  int ney = opts.numElementsY;

  MeshDatabase mesh(comm,nex,ney,procx,procy);

  if(opts.verbose) mesh.print(std::cout);

  // Build Tpetra Maps
  // -----------------
  // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
  RCP<const map_type> row_map =
    rcp(new map_type(GO_INVALID, mesh.getOwnedNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly),
                     0, comm));
  RCP<const map_type> owned_plus_shared_map =
    rcp(new map_type(GO_INVALID, mesh.getOwnedAndGhostNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly),
                     0, comm));

  if(opts.verbose) row_map->describe(out);

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

  RCP<fe_graph_type> fe_graph = rcp(new fe_graph_type(row_map, owned_plus_shared_map, 16));

  // Using 4 because we're using quads for this example, so there will be 4 nodes associated with each element.
  Teuchos::Array<global_ordinal_type> global_ids_in_row(4);

  // for each element in the mesh...
  Tpetra::beginAssembly(*fe_graph);
  for(size_t element_gidx=0; element_gidx<mesh.getNumOwnedElements(); element_gidx++)
  {
    // Populate global_ids_in_row:
    // - Copy the global node ids for current element into an array.
    // - Since each element's contribution is a clique, we can re-use this for
    //   each row associated with this element's contribution.
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_gids.extent(1); element_node_idx++)
    {
      global_ids_in_row[element_node_idx] = 
	owned_element_to_node_gids(element_gidx, element_node_idx);
    }

    // Add the contributions from the current row into the graph.
    // - For example, if Element 0 contains nodes [0,1,4,5] then we insert the nodes:
    //   - node 0 inserts [0, 1, 4, 5]
    //   - node 1 inserts [0, 1, 4, 5]
    //   - node 4 inserts [0, 1, 4, 5]
    //   - node 5 inserts [0, 1, 4, 5]
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_gids.extent(1); element_node_idx++)
    {
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
  if(opts.verbose) fe_graph->describe(out, Teuchos::VERB_EXTREME);

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
  const int nodesPerElem = mesh.getOwnedElementToNode().extent(1);

  RCP<fe_matrix_type> fe_matrix;
  RCP<fe_multivector_type> rhs;

  {
    TimeMonitor timerElementLoopMatrix
      (*TimeMonitor::getNewTimer ("3) ElementLoop  (Matrix)"));

    fe_matrix = rcp(new fe_matrix_type(fe_graph));
    rhs = rcp (new fe_multivector_type(domain_map, fe_graph->getImporter(), 1));

    Kokkos::View<local_ordinal_type[4][4], hostType>
      element_matrix("element_matrix");
    Teuchos::Array<Scalar> element_rhs(4);

    Teuchos::Array<global_ordinal_type> column_global_ids(4);     // global column ids list
    Teuchos::Array<Scalar> column_scalar_values(4);         // scalar values for each column

    // Loop over elements
    Tpetra::beginAssembly(*fe_matrix,*rhs);
    for (int element_gidx = 0; element_gidx < numOwnedElements;
         ++element_gidx) {
      // Get the contributions for the current element
      // shape info injected here in real life
      // "GetElementMatrix"
      ReferenceQuad4(element_matrix);
      ReferenceQuad4RHS(element_rhs);

      for (int element_node_idx=0;
	         element_node_idx < nodesPerElem;
	         ++element_node_idx) {
	      column_global_ids[element_node_idx] =
	        owned_element_to_node_gids(element_gidx, element_node_idx);
      }

      // For each node (row) on the current element:
      // - populate the values array
      // - add the values to the fe_matrix.
      // Note: hardcoded 4 here because we're using quads.
      for (int element_node_idx = 0; element_node_idx < 4;
	         ++element_node_idx) {
	      global_ordinal_type global_row_id =
	        owned_element_to_node_gids(element_gidx, element_node_idx);
	
        for(int col_idx = 0; col_idx < 4; col_idx++) {
          column_scalar_values[col_idx] = 
            element_matrix(element_node_idx, col_idx);
        }

        fe_matrix->sumIntoGlobalValues(
          global_row_id, column_global_ids, column_scalar_values);
        rhs->sumIntoGlobalValue(
          global_row_id, 0, element_rhs[element_node_idx]);
      }
    }
  } // timerElementLoopMatrix

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
  if(opts.verbose) fe_matrix->describe(out, Teuchos::VERB_EXTREME);

  // Save crs_matrix as a MatrixMarket file.
  if(opts.saveMM)
  {
    std::ofstream ofs("crsMatrix_InsertGlobalIndices_FESP.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparse(ofs, fe_matrix);
    std::ofstream ofs2("rhs_InsertGlobalIndices_FESP.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<multivector_type>::writeDense(ofs2, rhs);
  }

  return 0;
}


int executeInsertGlobalIndicesFESPKokkos_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const struct CmdLineOpts& opts)
{
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  const global_ordinal_type GO_INVALID = Teuchos::OrdinalTraits<global_ordinal_type>::invalid();
  using pair_type = Kokkos::pair<int,int>;

  // The output stream 'out' will ignore any output not from Process 0.
  RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
  Teuchos::FancyOStream& out = *pOut;

  // Processor decomp (only works on perfect squares)
  int numProcs  = comm->getSize();
  int sqrtProcs = sqrt(numProcs);

  if(sqrtProcs*sqrtProcs != numProcs)
  {
    if(0 == comm->getRank())
      std::cerr << "Error: Invalid number of processors provided, num processors must be a perfect square." << std::endl;
    return -1;
  }
  int procx = sqrtProcs;
  int procy = sqrtProcs;

  // Generate a simple 3x3 mesh
  int nex = opts.numElementsX;
  int ney = opts.numElementsY;

  MeshDatabase mesh(comm,nex,ney,procx,procy);

  if(opts.verbose) mesh.print(std::cout);

  // Build Tpetra Maps
  // -----------------
  // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
  RCP<const map_type> row_map =
    rcp (new map_type (GO_INVALID, mesh.getOwnedNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly),
                       0, comm));
  RCP<const map_type> owned_plus_shared_map =
    rcp (new map_type (GO_INVALID, mesh.getOwnedAndGhostNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly),
                       0, comm));

  if (opts.verbose) {
    row_map->describe (out);
  }

  // Graph Construction
  // ------------------
  // - Loop over every element in the mesh.
  //   - Get list of nodes associated with each element.
  //   - Insert the clique of nodes associated with each element into the graph.
  //
  auto domain_map = row_map;
  auto range_map  = row_map;


  Teuchos::TimeMonitor::getStackedTimer()->startBaseTimer();

  RCP<fe_graph_type> fe_graph =
    rcp (new fe_graph_type (row_map, owned_plus_shared_map, 16));

  // Using 4 because we're using quads for this example, so there will
  // be 4 nodes associated with each element.
  Teuchos::Array<global_ordinal_type> global_ids_in_row(4);

  {
    TimeMonitor timerElementLoopGraph(*TimeMonitor::getNewTimer("1) ElementLoop  (Graph)"));
    auto owned_element_to_node_gids = mesh.getOwnedElementToNode().getHostView(Tpetra::Access::ReadOnly);

    // for each element in the mesh...
    Tpetra::beginAssembly(*fe_graph);
    for (size_t element_gidx = 0;
         element_gidx < mesh.getNumOwnedElements ();
         ++element_gidx) {
      // Populate global_ids_in_row:
      // - Copy the global node ids for current element into an array.
      // - Since each element's contribution is a clique, we can re-use this for
      //   each row associated with this element's contribution.
      for (size_t element_node_idx = 0;
           element_node_idx < owned_element_to_node_gids.extent(1);
           ++element_node_idx) {
        global_ids_in_row[element_node_idx] =
          owned_element_to_node_gids(element_gidx, element_node_idx);
      }

      // Add the contributions from the current row into the graph.
      // - For example, if Element 0 contains nodes [0,1,4,5] then we insert the nodes:
      //   - node 0 inserts [0, 1, 4, 5]
      //   - node 1 inserts [0, 1, 4, 5]
      //   - node 4 inserts [0, 1, 4, 5]
      //   - node 5 inserts [0, 1, 4, 5]
      for (size_t element_node_idx = 0;
           element_node_idx < owned_element_to_node_gids.extent(1); ++element_node_idx) {
        fe_graph->insertGlobalIndices (global_ids_in_row[element_node_idx],
                                       global_ids_in_row());
      }
    }
  }


  // Call fillComplete on the fe_graph to 'finalize' it.
  {
    TimeMonitor timer(*TimeMonitor::getNewTimer("2) FillComplete (Graph)"));
    Tpetra::endAssembly(*fe_graph);
  }


  // Print out verbose information about the fe_graph.
  if(opts.verbose) fe_graph->describe(out, Teuchos::VERB_EXTREME);

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
  const int nodesPerElem = mesh.getOwnedElementToNode().extent(1);

  RCP<fe_matrix_type> fe_matrix = rcp(new fe_matrix_type(fe_graph));
  RCP<fe_multivector_type> rhs =
    rcp (new fe_multivector_type(domain_map, fe_graph->getImporter(), 1));


  
  auto localMap     = owned_plus_shared_map->getLocalMap();
  auto localColMap  = fe_matrix->getColMap()->getLocalMap();
 
  // no worksetting in this example
  scalar_2d_array_type all_element_matrix(
    "all_element_matrix",nodesPerElem*numOwnedElements);
  scalar_1d_array_type all_element_rhs(
    "all_element_rhs",nodesPerElem*numOwnedElements);
  local_ordinal_single_view_type  all_lcids(
    "all_lids",nodesPerElem*numOwnedElements);

  timerElementLoopMemory=Teuchos::null;
  {
    TimeMonitor timerElementLoopMatrix(*TimeMonitor::getNewTimer ("3.2) ElementLoop  (Matrix)"));
    Tpetra::beginAssembly(*fe_matrix,*rhs);

    auto owned_element_to_node_gids = 
      mesh.getOwnedElementToNode().getDeviceView(
        Tpetra::Access::ReadOnly);

    // Loop over elements
    auto localRHS = 
      rhs->getLocalViewDevice(Tpetra::Access::OverwriteAll);
    auto localMatrix = fe_matrix->getLocalMatrixDevice();
    // Work around subview of managed views being slower than unmanaged
    auto all_element_rhs_unmanaged = makeUnmanaged(all_element_rhs);
    auto all_element_matrix_unmanaged = makeUnmanaged(all_element_matrix);
    auto all_lcids_unmanaged = makeUnmanaged(all_lcids);
    Kokkos::parallel_for
      ("Assemble FE matrix and right-hand side",
       Kokkos::RangePolicy<execution_space, int> (0, numOwnedElements),
       KOKKOS_LAMBDA (const size_t element_idx) {
      const pair_type location_pair(
        nodesPerElem*element_idx, nodesPerElem*(element_idx+1));

      // this thread's piece of device allocation
      auto element_matrix = Kokkos::subview(
        all_element_matrix_unmanaged, location_pair, Kokkos::ALL);
      auto element_lcids = 
        Kokkos::subview(all_lcids_unmanaged,location_pair);
      auto element_rhs = 
        Kokkos::subview(all_element_rhs_unmanaged, location_pair);

      // Get the contributions for the current element
      ReferenceQuad4(element_matrix);
      ReferenceQuad4RHS(element_rhs);

      // Get the local column ids array for this element
      for (int element_node_idx = 0; element_node_idx < nodesPerElem;
            ++element_node_idx) {
        element_lcids(element_node_idx) =
          localColMap.getLocalElement(
            owned_element_to_node_gids(element_idx, element_node_idx));
      }
        
      // For each node (row) on the current element:
      // - populate the values array
      // - add the values to the fe_matrix.
      for (int element_node_idx = 0; element_node_idx < nodesPerElem; 
            ++element_node_idx) {
        const local_ordinal_type local_row_id =
          localMap.getLocalElement(owned_element_to_node_gids(
            element_idx, element_node_idx));

        // Atomically contribute for sums: parallel elements may be contributing
        // to the same node at the same time
        for (int col_idx = 0; col_idx < nodesPerElem; ++col_idx) {
          localMatrix.sumIntoValues(local_row_id, &element_lcids(col_idx), 1,
                                    &(element_matrix(element_node_idx,col_idx)),
                                    true, true);
        }
        Kokkos::atomic_add(
          &(localRHS(local_row_id,0)), element_rhs[element_node_idx]);
      }
    });
  }


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
  if(opts.verbose) fe_matrix->describe(out, Teuchos::VERB_EXTREME);

  // Save crs_matrix as a MatrixMarket file.
  if(opts.saveMM)
  {
    std::ofstream ofs("crsMatrix_InsertGlobalIndices_FESPKokkos.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparse(ofs, fe_matrix);
    std::ofstream ofs2("rhs_InsertGlobalIndices_FESPKokkos.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<multivector_type>::writeDense(ofs2, rhs);
  }

  return 0;
}

} // namespace TpetraExamples

#endif // TPETRAEXAMPLES_FEM_ASSEMBLY_INSERTGLOBALINDICES_FE_SP_HPP
