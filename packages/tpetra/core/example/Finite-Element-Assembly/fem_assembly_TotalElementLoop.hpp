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

#include "fem_assembly_typedefs.hpp"
#include "fem_assembly_MeshDatabase.hpp"
#include "fem_assembly_Element.hpp"
#include "fem_assembly_utility.hpp"
#include "fem_assembly_commandLineOpts.hpp"

namespace TpetraExamples
{

  int executeTotalElementLoopSP_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                 const struct CmdLineOpts& opts);

  int executeTotalElementLoopSPKokkos_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                       const struct CmdLineOpts& opts);

  int executeTotalElementLoopSP(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                const struct CmdLineOpts & opts)
  {
    using Teuchos::RCP;

    // The output stream 'out' will ignore any output not from Process 0.
    RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
    Teuchos::FancyOStream& out = *pOut;

    std::string useKokkos = opts.useKokkosAssembly ? "Kokkos Assembly" : "Serial Assembly";
    out << "================================================================================" << std::endl
        << "=  Total Element Loop (Static Profile; "<<useKokkos<<")"    << std::endl
        << "================================================================================" << std::endl
        << std::endl;

    int status = 0;
    for(size_t i=0; i<opts.repetitions; ++i)
      {
        if(opts.useKokkosAssembly)
          status += executeTotalElementLoopSPKokkos_(comm, opts);
        else
          status += executeTotalElementLoopSP_(comm, opts);
      }
    return status;
  }



  int executeTotalElementLoopSP_(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
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

    int maxEntriesPerRow = 16;

    // Build Tpetra Maps
    // -----------------
    // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
    RCP<const map_type> row_map = rcp(new map_type(GO_INVALID, mesh.getOwnedNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
    RCP<const map_type> owned_element_map = rcp(new map_type(GO_INVALID, mesh.getOwnedElementGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
    RCP<const map_type> ghost_element_map = rcp(new map_type(GO_INVALID, mesh.getGhostElementGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
    RCP<const import_type> elementImporter = rcp(new import_type(owned_element_map,ghost_element_map));

    if(opts.verbose) row_map->describe(out);

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
    for(size_t element_gidx=0; element_gidx<mesh.getNumOwnedElements(); element_gidx++)
      {
        // Populate global_ids_in_row:
        // - Copy the global node ids for current owned element into an array.
        // - Since each element's contribution is a clique, we can re-use this for
        //   each row associated with this element's contribution.
        for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
          {
            global_ids_in_row[element_node_idx] = owned_element_to_node_ids(element_gidx, element_node_idx);
          }

        // Add the contributions from the current row into the graph if the node is owned.
        // - For example, if Element 0 contains nodes [0,1,4,5] and nodes 0 and 4 are owned
        //   by the current processor, then:
        //   - node 0 inserts [0, 1, 4, 5]
        //   - node 1 <skip>
        //   - node 4 inserts [0, 1, 4, 5]
        //   - node 5 <skip>
        for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
          {
            if(mesh.nodeIsOwned(global_ids_in_row[element_node_idx]))
              {
                crs_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
              }
          }
      }

    // Insert the node contributions for every GHOST element:
    for(size_t element_gidx=0; element_gidx<mesh.getNumGhostElements(); element_gidx++)
      {
        for(size_t element_node_idx=0; element_node_idx<ghost_element_to_node_ids.extent(1); element_node_idx++)
          {
            global_ids_in_row[element_node_idx] = ghost_element_to_node_ids(element_gidx, element_node_idx);
          }
        for(size_t element_node_idx=0; element_node_idx<ghost_element_to_node_ids.extent(1); element_node_idx++)
          {
            if(mesh.nodeIsOwned(global_ids_in_row[element_node_idx]))
              {
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
    if(opts.verbose) crs_graph->describe(out, Teuchos::VERB_EXTREME);


    // Simulated Ghosting of Material State
    // -------------------
    {
      GhostState state(elementImporter,opts.numStateDoublesPerElement);
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

    crs_matrix_type crs_matrix (crs_graph);
    multivector_type rhs (crs_graph->getRowMap(), 1);

    Kokkos::View<Scalar[4][4], hostType> element_matrix ("element_matrix");
    Teuchos::Array<Scalar> element_rhs(4);

    Teuchos::Array<global_ordinal_type> column_global_ids(4);     // global column ids list
    Teuchos::Array<Scalar> column_scalar_values(4);         // scalar values for each column


    // Loop over owned elements:
    for (size_t element_gidx = 0;
         element_gidx < mesh.getNumOwnedElements (); ++element_gidx) {
      // Get the stiffness matrix for this element
      ReferenceQuad4(element_matrix);
      ReferenceQuad4RHS(element_rhs);

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
        if (mesh.nodeIsOwned (global_row_id)) {
          for (size_t col_idx = 0; col_idx < 4; ++col_idx) {
            column_scalar_values[col_idx] =
              element_matrix(element_node_idx, col_idx);
          }
          crs_matrix.sumIntoGlobalValues (global_row_id,
                                          column_global_ids,
                                          column_scalar_values);
          rhs.sumIntoGlobalValue (global_row_id, 0,
                                  element_rhs[element_node_idx]);
        }
      }
    }

    // Loop over ghost elements:
    // - This loop is the same as the element loop for owned elements, but this one
    //   is for ghost elements.
    for(size_t element_gidx=0; element_gidx<mesh.getNumGhostElements(); element_gidx++)
      {
        ReferenceQuad4(element_matrix);
        ReferenceQuad4RHS(element_rhs);

        for(size_t element_node_idx=0; element_node_idx<ghost_element_to_node_ids.extent(1); element_node_idx++)
          {
            column_global_ids[element_node_idx] = ghost_element_to_node_ids(element_gidx, element_node_idx);
          }

        for(size_t element_node_idx=0; element_node_idx<4; element_node_idx++)
          {
            global_ordinal_type global_row_id = ghost_element_to_node_ids(element_gidx, element_node_idx);
            if(mesh.nodeIsOwned(global_row_id))
              {
                for(size_t col_idx=0; col_idx<4; col_idx++)
                  {
                    column_scalar_values[col_idx] = element_matrix(element_node_idx, col_idx);
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
      crs_matrix.describe (out, Teuchos::VERB_EXTREME);
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

  int
  executeTotalElementLoopSPKokkos_
  (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
   const struct CmdLineOpts& opts)
  {
    using Teuchos::RCP;
    using Teuchos::TimeMonitor;

    const global_ordinal_type GO_INVALID = Teuchos::OrdinalTraits<global_ordinal_type>::invalid();
    const local_ordinal_type LO_INVALID = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();
    using pair_type = Kokkos::pair<int,int>;

    // The output stream 'out' will ignore any output not from Process 0.
    RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
    Teuchos::FancyOStream& out = *pOut;

    // Processor decomp (only works on perfect squares)
    int numProcs  = comm->getSize();
    int sqrtProcs = sqrt(numProcs);

    if(sqrtProcs*sqrtProcs != numProcs) {
      if(0 == comm->getRank()) {
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

    MeshDatabase mesh(comm,nex,ney,procx,procy);

    if(opts.verbose) mesh.print(std::cout);

    int maxEntriesPerRow = 16;

    // Build Tpetra Maps
    // -----------------
    // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
    RCP<const map_type> row_map =
      rcp(new map_type(GO_INVALID, mesh.getOwnedNodeGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
    RCP<const map_type> owned_element_map = rcp(new map_type(GO_INVALID, mesh.getOwnedElementGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
    RCP<const map_type> ghost_element_map = rcp(new map_type(GO_INVALID, mesh.getGhostElementGlobalIDs().getDeviceView(Tpetra::Access::ReadOnly), 0, comm));
    RCP<const import_type> elementImporter = rcp(new import_type(owned_element_map,ghost_element_map));

    if(opts.verbose) row_map->describe(out);

    // Graph Construction
    // ------------------
    // - Loop over every element in the mesh.
    //   - Get list of nodes associated with each element.
    //   - Insert node contributions if the node (row) is owned locally.
    //
    auto domain_map = row_map;
    auto range_map  = row_map;


    Teuchos::TimeMonitor::getStackedTimer()->startBaseTimer(); 
    RCP<TimeMonitor> timerElementLoopGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("1) ElementLoop  (Graph)")));
    RCP<crs_graph_type> crs_graph = rcp(new crs_graph_type(row_map, maxEntriesPerRow));
    {
      auto owned_element_to_node_ids = mesh.getOwnedElementToNode().getHostView(Tpetra::Access::ReadOnly);
      auto ghost_element_to_node_ids = mesh.getGhostElementToNode().getHostView(Tpetra::Access::ReadOnly);
   
      // Using 4 because we're using quads for this example, so there will be 4 nodes associated with each element.
      Teuchos::Array<global_ordinal_type> global_ids_in_row(4);
  
      // Insert node contributions for every OWNED element:
      for(size_t element_gidx=0; element_gidx<mesh.getNumOwnedElements(); element_gidx++)
        {
          // Populate global_ids_in_row:
          // - Copy the global node ids for current owned element into an array.
          // - Since each element's contribution is a clique, we can re-use this for
          //   each row associated with this element's contribution.
          for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
            {
              global_ids_in_row[element_node_idx] = owned_element_to_node_ids(element_gidx, element_node_idx);
            }
      
          // Add the contributions from the current row into the graph if the node is owned.
          // - For example, if Element 0 contains nodes [0,1,4,5] and nodes 0 and 4 are owned
          //   by the current processor, then:
          //   - node 0 inserts [0, 1, 4, 5]
          //   - node 1 <skip>
          //   - node 4 inserts [0, 1, 4, 5]
          //   - node 5 <skip>
          for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
            {
              if(mesh.nodeIsOwned(global_ids_in_row[element_node_idx]))
                {
                  crs_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
                }
            }
        }
  
      // Insert the node contributions for every GHOST element:
      for(size_t element_gidx=0; element_gidx<mesh.getNumGhostElements(); element_gidx++)
        {
          for(size_t element_node_idx=0; element_node_idx<ghost_element_to_node_ids.extent(1); element_node_idx++)
            {
              global_ids_in_row[element_node_idx] = ghost_element_to_node_ids(element_gidx, element_node_idx);
            }
          for(size_t element_node_idx=0; element_node_idx<ghost_element_to_node_ids.extent(1); element_node_idx++)
            {
              if(mesh.nodeIsOwned(global_ids_in_row[element_node_idx]))
                {
                  crs_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
                }
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
    if(opts.verbose) crs_graph->describe(out, Teuchos::VERB_EXTREME);

    // Simulated Ghosting of Material State
    // -------------------
    {
      GhostState state(elementImporter,opts.numStateDoublesPerElement);
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
    RCP<multivector_type> rhs = rcp(new multivector_type(crs_graph->getRowMap(), 1));
    {

      auto owned_element_to_node_ids = mesh.getOwnedElementToNode().getDeviceView(Tpetra::Access::ReadOnly);
      auto ghost_element_to_node_ids = mesh.getGhostElementToNode().getDeviceView(Tpetra::Access::ReadOnly);    
      auto localMatrix  = crs_matrix->getLocalMatrixDevice();
      auto localRHS     = rhs->getLocalViewDevice(Tpetra::Access::OverwriteAll);
      auto localRowMap  = crs_matrix->getRowMap()->getLocalMap();
      auto localColMap  = crs_matrix->getColMap()->getLocalMap();
    
      // Because we're processing elements in parallel, we need storage for all of them
      int numOwnedElements = mesh.getNumOwnedElements();
      int numGhostElements = mesh.getNumGhostElements();
      int nperel = owned_element_to_node_ids.extent(1);
      pair_type alln = pair_type(0,nperel);
      scalar_2d_array_type all_element_matrix("all_element_matrix",nperel*std::max(numOwnedElements,numGhostElements));
      scalar_1d_array_type all_element_rhs("all_element_rhs",nperel*std::max(numOwnedElements,numGhostElements));
      local_ordinal_single_view_type  all_lcids("all_lids",nperel*std::max(numOwnedElements,numGhostElements));
    
    
      timerElementLoopMemory = Teuchos::null;
      RCP<TimeMonitor> timerElementLoopMatrix = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("3.3) ElementLoop  (Matrix)")));

      // Work around subview of managed views being slower than unmanaged
      auto all_element_rhs_unmanaged = makeUnmanaged(all_element_rhs);
      auto all_element_matrix_unmanaged = makeUnmanaged(all_element_matrix);
      auto all_lcids_unmanaged = makeUnmanaged(all_lcids);
    
      // Loop over owned elements:
      Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, numOwnedElements),KOKKOS_LAMBDA(const size_t& element_gidx) {
          // Get subviews
          pair_type location_pair = pair_type(nperel*element_gidx,nperel*(element_gidx+1));
          auto element_rhs    = Kokkos::subview(all_element_rhs_unmanaged,location_pair);
          auto element_matrix = Kokkos::subview(all_element_matrix_unmanaged,location_pair,alln);
          auto element_lcids  = Kokkos::subview(all_lcids_unmanaged,location_pair);
      
          // Get the contributions for the current element
          ReferenceQuad4(element_matrix);
          ReferenceQuad4RHS(element_rhs);
      
          // Get the local column ids array for this element
          for(int element_node_idx=0; element_node_idx<nperel; element_node_idx++) {
            element_lcids(element_node_idx) = localColMap.getLocalElement(owned_element_to_node_ids(element_gidx, element_node_idx));
          }

          // For each node (row) on the current element:
          // - populate the values array
          // - add the values to the fe_matrix.
          for(int element_node_idx=0; element_node_idx<nperel; element_node_idx++)
            {
              global_ordinal_type global_row_id = owned_element_to_node_ids(element_gidx, element_node_idx);
              local_ordinal_type local_row_id = localRowMap.getLocalElement(global_row_id);
              if(local_row_id != LO_INVALID) {
                // Force atomics on sums
                for(int col_idx=0; col_idx<nperel; col_idx++)
                  localMatrix.sumIntoValues(local_row_id,&element_lcids(col_idx),1,&(element_matrix(element_node_idx,col_idx)),true,true);
                Kokkos::atomic_add(&(localRHS(local_row_id,0)),element_rhs[element_node_idx]);
              }
            }
        });
      execution_space ().fence ();

      // Loop over ghost elements:
      // - This loop is the same as the element loop for owned elements, but this one
      //   is for ghost elements.
      Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, numGhostElements),KOKKOS_LAMBDA(const size_t& element_gidx) {
          // Get subviews
          pair_type location_pair = pair_type(nperel*element_gidx,nperel*(element_gidx+1));
          auto element_rhs    = Kokkos::subview(all_element_rhs_unmanaged,location_pair);
          auto element_matrix = Kokkos::subview(all_element_matrix_unmanaged,location_pair,alln);
          auto element_lcids  = Kokkos::subview(all_lcids_unmanaged,location_pair);

          // Get the contributions for the current element
          ReferenceQuad4(element_matrix);
          ReferenceQuad4RHS(element_rhs);

          // Get the local column ids array for this element
          for(int element_node_idx=0; element_node_idx<nperel; element_node_idx++) {
            element_lcids(element_node_idx) = localColMap.getLocalElement(ghost_element_to_node_ids(element_gidx, element_node_idx));
          }

          for(int element_node_idx=0; element_node_idx<nperel; element_node_idx++)
            {
              global_ordinal_type global_row_id = ghost_element_to_node_ids(element_gidx, element_node_idx);
              local_ordinal_type local_row_id = localRowMap.getLocalElement(global_row_id);
              if(local_row_id != LO_INVALID) {
                // Force atomics on sums
                for(int col_idx=0; col_idx<nperel; col_idx++)
                  localMatrix.sumIntoValues(local_row_id,&element_lcids(col_idx),1,&(element_matrix(element_node_idx,col_idx)),true,true);
                Kokkos::atomic_add(&(localRHS(local_row_id,0)),element_rhs[element_node_idx]);
              }
            }
        });
      execution_space ().fence ();    
    timerElementLoopMatrix = Teuchos::null;
    }

    // After the contributions are added, 'finalize' the matrix using fillComplete()
    {
      TimeMonitor timer(*TimeMonitor::getNewTimer("4) FillComplete (Matrix)"));
      crs_matrix->fillComplete();
    }

    // Print out crs_matrix details.
    if(opts.verbose) crs_matrix->describe(out, Teuchos::VERB_EXTREME);

    Teuchos::TimeMonitor::getStackedTimer()->stopBaseTimer();

    // Save crs_matrix as a MatrixMarket file.
    if(opts.saveMM)
      {
        std::ofstream ofs("crsMatrix_TotalElementLoop_SPKokkos.out", std::ofstream::out);
        Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparse(ofs, crs_matrix);
        std::ofstream ofs2("rhs_TotalElementLoop_SPKokkos.out", std::ofstream::out);
        Tpetra::MatrixMarket::Writer<multivector_type>::writeDense(ofs2, rhs);
      }

    return 0;
  }


} // namespace TpetraExamples

#endif // TPETRAEXAMPLES_FEM_ASSEMBLY_TOTALELEMENTLOOP_HPP
