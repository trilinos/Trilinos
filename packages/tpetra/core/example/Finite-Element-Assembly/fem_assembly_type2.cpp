// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "typedefs.hpp"
#include "MeshDatabase.hpp"
#include "Element.hpp"

#define PRINT_VERBOSE 0



int main (int argc, char *argv[]) 
{
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  const GlobalOrdinal GO_INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

  auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  
  // MPI boilerplate
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm();

  // Initialize Kokkos
  Kokkos::initialize();

  // Processor decomp (only works on perfect squares)
  int numProcs  = comm->getSize();
  int sqrtProcs = sqrt(numProcs); 

  if(sqrtProcs*sqrtProcs != numProcs) 
  {
    if(0 == mpiSession.getRank())
      std::cerr << "Error: Invalid number of processors provided, num processors must be a perfect square." << std::endl;
    return -1;
  }
  int procx = sqrtProcs;
  int procy = sqrtProcs;

  // Generate a simple 3x3 mesh
  int nex = 3;
  int ney = 3;
  MeshDatabase mesh(comm,nex,ney,procx,procy);
  mesh.print(std::cout);

  // Build Tpetra Maps
  // -----------------
  // - Doxygen: https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
  RCP<const MapType> owned_row_map       = rcp(new MapType(GO_INVALID, mesh.getOwnedNodeGlobalIDs(), 0, comm));
  RCP<const MapType> overlapping_row_map = rcp(new MapType(GO_INVALID, mesh.getOwnedAndGhostNodeGlobalIDs(), 0, comm));
  ExportType exporter(overlapping_row_map, owned_row_map); 

  #if PRINT_VERBOSE
  owned_row_map->describe(*out);
  overlapping_row_map->describe(*out);
  #endif

  // Type-2: Graph Construction
  // --------------------------
  auto domain_map = owned_row_map;
  auto range_map  = owned_row_map;

  auto owned_element_to_node_ids = mesh.getOwnedElementToNode();

  RCP<TimeMonitor> timerGlobal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("X) Global")));
  RCP<TimeMonitor> timerElementLoopGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("1) ElementLoop  (All Graph)")));

  // Type-2 Assembly distinguishes owned and overlapping nodes.
  // - Owned nodes are nodes that only touch elements owned by the same process.
  // - Overlapping nodes are nodes which touch elements owned by a different process.
  //
  // In Type-2 assembly, the graph construction loop looks similar to Type-1 but
  // in this case we insert rows into crs_graph_overlapping, then we fillComplete
  // the overlapping graph.  Next we export contributions from overlapping graph 
  // to the owned graph and call fillComplete on the owned graph.
  //
  RCP<GraphType> crs_graph_owned = rcp(new GraphType(owned_row_map, 0));
  RCP<GraphType> crs_graph_overlapping = rcp(new GraphType(overlapping_row_map, 0));

  // Note: Using 4 because we're using quads for this example, so there will be 4 nodes associated with each element.
  Teuchos::Array<GlobalOrdinal> global_ids_in_row(4);

  // for each element in the mesh...
  for(size_t element_gidx=0; element_gidx<mesh.getNumOwnedElements(); element_gidx++)
  {
    // Populate global_ids_in_row:
    // - Copy the global node ids for current element into an array.
    // - Since each element's contribution is a clique, we can re-use this for 
    //   each row associated with this element's contribution.
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
      global_ids_in_row[element_node_idx] = owned_element_to_node_ids(element_gidx, element_node_idx);
    }

    // Add the contributions from the current row into the overlapping graph.
    // - For example, if Element 0 contains nodes [0,1,4,5] then we insert the nodes:
    //   - node 0 inserts [0, 1, 4, 5]
    //   - node 1 inserts [0, 1, 4, 5]
    //   - node 4 inserts [0, 1, 4, 5]
    //   - node 5 inserts [0, 1, 4, 5]
    //
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
       crs_graph_overlapping->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
    }
  }
  timerElementLoopGraph = Teuchos::null;

  // Call fillComplete on the crs_graph_owned to 'finalize' it.
  {
    RCP<TimeMonitor> timerFillCompleteOverlappingGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("2) FillComplete (Overlapping Graph)")));
    crs_graph_overlapping->fillComplete();
  }

  // Need to Export and fillComplete the crs_graph_owned structure...
  // NOTE: Need to implement a graph transferAndFillComplete() method.
  {
    RCP<TimeMonitor> timerExportOwnedGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("3) Export       (Owned Graph)")));
    crs_graph_owned->doExport(*crs_graph_overlapping, exporter, Tpetra::INSERT);
  }

  {
    RCP<TimeMonitor> timerFillCompleteOwnedGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("4) FillComplete (Owned Graph)")));
    crs_graph_owned->fillComplete();
  }

  // Let's see what we have
  #if PRINT_VERBOSE
  crs_graph_owned->describe(*out, Teuchos::VERB_EXTREME);
  crs_graph_overlapping->describe(*out, Teuchos::VERB_EXTREME);
  #endif

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
  // For Type 2 matrix fill, we create a crs_matrix object for both owned
  // and overlapping rows.  We will only fill the overlapping graph using
  // the same method as we filled the graph but in this case, nodes 
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
  RCP<TimeMonitor> timerElementLoopMatrix = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("5) ElementLoop  (All Matrix)")));

  // Create owned and overlapping CRS Matrices
  RCP<MatrixType> crs_matrix_owned       = rcp(new MatrixType(crs_graph_owned));
  RCP<MatrixType> crs_matrix_overlapping = rcp(new MatrixType(crs_graph_overlapping));

  scalar_2d_array_type element_matrix;
  Kokkos::resize(element_matrix, 4, 4);

  Teuchos::Array<GlobalOrdinal> column_global_ids(4);     // global column ids list
  Teuchos::Array<Scalar> column_scalar_values(4);         // scalar values for each column

  // Loop over elements
  for(size_t element_gidx=0; element_gidx<mesh.getNumOwnedElements(); element_gidx++)
  {
    // Get the contributions for the current element
    ReferenceQuad4(element_matrix);

    // Fill the global column ids array for this element
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
      column_global_ids[element_node_idx] = owned_element_to_node_ids(element_gidx, element_node_idx);
    }
    
    // For each node (row) on the current element:
    // - populate the values array
    // - add the values to the crs_matrix_owned.
    // Note: hardcoded 4 here because we're using quads.
    for(size_t element_node_idx=0; element_node_idx<4; element_node_idx++)
    { 
      GlobalOrdinal global_row_id = owned_element_to_node_ids(element_gidx, element_node_idx);

      for(size_t col_idx=0; col_idx<4; col_idx++)
      {
        column_scalar_values[col_idx] = element_matrix(element_node_idx, col_idx);
      }

      // For Type-2 Assembly, we only sumInot the overlapping crs_matrix.
      crs_matrix_overlapping->sumIntoGlobalValues(global_row_id, column_global_ids, column_scalar_values);
    }
  }
  timerElementLoopMatrix = Teuchos::null;

  // After contributions are added, we finalize the crs_matrix in 
  // the same manner that we did the crs_graph.  
  // On Type-2 assembly, we fillComplete the overlapping matrix, then
  // export contributions to the owned matrix using the exporter, then
  // fillComplete the owned matrix.
  {
    RCP<TimeMonitor> timerFillCompleteOverlappingMatrix = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("6) FillComplete (Overlapping Matrix)")));
    crs_matrix_overlapping->fillComplete();
  }

  {
    RCP<TimeMonitor> timerExportOwnedMatrix = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("7) Export       (Owned Matrix)")));
    crs_matrix_owned->doExport(*crs_matrix_overlapping, exporter, Tpetra::ADD);
  }
  
  {
    RCP<TimeMonitor> timerFillCompleteOwnedMatrix = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("8) FillComplete (Owned Matrix)")));
    crs_matrix_owned->fillComplete();
  }

  timerGlobal = Teuchos::null;

  // Print out crs_matrix_owned and crs_matrix_overlapping details.
  #if PRINT_VERBOSE
  crs_matrix_owned->describe(*out, Teuchos::VERB_EXTREME);
  crs_matrix_overlapping->describe(*out, Teuchos::VERB_EXTREME);
  #endif

  // Save crs_matrix as a MatrixMarket file.
  std::ofstream ofs("Finite-Element-Matrix-Assembly_Type2.out", std::ofstream::out);
  Tpetra::MatrixMarket::Writer<MatrixType>::writeSparse(ofs, crs_matrix_owned);
  ofs.close();

  // Print out timing results.
  TimeMonitor::report(comm.ptr(), std::cout, "");

  // Finalize Kokkos
  Kokkos::finalize();
 
  // This tells the Trilinos test framework that the test passed.
  if(0 == mpiSession.getRank())
  {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}  // END main()

