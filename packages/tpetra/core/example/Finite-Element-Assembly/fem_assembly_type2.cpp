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


#define DO_GRAPH_ASSEMBLY  1
#define DO_MATRIX_ASSEMBLY 1


int main (int argc, char *argv[]) 
{
  using Teuchos::RCP;

  const GlobalOrdinal GO_INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
  
  // MPI boilerplate
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm();

  // Initialize Kokkos
  Kokkos::initialize();

  auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));

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

  // RCP<const MapType> row_map_overlapping = rcp(new MapType(GO_INVALID, mesh.getGhostNodeGlobalIDs(), 0, comm));        // wcmclen

  // TODO: [Type-2 Assembly]
  //       - Rename row_map to owned_row_map
  //       - owned_row_map is exactly the same as type-1 row map.
  //       - overlap map
  //         - contains all the OwnedNodeGlobalIDs followed by GhostNodeGlobalIDs
  //           in a single array containing both concatenated and as 2nd arg to map c'tor.

  owned_row_map->describe(*out);
  overlapping_row_map->describe(*out);

#if DO_GRAPH_ASSEMBLY 
  // Type-2: Graph Construction
  // --------------------------
  // - Loop over every element in the mesh.
  //   - Get list of nodes associated with each element.
  //   - Insert the clique of nodes associated with each element into the graph. 
  //   
  auto domain_map = owned_row_map;
  auto range_map  = owned_row_map;

  auto owned_element_to_node_ids = mesh.getOwnedElementToNode();

  RCP<GraphType> crs_graph_owned = rcp(new GraphType(owned_row_map, 0));
  RCP<GraphType> crs_graph_overlapping = rcp(new GraphType(overlapping_row_map, 0));

  // TODO: [Type-2 Assembly]
  //       - Create overlapping_graph and owned_graph.
  //         - overlapping_graph uses overlapping_map
  //         - owned_Graph uses owned_map

  // This is 4 because we're using quads for this example,
  // so there will be 4 nodes associated with each element.
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

    // Add the contributions from the current row into the graph.
    // - For example, if Element 0 contains nodes [0,1,4,5] then we insert the nodes:
    //   - node 0 inserts [0, 1, 4, 5]
    //   - node 1 inserts [0, 1, 4, 5]
    //   - node 4 inserts [0, 1, 4, 5]
    //   - node 5 inserts [0, 1, 4, 5]
    //
    // For Type-2 Graph Construction we only insert into crs_graph_overlapping
    //
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
       crs_graph_overlapping->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
    }
  }

  // Call fillComplete on the crs_graph_owned to 'finalize' it.
  crs_graph_overlapping->fillComplete();

  // Need to Export and fillComplete the crs_graph_owned structure...
  // NOTE: Need to implement a graph transferAndFillComplete() method.
  ExportType exporter(overlapping_row_map, owned_row_map); 

  //crs_graph->transferAndFillComplete();
  crs_graph_owned->doExport(*crs_graph_overlapping, exporter, Tpetra::INSERT);
  crs_graph_owned->fillComplete();

  // Let's see what we have
  crs_graph_owned->describe(*out, Teuchos::VERB_EXTREME);
  crs_graph_overlapping->describe(*out, Teuchos::VERB_EXTREME);

#endif   // DO_GRAPH_ASSEMBLY



#if DO_MATRIX_ASSEMBLY

  // TODO: [Type-2 Assembly]
  //       - Need to call fillComplete for only the overlapping_graph.
  //       - Siefert thinks there's a chance this could fail... so we'll see

  // Type-1: Matrix Fill
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
  // For Type 1 matrix fill, we create the crs_matrix object and will fill it 
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

  RCP<MatrixType> crs_matrix_owned       = rcp(new MatrixType(crs_graph_owned));
  RCP<MatrixType> crs_matrix_overlapping = rcp(new MatrixType(crs_graph_overlapping));

  // TODO: [Type-2 Assembly]
  //       - Need to create a crs_matrix_owned and crs_matrix_overlapping

  scalar_2d_array_type element_matrix;
  Kokkos::resize(element_matrix, 4, 4);

  Teuchos::Array<GlobalOrdinal> column_global_ids(4);     // global column ids list
  Teuchos::Array<Scalar> column_scalar_values(4);         // scalar values for each column

  // Loop over elements
  for(size_t element_gidx=0; element_gidx<mesh.getNumOwnedElements(); element_gidx++)
  {
#if 1
    if(1 == mpiSession.getNProc()) 
      std::cout << "Element " << element_gidx << std::endl;
#endif

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

      crs_matrix_overlapping->sumIntoGlobalValues(global_row_id, column_global_ids, column_scalar_values);

      // TODO: [Type-2 Assembly]
      //       - Only sumInto the overlapping crs_matrix_owned

#if 1
      // Print out the row's contribution...
      if(1 == mpiSession.getNProc()) 
      {
        std::cout << "- sumIntoGlobalValues(" << std::setw(2) << global_row_id;
        std::cout << ",  [  ";
        for(size_t i=0; i<4; i++)
        {
          std::cout << column_global_ids[i] << "  ";
        }
        std::cout << "]";
        std::cout << ",  [  ";
        for(size_t i=0; i<4; i++)
        {
          std::cout << column_scalar_values[i] << "  ";
        }
        std::cout << "])" << std::endl;
      }
#endif
    }
  }

  // TODO: [Type-2 Assembly]
  // - 

  // After the contributions are added, 'finalize' the matrix using fillComplete()
  crs_matrix_overlapping->fillComplete();

  crs_matrix_owned->doExport(*crs_matrix_overlapping, exporter, Tpetra::ADD);
  crs_matrix_owned->fillComplete();

  crs_matrix_owned->describe(*out, Teuchos::VERB_EXTREME);
  crs_matrix_overlapping->describe(*out, Teuchos::VERB_EXTREME);


  std::ofstream ofs("type2.out", std::ofstream::out);

  Tpetra::MatrixMarket::Writer<MatrixType>::writeSparse(ofs, crs_matrix_owned);

  ofs.close();

#endif  // DO_MATRIX_ASSEMBLY
  
  // Finalize Kokkos
  Kokkos::finalize();

  return 0;
}  // END main()


#if 0

  // Call fillComplete on the crs_graph_owned to 'finalize' it.
  crs_graph_overlapping->fillComplete();

  // Need to Export and fillComplete the crs_graph_owned structure...
  // NOTE: Need to implement a graph transferAndFillComplete() method.
  ExportType exporter(overlapping_row_map, owned_row_map); 

  //crs_graph->transferAndFillComplete();
  crs_graph_owned->doExport(*crs_graph_overlapping, exporter, Tpetra::INSERT);
  crs_graph_owned->fillComplete();

  // Let's see what we have
  crs_graph_owned->describe(*out, Teuchos::VERB_EXTREME);
  crs_graph_overlapping->describe(*out, Teuchos::VERB_EXTREME);



// From Lesson05
  RCP<crs_matrix_type> B;
  {
    // We created exportTimer in main().  It's a global timer.
    // Actually starting and stopping the timer is local, but
    // computing timer statistics (e.g., in TimeMonitor::summarize(),
    // called in main()) is global.  There are ways to restrict the
    // latter to any given MPI communicator; the default is
    // MPI_COMM_WORLD.
    TimeMonitor monitor (*exportTimer); // Time the redistribution

    // Make an export object with procZeroMap as the source Map, and
    // globalMap as the target Map.  The Export type has the same
    // template parameters as a Map.  Note that Export does not depend
    // on the Scalar template parameter of the objects it
    // redistributes.  You can reuse the same Export for different
    // Tpetra object types, or for Tpetra objects of the same type but
    // different Scalar template parameters (e.g., Scalar=float or
    // Scalar=double).
    typedef Tpetra::Export<> export_type;
    export_type exporter (procZeroMap, globalMap);

    // Make a new sparse matrix whose row map is the global Map.
    B = rcp (new crs_matrix_type (globalMap, 0));

    // Redistribute the data, NOT in place, from matrix A (which lives
    // entirely on Proc 0) to matrix B (which is distributed evenly over
    // the processes).
    B->doExport (*A, exporter, Tpetra::INSERT);
  }

  // We time redistribution of B separately from fillComplete().
  B->fillComplete ();
#endif


#if 0
Update
------
  1) Call constructors for both overlapped and non-overlapped graphs [same as we discussed]

  2) Loop over owned elements and insert rows into *only* the overlapped graph.

  3) Call fillComplete *only* on the overlapped graph.

  4) Build an exporter (source=overlapped map, target = nonoverlapped map)
     [there is a placeholder for this at the end of matrix; bump that up here.

  5) Export from the overlapped graph to the non-overlapped graph and fillComplete 
     (this can be done with the transferAndFillComplete method).

#endif