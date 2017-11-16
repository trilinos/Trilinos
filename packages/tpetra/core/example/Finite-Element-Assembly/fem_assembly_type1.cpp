#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "typedefs.hpp"
#include "MeshDatabase.hpp"
#include "Element.hpp"



int main (int argc, char *argv[]) 
{
  using Teuchos::RCP;

  const GlobalOrdinal GO_INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
  
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
  // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
  RCP<const MapType> row_map = rcp(new MapType(GO_INVALID, mesh.getOwnedNodeGlobalIDs(), 0, comm));

  auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  row_map->describe(*out);

  // Type-1: Graph Construction
  // --------------------------
  // - Loop over every element in the mesh.
  //   - Get list of nodes associated with each element.
  //   - Insert the clique of nodes associated with each element into the graph. 
  //   
  auto domain_map = row_map;
  auto range_map  = row_map;

  RCP<GraphType> crs_graph = rcp(new GraphType(row_map, 0));
  auto owned_element_to_node_ids = mesh.getOwnedElementToNode(); 
  
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
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
       crs_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
    }
  }

  // Call fillComplete on the crs_graph to 'finalize' it.
  crs_graph->fillComplete();

  // Let's see what we have
  crs_graph->describe(*out, Teuchos::VERB_EXTREME);


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

  RCP<MatrixType> crs_matrix = rcp(new MatrixType(crs_graph));

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
    // - add the values to the crs_matrix.
    // Note: hardcoded 4 here because we're using quads.
    for(size_t element_node_idx=0; element_node_idx<4; element_node_idx++)
    { 
      GlobalOrdinal global_row_id = owned_element_to_node_ids(element_gidx, element_node_idx);

      for(size_t col_idx=0; col_idx<4; col_idx++)
      {
        column_scalar_values[col_idx] = element_matrix(element_node_idx, col_idx);
      }

      crs_matrix->sumIntoGlobalValues(global_row_id, column_global_ids, column_scalar_values);

#if 1
      // Print out the row's contribution...
      if(1 == mpiSession.getNProc()) 
      {
        std::cout << "- insertGlobalValues(" << std::setw(2) << global_row_id;
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

  // After the contributions are added, 'finalize' the matrix using fillComplete()
  crs_matrix->fillComplete();
  crs_matrix->describe(*out, Teuchos::VERB_EXTREME);
  
  // Finalize Kokkos
  Kokkos::finalize();

  return 0;
}  // END main()
