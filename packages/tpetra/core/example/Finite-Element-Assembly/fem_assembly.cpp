#include <cmath>
#include <iostream>
#include <iomanip>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "typedefs.hpp"
#include "MeshDatabase.hpp"


int main (int argc, char *argv[]) {
  using Teuchos::RCP;

  const GlobalOrdinal GO_INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
  
  // MPI boilerplate
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm();

  // Processor decomp (only works on perfect squares)
  int numProcs  = comm->getSize();
  int sqrtProcs = sqrt(numProcs); 
  if(sqrtProcs*sqrtProcs !=numProcs) return -1;
  int procx = sqrtProcs;
  int procy = sqrtProcs;

  // Type 1 & 2 modes supported so far by mesh database


  // Generate the mesh
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

  std::cout << "---------------------------------------" << std::endl;

  // Build graphs multiple ways
  // --------------------------

  // - Type 1 Graph Construction
  //   - Loop over owned elements only.
  //   row = domain = range
  auto domain_map = row_map;                                        // this right for type-1 assembly?
  auto range_map  = row_map;                                        // this right for type-1 assembly?
  RCP<GraphType> crs_graph = rcp(new GraphType(row_map, 0));        // TODO: maxNumEntriesPerRow should be set properly (9 in this case b/c 2D finite element)

  // Graph Construction
  // - Loop over every element in the mesh. Mesh->getNumOwnedElements()
  // - Get list of nodes associated with the element.   
  //     getOwnedElementToNode => a kokkos multidimensional array (don't look @ doxygen, the PDF in the trilinos repo is a better place to look) 
  //         - Note: local_ordinal_2d_array_type is typedef'ed from Kokkos::View<LocalOrdinal*[4],ExecutionSpace> 
  //     2D array::rows are elements & columns are nodes associated with the element
  //     (check ???)
  // - Insert the clique into the graph
  //    think crs_graph->insertGlobalIndices();
  auto owned_element_to_node_ids = mesh.getOwnedElementToNode();        // <-- owned_element_to_node_ids is a Kokkos::View
  //global_ordinal_view_type global_ids_in_row("global_ids_in_row", owned_element_to_node_ids.extent(1) );
  
  Teuchos::Array<GlobalOrdinal> global_ids_in_row(4);   // Teuchos Array of size 4 because this example uses quads.

  std::cout << "  # owned elements: " << owned_element_to_node_ids.extent(0) << std::endl;   // Number of owned elements  (extent() == dimension())
  std::cout << "  # nodes/element : " << owned_element_to_node_ids.extent(1) << std::endl;   // Number of nodes per element

  // for each element:
  for(size_t element_gidx=0; element_gidx<mesh.getNumOwnedElements(); element_gidx++)
  {
    std::cout << "element: " << element_gidx << std::endl;

    // Fill global_ids_in_row   
    // - Copy the node ids for the current element into an array containing the global
    //   ids for the row.  Since each element's contribution is a clique, we can use
    //   this for the contribution of every row associated with this element.
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
      global_ids_in_row[element_node_idx] = owned_element_to_node_ids(element_gidx, element_node_idx);
      std::cout << global_ids_in_row[element_node_idx] << " ";
    }
    std::cout << std::endl;

    // Add the contributions from the current row into the graph.
    // - The nodes in each element form a clique.
    // Loop over the nodes for each element and insert the adjacencies into the graph.
    // For example, if Element 0 contains nodes [0,1,4,5] then we insert the nodes:
    //   - node 0 inserts [0, 1, 4, 5]
    //   - node 1 inserts [0, 1, 4, 5]
    //   - node 4 inserts [0, 1, 4, 5]
    //   - node 5 inserts [0, 1, 4, 5]
    std::cout << "  ";
    // for each node_id in element i
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
       // std::cout << std::setw(2) << owned_element_to_node_ids(id, element_node_idx) << " ";
       crs_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
    }
    std::cout << std::endl;
    

    // loop over global_ids_in_row
    // for ...
    //    insert global_ids_in_row[] into matrix of row global_ids_in_row[j] using insert global ids.
    // 
    //    crs_graph->insertGlobalIndices(global_row, idx());
    //    signature:  insertGlobalIndices(const GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> & indices)
    //    NOTE:  Need a kokkos view version of insertGlobalIndices() // doesn't exist //
    //                
  }





//  crs_graph->fillComplete(domain_map, range_map);
  crs_graph->fillComplete();
  crs_graph->describe(*out, Teuchos::VERB_EXTREME);


  // - Type 2 Graph Construction
  //


  // Build matrices 
  // 


  // Build RHS vectors



  return 0;
}


