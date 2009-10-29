
#include <iostream>
#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

void fill_utest_mesh_meta_data(stk::mesh::MetaData& meta_data, bool use_temperature)
{
  stk::mesh::Part& elem_block = meta_data.declare_part( "block_1" , stk::mesh::Element );

  stk::mesh::set_cell_topology< shards::Hexahedron<8> >( elem_block );

  const unsigned number_of_states = 1;

  if (use_temperature) {
    stk::mesh::ScalarField& temperature_field =
      meta_data.declare_field<stk::mesh::ScalarField>( "temperature", number_of_states );
    stk::mesh::put_field( temperature_field, stk::mesh::Node,    elem_block );
  }

  stk::mesh::ScalarField& pressure_field =
    meta_data.declare_field<stk::mesh::ScalarField>( "pressure", number_of_states );
  stk::mesh::VectorField& velocity_field =
    meta_data.declare_field<stk::mesh::VectorField>( "velocity",    number_of_states );

  stk::mesh::put_field( pressure_field,    stk::mesh::Element, elem_block );
  stk::mesh::put_field( velocity_field,    stk::mesh::Node,    elem_block );
}

// Assume the following mesh of 4 hex8 elements on the first processor (proc 0).
//
// Global node and element numbering
//     3       7      11      15      19
//     +-------+-------+-------+-------+
//    /       /       /       /       /|
//  4/      8/     12/     16/     20/ |
//  +-------+-------+-------+-------+  |
//  |       |       |       |       |  +18
//  |  e1   |  e2   |  e3   |  e4   | /
//  |       |       |       |       |/
//  +-------+-------+-------+-------+
//  1       5      9       13      17
//
// Local node numbering
//     8       7
//     +-------+
//    /       /|
//  5/      6/ |
//  +-------+  |
//  |       |  +3
//  |  e1   | /
//  |       |/
//  +-------+
//  1       2
//
// A similar mesh will be created on the other processors, but using different elem-ids.
// i.e., proc 0 has elements 1 - 4, proc 1 has elements 5 - 8, and so on. 4 nodes will be
// shared between each neighboring pair of processors.
//----------------------------------------------------------------------

void elem_node_ids( stk::mesh::EntityId elem_id , stk::mesh::EntityId node_ids[] )
{
  if ( elem_id == 0 ) {
    std::cout << "use_case_1, elem_node_ids: ERROR, elem_id ("
        << elem_id << ") must be greater than 0." << std::endl;
    return;
  }

  const unsigned base = ( elem_id - 1 ) * 4 ;
  node_ids[0] = base + 1 ;
  node_ids[1] = base + 5 ;
  node_ids[2] = base + 6 ;
  node_ids[3] = base + 2 ;
  node_ids[4] = base + 4 ;
  node_ids[5] = base + 8 ;
  node_ids[6] = base + 7 ;
  node_ids[7] = base + 3 ;
}

void fill_utest_mesh_bulk_data(stk::mesh::BulkData& bulk_data)
{
  const unsigned num_elems = 4;

  int numProcs = 1, myProc = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);

  stk::mesh::EntityId elem_id = 1 + myProc*num_elems;
  stk::mesh::EntityId node_ids[ shards::Hexahedron<8>::node_count ];

  stk::mesh::Part& elem_block = *(bulk_data.mesh_meta_data().get_part( "block_1" ));

  for(unsigned i = 0; i<num_elems; ++i) {
    elem_node_ids( elem_id+i, node_ids );
    stk::mesh::declare_element( bulk_data, elem_block, elem_id+i, node_ids );
  }
}

