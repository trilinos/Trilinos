/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <iostream>
#include <Shards_BasicTopologies.hpp>

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
#include <stk_mesh/fem/DefaultFEM.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

using stk::mesh::fem::NODE_RANK;

typedef stk::mesh::Field<double>                          ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>    VectorField ;

void fill_utest_mesh_meta_data(stk::mesh::MetaData& meta_data)
{
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const stk::mesh::EntityRank element_rank = stk::mesh::Element;
#else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
  stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(meta_data);
  const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(fem);
#endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */

  stk::mesh::Part& elem_block = meta_data.declare_part( "block_1" , element_rank );


#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  stk::mesh::set_cell_topology< shards::Hexahedron<8> >( elem_block );
#else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
  stk::mesh::fem::set_cell_topology< shards::Hexahedron<8> >( elem_block );
#endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */

  const unsigned number_of_states = 1;

  ScalarField& temperature_field = meta_data.declare_field<ScalarField>( "temperature", number_of_states );
  ScalarField& pressure_field = meta_data.declare_field<ScalarField>( "pressure", number_of_states );
  VectorField& velocity_field = meta_data.declare_field<VectorField>( "velocity",    number_of_states );

  stk::mesh::put_field( temperature_field, NODE_RANK,     elem_block );
  stk::mesh::put_field( pressure_field,    element_rank, elem_block );
  stk::mesh::put_field( velocity_field,    NODE_RANK,     elem_block );

  meta_data.commit();
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
  bulk_data.modification_begin();
  const unsigned num_elems = 4;

  const unsigned myProc = bulk_data.parallel_rank();

  stk::mesh::EntityId elem_id = 1 + myProc*num_elems;
  stk::mesh::EntityId node_ids[ shards::Hexahedron<8>::node_count ];

  stk::mesh::Part& elem_block = *(bulk_data.mesh_meta_data().get_part( "block_1" ));

  for(unsigned i = 0; i<num_elems; ++i) {
    elem_node_ids( elem_id+i, node_ids );
    stk::mesh::declare_element( bulk_data, elem_block, elem_id+i, node_ids );
  }
  bulk_data.modification_end();
}

