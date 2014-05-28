/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <iostream>
#include <Shards_BasicTopologies.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

static const size_t NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

typedef stk::mesh::Field<double>                          ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>    VectorField ;

void fill_utest_mesh_meta_data(stk::mesh::fem::FEMMetaData& fem_meta)
{
  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

  stk::mesh::Part& elem_block = fem_meta.declare_part( "block_1" , element_rank );

  stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::fem::set_cell_topology( elem_block, hex_top );

  const unsigned number_of_states = 1;

  ScalarField& temperature_field = fem_meta.declare_field<ScalarField>( "temperature", number_of_states );
  ScalarField& pressure_field = fem_meta.declare_field<ScalarField>( "pressure", number_of_states );
  VectorField& velocity_field = fem_meta.declare_field<VectorField>( "velocity",    number_of_states );

  stk::mesh::put_field( temperature_field, NODE_RANK,     elem_block );
  stk::mesh::put_field( pressure_field,    element_rank, elem_block );
  stk::mesh::put_field( velocity_field,    NODE_RANK,     elem_block );

  fem_meta.commit();
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

  stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(bulk_data);
  stk::mesh::Part& elem_block = *(fem_meta.get_part( "block_1" ));

  for(unsigned i = 0; i<num_elems; ++i) {
    elem_node_ids( elem_id+i, node_ids );
    stk::mesh::fem::declare_element( bulk_data, elem_block, elem_id+i, node_ids );
  }
  bulk_data.modification_end();
}

