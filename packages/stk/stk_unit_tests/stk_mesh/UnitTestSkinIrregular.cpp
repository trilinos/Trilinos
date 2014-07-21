/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>                    // for max
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_mesh/base/Types.hpp>      // for EntityId, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Bucket; } }








using stk::mesh::EntityId;
using stk::mesh::EntityRank;
using stk::mesh::MetaData;

//---------------------------------------------------------------------------------------

TEST( UnitTestSkin, SkinPocket)
{
  enum { SpatialDim = 3 };

  stk::ParallelMachine pm = MPI_COMM_WORLD ;
  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::MetaData fem_meta;
  fem_meta.initialize(SpatialDim);
  stk::mesh::BulkData bulk_data( fem_meta , pm );
  stk::mesh::Part & hex_part = fem_meta.declare_part_with_topology( "hex_part", stk::topology::HEX_8 );
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;
  const EntityRank side_rank    = fem_meta.side_rank();

  //create and skin a 2 hex-element mesh with a pocket
  //in a normal mesh 6 and 13 would be the same node
  //
  //    8-------7-------12
  //   /|      /|\     /|
  //  / |     / | \   / |
  // 5-------6  | 13-11 |
  // |  4----|--3/---|--10
  // | /     | //    | /
  // |/      |//     |/
  // 1-------2-------9
  //

  fem_meta.commit();

  bulk_data.modification_begin();

  // declare left element on first process
  if (p_rank == 0)
  {
    EntityId element_id = 1;
    EntityId node_ids[8] = { 1, 2, 3, 4, 5, 6, 7, 8};

    stk::mesh::declare_element( bulk_data, hex_part, element_id, node_ids);

  }

  // declare right element on last process
  if (p_rank == p_size -1)
  {
    EntityId element_id = 2;
    EntityId node_ids[8] = { 2, 9, 10, 3, 13, 11, 12, 7};

    stk::mesh::declare_element( bulk_data, hex_part, element_id, node_ids);
  }

  bulk_data.modification_end();

  //skin the mesh
  stk::mesh::skin_mesh(bulk_data);

  //each element should have 6 faces attached to it
  for (EntityId element_id = 1; element_id < 3; ++element_id) {
    stk::mesh::Entity element = bulk_data.get_entity( element_rank, element_id);
    if ( bulk_data.is_valid(element) ) {
      EXPECT_EQ( bulk_data.num_connectivity(element, side_rank), 6u);
    }
  }
}

