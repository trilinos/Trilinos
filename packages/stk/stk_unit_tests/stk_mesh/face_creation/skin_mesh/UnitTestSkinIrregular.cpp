// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <gtest/gtest.h>                             // for AssertHelper, etc
#include <stk_mesh/base/BulkData.hpp>                // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Entity.hpp>                  // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>              // for declare_element
#include <stk_mesh/base/MetaData.hpp>                // for MetaData
#include <stk_mesh/base/SkinMesh.hpp>                // for skin_mesh
#include <stk_mesh/base/Types.hpp>                   // for EntityId, etc
#include <stk_util/parallel/Parallel.hpp>
#include "mpi.h"                                     // for MPI_COMM_WORLD, etc
#include "stk_topology/topology.hpp"                 // for topology, etc
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
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

  stk::mesh::MeshBuilder builder(pm);
  builder.set_spatial_dimension(SpatialDim);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk_data = *bulkPtr;
  stk::mesh::MetaData& fem_meta = bulk_data.mesh_meta_data();
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
    stk::mesh::EntityIdVector node_ids { 1, 2, 3, 4, 5, 6, 7, 8};

    stk::mesh::declare_element( bulk_data, hex_part, element_id, node_ids);

  }

  // declare right element on last process
  if (p_rank == p_size -1)
  {
    EntityId element_id = 2;
    stk::mesh::EntityIdVector node_ids { 2, 9, 10, 3, 13, 11, 12, 7};

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

