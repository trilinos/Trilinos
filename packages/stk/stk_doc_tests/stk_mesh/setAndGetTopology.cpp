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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>   // for MeshBuilder
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include <stk_topology/topology.hpp>    // for topology, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityId
namespace stk { namespace mesh { class Part; } }

namespace {

void declare_element_nodes(stk::mesh::BulkData &mesh, stk::mesh::Entity elem1, stk::mesh::Entity elem2)
{
  for(unsigned node_ord = 0; node_ord < 4; ++node_ord)
  {
    stk::mesh::Entity new_node = mesh.declare_node(node_ord + 100 * mesh.identifier(elem1));
    mesh.declare_relation(elem1, new_node, node_ord);
  }

  for(unsigned node_ord = 0; node_ord < 8; ++node_ord)
  {
    stk::mesh::Entity new_node2 = mesh.declare_node(node_ord + 100 * mesh.identifier(elem2));
    mesh.declare_relation(elem2, new_node2, node_ord);
  }
}

//-BEGIN
TEST(stkMeshHowTo, setAndGetTopology)
{
  const unsigned spatialDimension = 3;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDimension);
  builder.set_entity_rank_names(stk::mesh::entity_rank_names());
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::Part &tetPart = metaData.declare_part_with_topology("tet part", stk::topology::TET_4);

  stk::mesh::Part &hexPart = metaData.declare_part("existing part with currently unknown topology");
  // . . . then later assigned
  stk::mesh::set_topology(hexPart, stk::topology::HEX_8);

  metaData.commit();
  stk::mesh::BulkData& bulkData = *bulkPtr;

  bulkData.modification_begin();
  stk::mesh::EntityId elem1Id = 1, elem2Id = 2;
  stk::mesh::Entity elem1 = bulkData.declare_element(elem1Id, stk::mesh::ConstPartVector{&tetPart});
  stk::mesh::Entity elem2 = bulkData.declare_element(elem2Id, stk::mesh::ConstPartVector{&hexPart});
  declare_element_nodes(bulkData, elem1, elem2);
  bulkData.modification_end();

  stk::topology elem1_topology = bulkData.bucket(elem1).topology();
  stk::topology elem2_topology = bulkData.bucket(elem2).topology();

  EXPECT_EQ(stk::topology::TET_4, elem1_topology);
  EXPECT_EQ(stk::topology::HEX_8, elem2_topology);
}
//-END
}
