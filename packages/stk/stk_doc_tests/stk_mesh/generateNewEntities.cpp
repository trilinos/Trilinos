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
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include <stk_topology/topology.hpp>    // for topology, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityId
namespace stk { namespace mesh { class Part; } }

namespace {

void check_connectivities_for_stkMeshHowTo_generateNewEntities(stk::mesh::BulkData &mesh,
                                                               stk::mesh::Entity tet_elem, stk::mesh::Entity hex_elem,
                                                               const std::vector<stk::mesh::Entity> &generated_entities)
{
  // Check that entities are connected as expected.
  unsigned node_i = 0;
  EXPECT_EQ(4u, mesh.num_nodes(tet_elem));
  stk::mesh::Entity const *node_rels = mesh.begin_nodes(tet_elem);
  for(unsigned node_ord = 0 ; node_ord < 4; ++node_ord, ++node_i)
  {
    EXPECT_TRUE(mesh.is_valid(node_rels[node_ord]));
    EXPECT_EQ(generated_entities[node_i], node_rels[node_ord]);
  }
  EXPECT_EQ(8u, mesh.num_nodes(hex_elem));
  node_rels = mesh.begin_nodes(hex_elem);
  for(unsigned node_ord = 0 ; node_ord < 8; ++node_ord, ++node_i)
  {
    EXPECT_TRUE(mesh.is_valid(node_rels[node_ord]));
    EXPECT_EQ(generated_entities[node_i], node_rels[node_ord]);
  }
}

//-BEGIN
TEST(stkMeshHowTo, generateNewEntities)
{
  const unsigned spatialDimension = 3;

  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDimension);
  builder.set_entity_rank_names(stk::mesh::entity_rank_names());
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::Part &tetPart = metaData.declare_part_with_topology("tetElementPart", stk::topology::TET_4);
  stk::mesh::Part &hexPart = metaData.declare_part_with_topology("hexElementPart", stk::topology::HEX_8);
  metaData.commit();

  // Parts vectors handy for setting topology later.
  std::vector<stk::mesh::Part *> add_tetPart(1);
  add_tetPart[0] = &tetPart;
  std::vector<stk::mesh::Part *> add_hexPart(1);
  add_hexPart[0] = &hexPart;

  stk::mesh::BulkData& mesh = *bulkPtr;
  mesh.modification_begin();

  std::vector<size_t> requests(metaData.entity_rank_count(), 0);
  const size_t num_nodes_requested = 12;
  const size_t num_elems_requested =  2;
  requests[stk::topology::NODE_RANK] = num_nodes_requested;
  requests[stk::topology::ELEMENT_RANK] = num_elems_requested;
  std::vector<stk::mesh::Entity> requested_entities;

  mesh.generate_new_entities(requests, requested_entities);

  // Set topologies of new entities with rank > stk::topology::NODE_RANK.
  stk::mesh::Entity elem1 = requested_entities[num_nodes_requested];
  mesh.change_entity_parts(elem1, add_tetPart);
  stk::mesh::Entity elem2 = requested_entities[num_nodes_requested + 1];
  mesh.change_entity_parts(elem2, add_hexPart);

  // Set downward relations of entities with rank > stk::topology::NODE_RANK
  unsigned node_i = 0;
  for(unsigned node_ord = 0 ; node_ord < 4; ++node_ord, ++node_i)
  {
    mesh.declare_relation( elem1 , requested_entities[node_i] , node_ord);
  }
  for(unsigned node_ord = 0 ; node_ord < 8; ++node_ord, ++node_i)
  {
    mesh.declare_relation( elem2 , requested_entities[node_i] , node_ord);
  }
  mesh.modification_end();

  check_connectivities_for_stkMeshHowTo_generateNewEntities(mesh, elem1, elem2, requested_entities);

  // Not setting topologies of new entities with rank > stk::topology::NODE_RANK causes throw
  mesh.modification_begin();
  std::vector<stk::mesh::Entity> more_requested_entities;
  mesh.generate_new_entities(requests, more_requested_entities);
#ifdef NDEBUG
  mesh.modification_end();
#else
  EXPECT_THROW(mesh.modification_end(), std::logic_error);
#endif
}
//-END
}
