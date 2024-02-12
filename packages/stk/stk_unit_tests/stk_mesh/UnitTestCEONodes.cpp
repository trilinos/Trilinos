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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE, etc
#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/SkinBoundary.hpp"
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "UnitTestTextMeshFixture.hpp"

namespace
{

void verify_shared_nodes_ownership(const stk::mesh::BulkData& bulk, int expectedOwnerProc)
{
  stk::mesh::for_each_entity_run(bulk, stk::topology::NODE_RANK, bulk.mesh_meta_data().globally_shared_part(),
  [&expectedOwnerProc](const stk::mesh::BulkData& mesh, stk::mesh::Entity entity)
  {
    EXPECT_EQ(expectedOwnerProc, mesh.parallel_owner_rank(entity));
  });
}

void ceo_shared_nodes(stk::mesh::BulkData& bulk)
{
  stk::mesh::EntityProcVec entityProcs;
  const int currentOwnerProc = 0;
  const int newOwnerProc = 1;
  if (bulk.parallel_rank() == 0) {
    stk::mesh::for_each_entity_run(bulk, stk::topology::NODE_RANK, bulk.mesh_meta_data().globally_shared_part(),
    [&currentOwnerProc, &newOwnerProc, &entityProcs](const stk::mesh::BulkData& mesh, stk::mesh::Entity entity)
    {
      EXPECT_EQ(currentOwnerProc, mesh.parallel_owner_rank(entity));
      entityProcs.emplace_back(entity, newOwnerProc);
    });
  }

  bulk.change_entity_owner(entityProcs);
}

class TestCEONodesAura2d : public TestTextMeshAura2d
{
public:
  TestCEONodesAura2d() {}
};

TEST_F(TestCEONodesAura2d, quadSharedNodes)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "0, 1, QUAD_4_2D, 1,2,5,4\n"
                         "0, 2, QUAD_4_2D, 2,3,6,5\n"
                         "1, 3, QUAD_4_2D, 4,5,8,7\n"
                         "1, 4, QUAD_4_2D, 5,6,9,8\n";

  setup_text_mesh(meshDesc);

  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part());

  int ownerProc = 0;
  verify_shared_nodes_ownership(get_bulk(), ownerProc);

  ceo_shared_nodes(get_bulk());

  ownerProc = 1;
  verify_shared_nodes_ownership(get_bulk(), ownerProc);
}

class TestCEONodesAura : public TestTextMeshAura
{
public:
  TestCEONodesAura() {}
};

TEST_F(TestCEONodesAura, hexSharedNodes)
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  std::string meshDesc = "generated:1x1x2";

  stk::io::fill_mesh(meshDesc, get_bulk());

  stk::mesh::create_edges(get_bulk());

  int ownerProc = 0;
  verify_shared_nodes_ownership(get_bulk(), ownerProc);

  ceo_shared_nodes(get_bulk());

  ownerProc = 1;
  verify_shared_nodes_ownership(get_bulk(), ownerProc);
}

} // empty namespace

