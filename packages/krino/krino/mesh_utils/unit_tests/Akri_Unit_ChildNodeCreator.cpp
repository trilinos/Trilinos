// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ChildNodeCreator.hpp>
#include <Akri_EntityIdPool.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>

#include <array>
#include <vector>

namespace krino
{

// TwoRightTets shares nodes 0, 2 and 4 between the two tets, so distributing
// tet 1 to proc 0 and tet 2 to proc 1 makes those nodes shared. Child nodes
// whose parents are shared must be created consistently on both procs.
class TwoTetChildNodeCreation : public StkMeshTetFixture
{
protected:
  TwoTetChildNodeCreation()
  {
    set_valid_proc_sizes_for_test({1, 2});
    if (is_valid_proc_size_for_test())
    {
      TwoRightTets meshSpec;
      // Put both tets on proc 0 in serial; split one tet per proc at np2 so the
      // nodes shared by both tets (0, 2, 4) become parallel-shared nodes.
      const std::vector<int> procOwners = (parallel_size() == 1) ? std::vector<int>{0, 0} : std::vector<int>{0, 1};
      StkMeshTetFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1, 1}, procOwners);
    }
  }

  // Mirrors the id-generation lambda used by the production NodeRefiner caller.
  GenerateNewEntityIdsFunction make_id_generator()
  {
    stk::mesh::BulkData & mesh = get_mesh();
    return [&mesh](stk::topology::rank_t rank, size_t numIdsNeeded, std::vector<stk::mesh::EntityId> & ids)
    {
      EntityIdPool::generate_new_ids(mesh, rank, numIdsNeeded, ids, false, false);
    };
  }
};

//--------------------------------------------------------------------------
// A single child-node request creates one valid node with a fresh id that is
// distinct from its parents.
//--------------------------------------------------------------------------
TEST_F(TwoTetChildNodeCreation, singleRequest_createsValidChildNodeWithFreshId)
{
  if (!is_valid_proc_size_for_test())
    return;

  stk::mesh::BulkData & mesh = get_mesh();

  // Nodes 0 and 2 are shared by both tets (present on every proc that owns a tet).
  std::array<stk::mesh::Entity, 2> parents{ get_assigned_node_for_index(0), get_assigned_node_for_index(2) };
  ASSERT_TRUE(mesh.is_valid(parents[0]));
  ASSERT_TRUE(mesh.is_valid(parents[1]));
  const stk::mesh::EntityId parentId0 = mesh.identifier(parents[0]);
  const stk::mesh::EntityId parentId1 = mesh.identifier(parents[1]);

  stk::mesh::Entity childNode;
  std::vector<ChildNodeRequest> requests;
  requests.emplace_back(std::vector<stk::mesh::Entity*>{&parents[0], &parents[1]}, &childNode);

  auto generate_new_ids = make_id_generator();

  mesh.modification_begin();
  batch_create_child_nodes(mesh, requests, {}, generate_new_ids);
  mesh.modification_end();

  ASSERT_TRUE(mesh.is_valid(childNode));
  const stk::mesh::EntityId childId = mesh.identifier(childNode);
  EXPECT_GT(childId, 0u);
  EXPECT_NE(childId, parentId0);
  EXPECT_NE(childId, parentId1);
  EXPECT_EQ(stk::topology::NODE_RANK, mesh.entity_rank(childNode));
}

//--------------------------------------------------------------------------
// Multiple requests in one batch each create a distinct node.
//--------------------------------------------------------------------------
TEST_F(TwoTetChildNodeCreation, multipleRequests_createDistinctChildNodes)
{
  if (!is_valid_proc_size_for_test())
    return;

  stk::mesh::BulkData & mesh = get_mesh();

  std::array<stk::mesh::Entity, 2> parentsA{ get_assigned_node_for_index(0), get_assigned_node_for_index(2) };
  std::array<stk::mesh::Entity, 2> parentsB{ get_assigned_node_for_index(2), get_assigned_node_for_index(4) };

  stk::mesh::Entity childA;
  stk::mesh::Entity childB;
  std::vector<ChildNodeRequest> requests;
  requests.emplace_back(std::vector<stk::mesh::Entity*>{&parentsA[0], &parentsA[1]}, &childA);
  requests.emplace_back(std::vector<stk::mesh::Entity*>{&parentsB[0], &parentsB[1]}, &childB);

  auto generate_new_ids = make_id_generator();

  mesh.modification_begin();
  batch_create_child_nodes(mesh, requests, {}, generate_new_ids);
  mesh.modification_end();

  ASSERT_TRUE(mesh.is_valid(childA));
  ASSERT_TRUE(mesh.is_valid(childB));
  EXPECT_NE(mesh.identifier(childA), mesh.identifier(childB));
}

//--------------------------------------------------------------------------
// When both parents are shared across procs, the child node must be created
// with the SAME id on every sharing proc and marked shared. This is the whole
// point of the parallel negotiation in batch_create_child_nodes.
//--------------------------------------------------------------------------
TEST_F(TwoTetChildNodeCreation, sharedParents_createConsistentSharedChildNode)
{
  if (!is_valid_proc_size_for_test())
    return;

  stk::mesh::BulkData & mesh = get_mesh();

  // Parents 0 and 2 are both shared between proc 0 and proc 1 at np2.
  std::array<stk::mesh::Entity, 2> parents{ get_assigned_node_for_index(0), get_assigned_node_for_index(2) };
  ASSERT_TRUE(mesh.is_valid(parents[0]));
  ASSERT_TRUE(mesh.is_valid(parents[1]));

  stk::mesh::Entity childNode;
  std::vector<ChildNodeRequest> requests;
  requests.emplace_back(std::vector<stk::mesh::Entity*>{&parents[0], &parents[1]}, &childNode);

  auto generate_new_ids = make_id_generator();

  mesh.modification_begin();
  batch_create_child_nodes(mesh, requests, {}, generate_new_ids);
  mesh.modification_end();

  ASSERT_TRUE(mesh.is_valid(childNode));

  // The identifier must agree across all procs (min == max).
  const uint64_t localId = mesh.identifier(childNode);
  uint64_t globalMinId = 0;
  uint64_t globalMaxId = 0;
  stk::all_reduce_min(mesh.parallel(), &localId, &globalMinId, 1);
  stk::all_reduce_max(mesh.parallel(), &localId, &globalMaxId, 1);
  EXPECT_EQ(globalMinId, globalMaxId) << "child node id must be identical on all sharing procs";
  EXPECT_EQ(localId, globalMinId);

  if (mesh.parallel_size() > 1)
  {
    EXPECT_TRUE(mesh.bucket(childNode).shared()) << "a child of shared parents must itself be shared";
  }
}

}
