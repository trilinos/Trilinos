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

#include <gtest/gtest.h>
#include <stk_ngp_test/ngp_test.hpp>
#include <string>
#include <vector>
#include "ngp/NgpDeviceCommMapUpdateTestUtils.hpp"

#ifdef STK_USE_DEVICE_MESH

NGP_TEST_F(NgpDeviceCommMapUpdate, changeParts_sharedNodes_twoProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  setup_hex_mesh(twoHexMesh);
  capture_nodes_and_gold_values();

  const int peerProc = 1 - bulkPtr->parallel_rank();

  change_parts_on_device(get_owned_shared_nodes(), {newPart->mesh_meta_data_ordinal()});

  check_node_comm_map_keys(peerProc, {5, 6, 7, 8}, {});
  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, changeParts_sharedNodes_threeProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  setup_hex_mesh(threeHexMesh);
  capture_nodes_and_gold_values();

  const int myRank = bulkPtr->parallel_rank();

  change_parts_on_device(get_owned_shared_nodes(), {newPart->mesh_meta_data_ordinal()});

  if (myRank == 0) {
    check_node_comm_map_keys(1, {5, 6, 7, 8}, {});
    check_node_comm_map_keys(2, {}, {});
  }
  else if (myRank == 1) {
    check_node_comm_map_keys(0, {5, 6, 7, 8}, {});
    check_node_comm_map_keys(2, {9, 10, 11, 12}, {});
  }
  else {
    check_node_comm_map_keys(0, {}, {});
    check_node_comm_map_keys(1, {9, 10, 11, 12}, {});
  }

  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, changeParts_sharedNodes_withAura_threeProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  setup_hex_mesh(threeHexMesh, stk::mesh::BulkData::AUTO_AURA);
  capture_nodes_and_gold_values();

  // With aura, elements are ghosted too, so the element comm maps have entries as well as the node ones, and
  // the node comm maps have ghost entries alongside their shared ones.
  change_parts_on_device(get_owned_shared_nodes(), {newPart->mesh_meta_data_ordinal()});

  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, changeParts_customGhosting_symmInfo_threeProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  setup_hex_mesh(threeHexMesh);

  // Node 1 is on P0 only and is custom-ghosted to both P1 and P2.  With symmetric ghost info, P1 and P2 each
  // get a SYMM_INFO entry naming the other, so node 1 is a ghost entry in every pair's comm maps.
  stk::mesh::unit_test::proc0_ghost_node_to_proc1_and_proc2(*bulkPtr, 1);
  capture_nodes_and_gold_values();

  const int myRank = bulkPtr->parallel_rank();

  change_parts_on_device(get_owned_shared_nodes(), {newPart->mesh_meta_data_ordinal()});

  if (myRank == 0) {
    check_node_comm_map_keys(1, {5, 6, 7, 8}, {1});
    check_node_comm_map_keys(2, {}, {1});
  }
  else if (myRank == 1) {
    check_node_comm_map_keys(0, {5, 6, 7, 8}, {1});
    check_node_comm_map_keys(2, {9, 10, 11, 12}, {1});
  }
  else {
    check_node_comm_map_keys(0, {}, {1});
    check_node_comm_map_keys(1, {9, 10, 11, 12}, {1});
  }

  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, changeParts_sharedAndGhostedNode_threeProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  setup_hex_mesh(threeHexMesh);

  // Node 5 is shared between P0 and P1 and also custom-ghosted from P0 to P2.  With symmetric ghost info,
  // P1 gets a SYMM_INFO entry naming P2 and P2 one naming P1, so node 5 is a shared entry in P0's and P1's
  // comm maps for each other and a ghost entry in every other pair's, and the refresh must keep all of them
  // in step.
  stk::mesh::unit_test::proc0_ghost_node5_to_proc2(*bulkPtr);
  capture_nodes_and_gold_values();

  const int myRank = bulkPtr->parallel_rank();

  change_parts_on_device(get_owned_shared_nodes(), {newPart->mesh_meta_data_ordinal()});

  if (myRank == 0) {
    check_node_comm_map_keys(1, {5, 6, 7, 8}, {});
    check_node_comm_map_keys(2, {}, {5});
  }
  else if (myRank == 1) {
    check_node_comm_map_keys(0, {5, 6, 7, 8}, {});
    check_node_comm_map_keys(2, {9, 10, 11, 12}, {5});
  }
  else {
    check_node_comm_map_keys(0, {}, {5});
    check_node_comm_map_keys(1, {9, 10, 11, 12}, {5});
  }

  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, declareEntities_addsLocalNodes_twoProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  setup_axis_hex_mesh(3);
  capture_nodes_and_gold_values();

  const int myRank = bulkPtr->parallel_rank();

  // Different ids on each rank, so the new nodes are local to their rank and nothing new is shared.
  declare_nodes_on_device({static_cast<unsigned>(2000 + myRank), static_cast<unsigned>(3000 + myRank)});

  // Sharing is unchanged, so the map must hold the same entities as before; only their mesh indices move.
  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, destroyRelations_sharedFaces_twoProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  setup_hex_mesh_with_faces(twoHexMesh);
  capture_nodes_and_gold_values();
  capture_faces_and_gold_values();

  const int peerProc = 1 - bulkPtr->parallel_rank();
  const std::vector<stk::mesh::EntityId> sharedFaceIds = get_shared_face_ids();
  ASSERT_EQ(1u, sharedFaceIds.size());

  // batch_destroy_relations is its own entry point and must refresh the map itself.  Dropping the
  // element-to-face relations moves the faces to other buckets but leaves sharing alone.
  destroy_relations_on_device(get_owned_elements(), stk::topology::FACE_RANK);

  check_comm_map_keys(stk::topology::FACE_RANK, peerProc, sharedFaceIds, {});
  check_node_comm_map_keys(peerProc, {5, 6, 7, 8}, {});
  check_device_comm_map_unchanged();
  check_face_parallel_sum_on_device();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, changeParts_removePart_twoProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  setup_hex_mesh(twoHexMesh);
  capture_nodes_and_gold_values();

  const int peerProc = 1 - bulkPtr->parallel_rank();
  const std::vector<stk::mesh::Entity> sharedNodes = get_owned_shared_nodes();

  change_parts_on_device(sharedNodes, {newPart->mesh_meta_data_ordinal()});

  // Then remove the part again, which takes batch_change_entity_parts through its part-removal path.
  // newNodePart has no rank, so no removal can be denied and the denied-part-removal exchange is skipped.
  auto entityView = make_device_view<DeviceEntityView>("entities", sharedNodes);
  DevicePartOrdinalView addPartView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "addPartOrdinals"), 0);
  auto removePartView =
      make_device_view<DevicePartOrdinalView>("removePartOrdinals",
                                              std::vector<stk::mesh::PartOrdinal>{
                                                  newPart->mesh_meta_data_ordinal()});
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
  ngpMesh.batch_change_entity_parts(entityView, addPartView, removePartView);

  check_node_comm_map_keys(peerProc, {5, 6, 7, 8}, {});
  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, changeParts_rankWithNoMeshAtAll_fourProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4) { GTEST_SKIP(); }

  // threeHexMesh puts elements on procs 0-2 only, so P3 has no mesh and an empty comm map.  It must still
  // take part in every collective in the part change.
  setup_hex_mesh(threeHexMesh);
  capture_nodes_and_gold_values();

  change_parts_on_device(get_owned_shared_nodes(), {newPart->mesh_meta_data_ordinal()});

  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, changeParts_onlyOneProcHasWork_twoProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  setup_hex_mesh(twoHexMesh);
  capture_nodes_and_gold_values();

  const int myRank = bulkPtr->parallel_rank();
  const int peerProc = 1 - myRank;

  std::vector<stk::mesh::Entity> entities;
  std::vector<stk::mesh::PartOrdinal> addPartOrdinals;
  if (myRank == 0) {
    entities = get_owned_shared_nodes();
    addPartOrdinals.push_back(newPart->mesh_meta_data_ordinal());
  }
  change_parts_on_device(entities, addPartOrdinals);

  check_node_comm_map_keys(peerProc, {5, 6, 7, 8}, {});
  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, destroyElements_sharedNodesSurvive_twoProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  setup_hex_mesh(twoHexMesh);
  capture_nodes_and_gold_values();

  const int peerProc = 1 - bulkPtr->parallel_rank();

  // Without aura no element is shared or ghosted, so destroying them is allowed and leaves the map's
  // entities unchanged.  It can still move the elements' nodes to other buckets, so the node map's mesh
  // indices must be refreshed.
  destroy_on_device(get_owned_elements());

  check_node_comm_map_keys(peerProc, {5, 6, 7, 8}, {});
  check_comm_map_keys(stk::topology::ELEM_RANK, peerProc, {}, {});
  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, destroyThenRedeclareRelations_twoProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  setup_hex_mesh_with_faces(twoHexMesh);
  capture_nodes_and_gold_values();
  capture_faces_and_gold_values();

  const std::vector<stk::mesh::Entity> elems = get_owned_elements();
  const RelationsToDeclare elemFaces = capture_relations(elems, stk::topology::FACE_RANK);

  destroy_relations_on_device(elems, stk::topology::FACE_RANK);
  check_device_comm_map_unchanged();

  // batch_declare_relations is also its own entry point and must refresh the map too.  Putting the relations
  // back gives the faces block_1 membership again, moving them between buckets once more.
  declare_relations_on_device(elemFaces);

  check_device_comm_map_unchanged();
  check_face_parallel_sum_on_device();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, destroyThenRedeclareRelations_withAura_threeProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  setup_hex_mesh_with_faces(threeHexMesh, stk::mesh::BulkData::AUTO_AURA);
  capture_nodes_and_gold_values();
  capture_faces_and_gold_values();

  const std::vector<stk::mesh::Entity> elems = get_owned_elements();
  const RelationsToDeclare elemFaces = capture_relations(elems, stk::topology::FACE_RANK);

  // With aura, the face and element comm maps have both shared and ghost entries for several procs, and
  // all of them must keep their entities through two modifications that move faces between buckets.
  destroy_relations_on_device(elems, stk::topology::FACE_RANK);
  check_device_comm_map_unchanged();

  declare_relations_on_device(elemFaces);

  check_device_comm_map_unchanged();
  check_face_parallel_sum_on_device();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, meshIndicesRefresh_partChangeSplitsSharedNodeBucket_twoProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  setup_hex_mesh(twoHexMesh);
  capture_nodes_and_gold_values();

  const int peerProc = 1 - bulkPtr->parallel_rank();
  const std::vector<stk::mesh::Entity> sharedNodes = get_owned_shared_nodes();

  // A first device modification, so the indices below are ones the device has already refreshed.
  change_parts_on_device(sharedNodes, {newPart->mesh_meta_data_ordinal()});
  const std::vector<stk::mesh::FastMeshIndex> before =
      device_comm_map_indices(stk::mesh::get_updated_ngp_mesh(*bulkPtr), stk::topology::NODE_RANK, peerProc);

  // Then put only half of the shared nodes in another part.  That splits their bucket, so the map's mesh
  // indices must change while its entities and their order must not.
  std::vector<stk::mesh::Entity> half;
  if (!sharedNodes.empty()) { half.assign(sharedNodes.begin(), sharedNodes.begin() + sharedNodes.size() / 2); }
  change_parts_on_device(half, {splitPart->mesh_meta_data_ordinal()});
  const std::vector<stk::mesh::FastMeshIndex> after =
      device_comm_map_indices(stk::mesh::get_updated_ngp_mesh(*bulkPtr), stk::topology::NODE_RANK, peerProc);

  ASSERT_EQ(before.size(), after.size());
  bool anyMoved = false;
  for (size_t i = 0; i < before.size(); ++i) {
    anyMoved |= before[i].bucket_id != after[i].bucket_id || before[i].bucket_ord != after[i].bucket_ord;
  }
  // P0 owns the shared nodes and moves half of them, and P1 receives the same part change for its copies.
  EXPECT_TRUE(anyMoved) << "P" << bulkPtr->parallel_rank() << ": no comm entity moved, so this test cannot "
                        << "distinguish refreshed mesh indices from stale ones";

  check_node_comm_map_keys(peerProc, {5, 6, 7, 8}, {});
  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, meshIndicesRefresh_relationChangesMoveSharedFaces_withAura_threeProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  setup_hex_mesh_with_faces(threeHexMesh, stk::mesh::BulkData::AUTO_AURA);
  capture_nodes_and_gold_values();
  capture_faces_and_gold_values();

  const std::vector<stk::mesh::Entity> elems = get_owned_elements();
  const RelationsToDeclare elemFaces = capture_relations(elems, stk::topology::FACE_RANK);
  // A host BulkData query, so make it before the device modifications.
  const std::vector<stk::mesh::Entity> ownedShared = get_owned_shared_nodes();

  // Detaching the faces from their elements drops their induced block_1 membership, moving every face --
  // shared and aura-ghosted alike -- to another bucket; redeclaring the relations moves them back.  Each of
  // these modifications, and the part change after them, must refresh the map's mesh indices without
  // changing which entities it holds.
  destroy_relations_on_device(elems, stk::topology::FACE_RANK);
  declare_relations_on_device(elemFaces);
  change_parts_on_device(ownedShared, {newPart->mesh_meta_data_ordinal()});

  check_device_comm_map_unchanged();
  check_face_parallel_sum_on_device();
  check_parallel_sum_on_device();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, changeParts_customGhosting_noSymmInfo_threeProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) {
    GTEST_SKIP();
  }

  useSymmetricGhostInfo = false;
  setup_hex_mesh(threeHexMesh);
  stk::mesh::unit_test::proc0_ghost_node_to_proc1_and_proc2(*bulkPtr, 1);
  capture_nodes_and_gold_values();

  const int myRank = bulkPtr->parallel_rank();
  change_parts_on_device(get_owned_shared_nodes(), {newPart->mesh_meta_data_ordinal()});

  // Unlike changeParts_customGhosting_symmInfo_threeProcs, P1 and P2 have no ghost entry for each other.
  if (myRank == 0) {
    check_node_comm_map_keys(1, {5, 6, 7, 8}, {1});
    check_node_comm_map_keys(2, {}, {1});
  } else if (myRank == 1) {
    check_node_comm_map_keys(0, {5, 6, 7, 8}, {1});
    check_node_comm_map_keys(2, {9, 10, 11, 12}, {});
  } else {
    check_node_comm_map_keys(0, {}, {1});
    check_node_comm_map_keys(1, {9, 10, 11, 12}, {});
  }

  check_device_comm_map_unchanged();
}

NGP_TEST_F(NgpDeviceCommMapUpdate, destroyAuraGhostedElement_throws_withAura_twoProcs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  setup_hex_mesh(twoHexMesh, stk::mesh::BulkData::AUTO_AURA);
  capture_nodes_and_gold_values();

  const int peerProc = 1 - bulkPtr->parallel_rank();

  // Each element is aura-ghosted to the other proc, so destroying it on device would change the comm map on
  // one side only.  That is refused before anything is modified, so the map and field exchange still work.
  check_comm_map_keys(stk::topology::ELEM_RANK, peerProc, {}, {1, 2});

  EXPECT_ANY_THROW(destroy_on_device(get_owned_elements()));

  check_device_comm_map_unchanged();
  check_parallel_sum_on_device();
}

#endif  // STK_USE_DEVICE_MESH
