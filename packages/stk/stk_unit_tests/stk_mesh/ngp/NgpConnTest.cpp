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

#include <Kokkos_Core.hpp>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/DestroyElements.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/baseImpl/MeshConnectivity.hpp>
#include <stk_mesh/baseImpl/MeshConnUtils.hpp>
#include <stk_io/FillMesh.hpp>

void test_basic_access_on_device()
{
  EXPECT_LT(sizeof(stk::mesh::impl::EntityConnectivity), 90u);
  EXPECT_LT(sizeof(stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace>),106u);

  stk::mesh::impl::EntityConnectivity entConn;

  size_t totalAlloc = 1024, minAlloc = 32, maxAlloc = 64;
  stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> meshConnDevice(10, totalAlloc, minAlloc, maxAlloc);
  stk::mesh::Entity ent1(1);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {
      NGP_EXPECT_EQ(0u, entConn.get_capacity());
      auto conn = meshConnDevice.get_connected_entities(ent1,stk::topology::NODE_RANK);
      auto conn2 = entConn.get_connected_entities(stk::topology::NODE_RANK);
      NGP_EXPECT_EQ(0u, conn.size());
      NGP_EXPECT_EQ(0u, conn2.size());
  });
}

NGP_TEST(NgpConnectivity, basic)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  test_basic_access_on_device();
}

void test_add_connectivity_on_host()
{
  size_t totalAlloc = 1024, minAlloc = 32, maxAlloc = 64;
  stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> meshConnHost(2, totalAlloc, minAlloc, maxAlloc);

  stk::mesh::Entity ent1(1), ent2(2), ent3(3);
  stk::mesh::ConnectivityOrdinal ord2(2), ord3(3);
  std::vector<stk::mesh::Entity> ents = { ent2, ent3 };
  std::vector<stk::mesh::ConnectivityOrdinal> ords = { ord2, ord3 };

  meshConnHost.add_connectivity(ent1, stk::topology::FACE_RANK, ents.size(), ents.data(), ords.data(), nullptr, false);
  meshConnHost.add_connectivity(ent1, stk::topology::NODE_RANK, ents.size(), ents.data(), ords.data(), nullptr, false);
  auto elems = meshConnHost.get_connected_entities(ent1, stk::topology::ELEM_RANK);
  auto faces = meshConnHost.get_connected_entities(ent1, stk::topology::FACE_RANK);
  auto edges = meshConnHost.get_connected_entities(ent1, stk::topology::EDGE_RANK);
  auto nodes = meshConnHost.get_connected_entities(ent1, stk::topology::NODE_RANK);

  auto elemOrds = meshConnHost.get_connected_ordinals(ent1, stk::topology::ELEM_RANK);
  auto faceOrds = meshConnHost.get_connected_ordinals(ent1, stk::topology::FACE_RANK);
  auto edgeOrds = meshConnHost.get_connected_ordinals(ent1, stk::topology::EDGE_RANK);
  auto nodeOrds = meshConnHost.get_connected_ordinals(ent1, stk::topology::NODE_RANK);

  auto elemPerms = meshConnHost.get_connected_permutations(ent1, stk::topology::ELEM_RANK);
  auto facePerms = meshConnHost.get_connected_permutations(ent1, stk::topology::FACE_RANK);
  auto edgePerms = meshConnHost.get_connected_permutations(ent1, stk::topology::EDGE_RANK);
  auto nodePerms = meshConnHost.get_connected_permutations(ent1, stk::topology::NODE_RANK);

  EXPECT_EQ(0u, elems.size());
  EXPECT_EQ(2u, faces.size());
  EXPECT_EQ(0u, edges.size());
  EXPECT_EQ(2u, nodes.size());

  EXPECT_EQ(faces[0], ent2);
  EXPECT_EQ(faces[1], ent3);
  EXPECT_EQ(nodes[0], ent2);
  EXPECT_EQ(nodes[1], ent3);

  EXPECT_EQ(0u, elemOrds.size());
  EXPECT_EQ(2u, faceOrds.size());
  EXPECT_EQ(0u, edgeOrds.size());
  EXPECT_EQ(2u, nodeOrds.size());

  EXPECT_EQ(faceOrds[0], ord2);
  EXPECT_EQ(faceOrds[1], ord3);
  EXPECT_EQ(nodeOrds[0], ord2);
  EXPECT_EQ(nodeOrds[1], ord3);

  EXPECT_EQ(0u, elemPerms.size());
  EXPECT_EQ(0u, facePerms.size());
  EXPECT_EQ(0u, edgePerms.size());
  EXPECT_EQ(0u, nodePerms.size());
}

NGP_TEST(NgpConnectivity, add_connectivity_on_host)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  test_add_connectivity_on_host();
}

void test_add_connectivity_with_permutations_on_host()
{
  size_t totalAlloc = 1024, minAlloc = 32, maxAlloc = 64;
  stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> meshConnHost(2, totalAlloc, minAlloc, maxAlloc);

  stk::mesh::Entity ent1(1), ent2(2), ent3(3);
  stk::mesh::ConnectivityOrdinal ord2(2), ord3(3);
  std::vector<stk::mesh::Entity> ents = { ent2, ent3 };
  std::vector<stk::mesh::ConnectivityOrdinal> ords = { ord2, ord3 };
  stk::mesh::Permutation perm2=static_cast<stk::mesh::Permutation>(2), perm3=static_cast<stk::mesh::Permutation>(3);
  std::vector<stk::mesh::Permutation> permVec = { perm2, perm3 };
  const bool addPermutations = true;

  meshConnHost.set_capacity(ent1, 6, addPermutations);
  meshConnHost.add_connectivity(ent1, stk::topology::NODE_RANK, ents.size(), ents.data(), ords.data(), nullptr, addPermutations);
  meshConnHost.add_connectivity(ent1, stk::topology::FACE_RANK, ents.size(), ents.data(), ords.data(), permVec.data(), addPermutations);
  meshConnHost.add_connectivity(ent1, stk::topology::EDGE_RANK, ents.size(), ents.data(), ords.data(), nullptr, addPermutations);

  auto elems = meshConnHost.get_connected_entities(ent1, stk::topology::ELEM_RANK);
  auto faces = meshConnHost.get_connected_entities(ent1, stk::topology::FACE_RANK);
  auto edges = meshConnHost.get_connected_entities(ent1, stk::topology::EDGE_RANK);
  auto nodes = meshConnHost.get_connected_entities(ent1, stk::topology::NODE_RANK);

  auto elemOrds = meshConnHost.get_connected_ordinals(ent1, stk::topology::ELEM_RANK);
  auto faceOrds = meshConnHost.get_connected_ordinals(ent1, stk::topology::FACE_RANK);
  auto edgeOrds = meshConnHost.get_connected_ordinals(ent1, stk::topology::EDGE_RANK);
  auto nodeOrds = meshConnHost.get_connected_ordinals(ent1, stk::topology::NODE_RANK);

  auto elemPerms = meshConnHost.get_connected_permutations(ent1, stk::topology::ELEM_RANK);
  auto facePerms = meshConnHost.get_connected_permutations(ent1, stk::topology::FACE_RANK);
  auto edgePerms = meshConnHost.get_connected_permutations(ent1, stk::topology::EDGE_RANK);
  auto nodePerms = meshConnHost.get_connected_permutations(ent1, stk::topology::NODE_RANK);

  EXPECT_EQ(0u, elems.size());
  EXPECT_EQ(2u, faces.size());
  EXPECT_EQ(2u, edges.size());
  EXPECT_EQ(2u, nodes.size());

  EXPECT_EQ(faces[0], ent2);
  EXPECT_EQ(faces[1], ent3);
  EXPECT_EQ(nodes[0], ent2);
  EXPECT_EQ(nodes[1], ent3);

  EXPECT_EQ(0u, elemOrds.size());
  EXPECT_EQ(2u, faceOrds.size());
  EXPECT_EQ(2u, edgeOrds.size());
  EXPECT_EQ(2u, nodeOrds.size());

  EXPECT_EQ(faceOrds[0], ord2);
  EXPECT_EQ(faceOrds[1], ord3);
  EXPECT_EQ(nodeOrds[0], ord2);
  EXPECT_EQ(nodeOrds[1], ord3);

  EXPECT_EQ(0u, elemPerms.size());
  EXPECT_EQ(2u, facePerms.size());
  EXPECT_EQ(2u, edgePerms.size());
  EXPECT_EQ(2u, nodePerms.size());

  EXPECT_EQ(facePerms[0], perm2);
  EXPECT_EQ(facePerms[1], perm3);
  EXPECT_EQ(edgePerms[0], stk::mesh::INVALID_PERMUTATION);
  EXPECT_EQ(edgePerms[1], stk::mesh::INVALID_PERMUTATION);
  EXPECT_EQ(nodePerms[0], stk::mesh::INVALID_PERMUTATION);
  EXPECT_EQ(nodePerms[1], stk::mesh::INVALID_PERMUTATION);
}

NGP_TEST(NgpConnectivity, add_connectivity_with_permutations_on_host)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  test_add_connectivity_with_permutations_on_host();
}

void test_add_connectivity_access_on_device()
{
  size_t totalAlloc = 1024, minAlloc = 32, maxAlloc = 64;
  stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> meshConn(2, totalAlloc, minAlloc, maxAlloc);

  stk::mesh::Entity ent1(1), ent2(2), ent3(3);
  stk::mesh::ConnectivityOrdinal ord2(2), ord3(3);
  stk::mesh::Entity face4(4), face5(5), face6(6);
  stk::mesh::ConnectivityOrdinal faceOrd4(4), faceOrd5(5), faceOrd6(6);
  std::vector<stk::mesh::Entity> ents = { ent2, ent3 };
  std::vector<stk::mesh::ConnectivityOrdinal> ords = { ord2, ord3 };
  std::vector<stk::mesh::Entity> faceVec = { face4, face5, face6 };
  std::vector<stk::mesh::ConnectivityOrdinal> faceOrdVec = { faceOrd4, faceOrd5, faceOrd6 };

  meshConn.set_capacity(ent1, 5, false);
  meshConn.add_connectivity(ent1, stk::topology::NODE_RANK, ents.size(), ents.data(), ords.data(), nullptr, false);
  meshConn.add_connectivity(ent1, stk::topology::FACE_RANK, faceVec.size(), faceVec.data(), faceOrdVec.data(), nullptr, false);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {
      auto connNodes = meshConn.get_connected_entities(ent1,stk::topology::NODE_RANK);
      auto connFaces = meshConn.get_connected_entities(ent1,stk::topology::FACE_RANK);
      NGP_EXPECT_EQ(2u, connNodes.size());
      NGP_EXPECT_EQ(3u, connFaces.size());
  });
}

NGP_TEST(NgpConnectivity, add_connectivity_access_on_device)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  test_add_connectivity_access_on_device();
}

void test_fill_connectivity_from_bulkdata(const stk::mesh::BulkData& bulk)
{
  stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> meshConn;
  stk::mesh::impl::fill_mesh_connectivity(bulk, meshConn);

  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  stk::mesh::Entity face16 = bulk.get_entity(stk::topology::FACE_RANK, 16);

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  unsigned numElems = stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK, meta.universal_part());
  if (numElems == 1) {
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {
    NGP_EXPECT_EQ(8u, meshConn.get_connected_entities(elem1, stk::topology::NODE_RANK).size());
    NGP_EXPECT_EQ(6u, meshConn.get_connected_entities(elem1, stk::topology::FACE_RANK).size());
    NGP_EXPECT_EQ(4u, meshConn.get_connected_entities(face16, stk::topology::NODE_RANK).size());
    NGP_EXPECT_EQ(1u, meshConn.get_connected_entities(face16, stk::topology::ELEM_RANK).size());
    NGP_EXPECT_EQ(0u, meshConn.get_connected_entities(node1, stk::topology::NODE_RANK).size());
    NGP_EXPECT_EQ(3u, meshConn.get_connected_entities(node1, stk::topology::FACE_RANK).size());
  });
  }
}

NGP_TEST(NgpConnectivity, fill_from_bulkdata_1_hex)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:1x1x1|sideset:xXyYzZ", *bulkPtr);

  test_fill_connectivity_from_bulkdata(*bulkPtr);
}

NGP_TEST(NgpConnectivity, fill_from_bulkdata_6_tets)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:1x1x1|tets|sideset:xXyYzZ", *bulkPtr);

  test_fill_connectivity_from_bulkdata(*bulkPtr);
}

NGP_TEST(NgpConnectivity, fill_from_bulkdata_4_hexes)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:2x2x1|sideset:xXyYzZ", *bulkPtr);

  test_fill_connectivity_from_bulkdata(*bulkPtr);
}

NGP_TEST(NgpConnectivity, fill_from_bulkdata_1000_hexes)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:10x10x10|sideset:xXyYzZ", *bulkPtr);

  test_fill_connectivity_from_bulkdata(*bulkPtr);
}

NGP_TEST(NgpConnectivity, fill_from_bulkdata_27000_hexes)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:30x30x30|sideset:xXyYzZ", *bulkPtr);

  test_fill_connectivity_from_bulkdata(*bulkPtr);
}

NGP_TEST(NgpConnectivity, fill_from_bulkdata_6000_tets)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:10x10x10|tets|sideset:xXyYzZ", *bulkPtr);

  test_fill_connectivity_from_bulkdata(*bulkPtr);
}

void create_elem_side1(stk::mesh::BulkData& mesh)
{
  mesh.modification_begin();

  stk::mesh::Entity elem1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  unsigned zeroBasedSideOrdinal = 0;
  mesh.declare_element_side<stk::mesh::PartVector>(elem1, zeroBasedSideOrdinal);

  mesh.modification_end();
}

void delete_elem2(stk::mesh::BulkData& mesh)
{
  stk::mesh::Entity elem2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  stk::mesh::EntityVector elems = {elem2};
  stk::mesh::destroy_elements(mesh, elems);
  const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  EXPECT_EQ(1u, stk::mesh::count_entities(mesh, stk::topology::ELEM_RANK, meta.universal_part()));
  EXPECT_EQ(8u, stk::mesh::count_entities(mesh, stk::topology::NODE_RANK, meta.universal_part()));
}

void change_parts_elem2(stk::mesh::BulkData& mesh)
{
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  stk::mesh::Part& newPart = meta.declare_part("newElemPart", stk::topology::ELEM_RANK);
  stk::mesh::PartVector addParts = {&newPart};
  stk::mesh::PartVector noRemoveParts;
  stk::mesh::Entity elem2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  stk::mesh::EntityVector elems = {elem2};
  mesh.batch_change_entity_parts(elems, addParts, noRemoveParts);
}

NGP_TEST(NgpConnectivity, are_enough_blocks_available)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::vector<std::pair<size_t,size_t>> needed = {{16,4},{32,8},{64,2}};
  std::vector<std::pair<size_t,size_t>> avail_enough = {{32,12}, {64,2}};
  EXPECT_TRUE(stk::mesh::impl::blocks_are_available(needed, avail_enough));

  needed = {{16,4},{32,8},{64,2}};
  avail_enough = {{16,2}, {32,10}, {64,2}};
  EXPECT_TRUE(stk::mesh::impl::blocks_are_available(needed, avail_enough));

  needed = {{16,4},{32,8},{64,2}};
  std::vector<std::pair<size_t,size_t>> avail_not_enough = {{32,10}, {64,2}};
  EXPECT_FALSE(stk::mesh::impl::blocks_are_available(needed, avail_not_enough));

  needed = {{16,4},{32,8},{64,2}};
  avail_not_enough = {{16,4}, {32,8}, {64,1}};
  EXPECT_FALSE(stk::mesh::impl::blocks_are_available(needed, avail_not_enough));
}

NGP_TEST(NgpConnectivity, blocks_are_available_update_mesh_connectivity)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:1x1x1", *bulkPtr);

  stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> meshConn;
  stk::mesh::impl::fill_mesh_connectivity(*bulkPtr, meshConn);

  create_elem_side1(*bulkPtr);

  auto [allocBlocksNeeded, modifiedEntities] = stk::mesh::impl::dealloc_and_get_needed_allocs_using_entity_states(*bulkPtr, meshConn);

  EXPECT_TRUE(!allocBlocksNeeded.empty());
  {
    std::vector<std::pair<size_t,size_t>> gold = {{16,4},{64,2}};
    EXPECT_EQ(gold, allocBlocksNeeded);
  }
  EXPECT_EQ(6u, modifiedEntities.size());

  //empty mod cycle
  bulkPtr->modification_begin();
  bulkPtr->modification_end();

  //the empty mod cycle wiped out the mesh entity-states indicating what exactly
  //had been modified when 'create_elem_side1' was called.
  //So a call to 'dealloc_and_get_needed_allocs_using_entity_states' now, will
  //return incomplete information.
  auto [allocBlocksNeeded_wrongStates, modifiedEntities_wrongStates] = stk::mesh::impl::dealloc_and_get_needed_allocs_using_entity_states(*bulkPtr, meshConn);

  //but a call to 'dealloc_and_get_needed_allocs_full_mesh' will return the correct
  //information because it ignores entity-state info and checks the full mesh.
  auto [allocBlocksNeeded_fullMesh, modifiedEntities_fullMesh] = stk::mesh::impl::dealloc_and_get_needed_allocs_full_mesh(*bulkPtr, meshConn);

  EXPECT_NE(allocBlocksNeeded_wrongStates, allocBlocksNeeded);
  EXPECT_EQ(allocBlocksNeeded_fullMesh, allocBlocksNeeded);
  EXPECT_NE(modifiedEntities_wrongStates, modifiedEntities);
  EXPECT_EQ(modifiedEntities_fullMesh, modifiedEntities);

  {
    auto allocBlocksAvail = stk::mesh::impl::get_available_blocks(meshConn.get_memory_pool());
    EXPECT_TRUE(stk::mesh::impl::blocks_are_available(allocBlocksNeeded, allocBlocksAvail));
  }

  EXPECT_NO_THROW(stk::mesh::impl::update_mesh_connectivity(*bulkPtr, meshConn, modifiedEntities));

  {
    auto allocBlocksAvail = stk::mesh::impl::get_available_blocks(meshConn.get_memory_pool());
    std::vector<std::pair<size_t,size_t>> gold = {{16,328},{32,168},{64,86}};
    EXPECT_EQ(gold, allocBlocksAvail);
  }
}

NGP_TEST(NgpConnectivity, delete_entities_update_mesh_connectivity)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:1x1x2", *bulkPtr);

  stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> meshConn;
  stk::mesh::impl::fill_mesh_connectivity(*bulkPtr, meshConn);

  delete_elem2(*bulkPtr);

  auto [allocBlocksNeeded, modifiedEntities] = stk::mesh::impl::dealloc_and_get_needed_allocs_full_mesh(*bulkPtr, meshConn);
  EXPECT_TRUE(!allocBlocksNeeded.empty());
  {
    std::vector<std::pair<size_t,size_t>> gold = {{8,4}};
    EXPECT_EQ(gold, allocBlocksNeeded);
  }
  constexpr unsigned fourNodesLostElement = 4;
  EXPECT_EQ(fourNodesLostElement, modifiedEntities.size());

  auto allocBlocksAvail = stk::mesh::impl::get_available_blocks(meshConn.get_memory_pool());
  std::vector<std::pair<size_t,size_t>> gold = { {16, 332}, {32, 168}, {64, 87} };
  EXPECT_EQ(gold, allocBlocksAvail);

  EXPECT_TRUE(stk::mesh::impl::blocks_are_available(allocBlocksNeeded, allocBlocksAvail));

  EXPECT_NO_THROW(stk::mesh::impl::update_mesh_connectivity(*bulkPtr, meshConn, modifiedEntities));

  allocBlocksAvail = stk::mesh::impl::get_available_blocks(meshConn.get_memory_pool());
  gold = {{16,328},{32,168},{64,87}};
  EXPECT_EQ(gold, allocBlocksAvail);
}

NGP_TEST(NgpConnectivity, change_parts_connectivity_doesnt_need_updating)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:1x1x2", *bulkPtr);

  stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> meshConn;
  stk::mesh::impl::fill_mesh_connectivity(*bulkPtr, meshConn);

  change_parts_elem2(*bulkPtr);

  auto [allocBlocksNeeded, modifiedEntities] = stk::mesh::impl::dealloc_and_get_needed_allocs_full_mesh(*bulkPtr, meshConn);
  EXPECT_TRUE(allocBlocksNeeded.empty());
}

NGP_TEST(NgpConnectivity, blocks_not_available_cannot_update_mesh_connectivity)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:2x2x1", *bulkPtr);

  stk::mesh::impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> meshConn;
  stk::mesh::impl::fill_mesh_connectivity(*bulkPtr, meshConn);

  stk::mesh::PartVector parts = {&bulkPtr->mesh_meta_data().universal_part()};
  stk::mesh::create_all_block_boundary_sides(*bulkPtr, bulkPtr->mesh_meta_data().universal_part(), parts);

  auto [allocBlocksNeeded, modifiedEntities] = stk::mesh::impl::dealloc_and_get_needed_allocs_full_mesh(*bulkPtr, meshConn);
  EXPECT_TRUE(!allocBlocksNeeded.empty());
  std::vector<std::pair<size_t,size_t>> gold = { {32, 8}, {64, 26}, {128, 4} };
  EXPECT_EQ(gold, allocBlocksNeeded);
  EXPECT_EQ(38u, modifiedEntities.size());

  auto allocBlocksAvail = stk::mesh::impl::get_available_blocks(meshConn.get_memory_pool());
  gold = { {16, 336}, {32, 168}, {64, 88} };
  EXPECT_EQ(gold, allocBlocksAvail);
  EXPECT_FALSE(stk::mesh::impl::blocks_are_available(allocBlocksNeeded, allocBlocksAvail));
}

