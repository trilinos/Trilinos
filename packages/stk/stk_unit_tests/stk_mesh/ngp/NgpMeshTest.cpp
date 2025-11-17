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

#include "stk_ngp_test/ngp_test.hpp"
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpReductions.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/TestHexFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"

#include <limits>

using NgpMeshDefaultMemSpace = stk::mesh::NgpMeshDefaultMemSpace;

class NgpMeshTest : public stk::mesh::fixtures::TestHexFixture
{
public:
  void run_get_nodes_using_FastMeshIndex_test()
  {
    setup_mesh(1, 1, 4);

    stk::NgpVector<double> numNodesVec("numNodes", 1);

    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
                         KOKKOS_LAMBDA(const int /*i*/)
                         {
                           stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, stk::mesh::FastMeshIndex{0,0});
                           numNodesVec.device_get(0) = nodes.size();
                         });

    numNodesVec.copy_device_to_host();
    ASSERT_EQ(8u, numNodesVec[0]);
  }

  template <typename NgpMemSpace>
  void run_get_nodes_using_FastMeshIndex_test()
  {
    setup_mesh(1, 1, 4);

    stk::NgpVector<double> numNodesVec("numNodes", 1);

    stk::mesh::NgpMeshT<NgpMemSpace> & ngpMesh = stk::mesh::get_updated_ngp_mesh<NgpMemSpace>(get_bulk());
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
                         KOKKOS_LAMBDA(const int /*i*/)
                         {
                           stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, stk::mesh::FastMeshIndex{0,0});
                           numNodesVec.device_get(0) = nodes.size();
                         });

    numNodesVec.copy_device_to_host();
    ASSERT_EQ(8u, numNodesVec[0]);
  }

  void run_edge_check(unsigned numExpectedEdgesPerElem)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::ELEM_RANK, !stk::mesh::Selector());
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, bucketIds,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
        stk::mesh::ConnectedEntities edges = ngpMesh.get_edges(stk::topology::ELEM_RANK, entityIndex);
        NGP_EXPECT_EQ(numExpectedEdgesPerElem, edges.size());
      }, stk::ngp::ExecSpace()
    );
  }

  void delete_edge_on_each_element()
  {
    get_bulk().modification_begin();

    stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
    stk::mesh::ConnectedEntities edges = get_bulk().get_connected_entities(elem1, stk::topology::EDGE_RANK);
    stk::mesh::ConnectedEntities edgeElems = get_bulk().get_connected_entities(edges[0], stk::topology::ELEM_RANK);
    EXPECT_FALSE(edgeElems.empty());
    EXPECT_EQ(1u, edgeElems.size());
    EXPECT_EQ(elem1, edgeElems[0]);

    const stk::mesh::ConnectivityOrdinal* edgeElemOrds = get_bulk().begin_ordinals(edges[0], stk::topology::ELEM_RANK);
    stk::mesh::Entity edge = edges[0];
    EXPECT_TRUE(get_bulk().destroy_relation(elem1, edge, edgeElemOrds[0]));
    EXPECT_TRUE(get_bulk().destroy_entity(edge));

    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    edges = get_bulk().get_connected_entities(elem2, stk::topology::EDGE_RANK);
    EXPECT_EQ(12u, edges.size());
    edgeElems = get_bulk().get_connected_entities(edges[5], stk::topology::ELEM_RANK);
    EXPECT_EQ(1u, edgeElems.size());
    EXPECT_EQ(elem2, edgeElems[0]);
    edgeElemOrds = get_bulk().begin_ordinals(edges[5], stk::topology::ELEM_RANK);
    edge = edges[5];
    EXPECT_TRUE(get_bulk().destroy_relation(elem2, edge, edgeElemOrds[0]));
    EXPECT_TRUE(get_bulk().destroy_entity(edge));

    get_bulk().modification_end();
  }
};

NGP_TEST_F(NgpMeshTest, get_nodes_using_FastMeshIndex)
{
  run_get_nodes_using_FastMeshIndex_test();
}

NGP_TEST_F(NgpMeshTest, get_nodes_using_FastMeshIndex_custom_NgpMemSpace)
{
  run_get_nodes_using_FastMeshIndex_test<NgpMeshDefaultMemSpace>();
}

NGP_TEST_F(NgpMeshTest, hexes_with_edges_update_connectivity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  setup_mesh(1,1,2);
  stk::mesh::get_updated_ngp_mesh(get_bulk());

  stk::mesh::Part& edgePart = get_meta().declare_part("edges", stk::topology::EDGE_RANK);

  stk::mesh::create_edges(get_bulk(), get_meta().universal_part(), &edgePart);
  stk::mesh::get_updated_ngp_mesh(get_bulk());

  EXPECT_EQ(20u, stk::mesh::count_entities(get_bulk(), stk::topology::EDGE_RANK, edgePart));

  unsigned numExpectedEdgesPerElement = 12;
  run_edge_check(numExpectedEdgesPerElement);

  delete_edge_on_each_element();
  EXPECT_EQ(18u, stk::mesh::count_entities(get_bulk(), stk::topology::EDGE_RANK, edgePart));

  numExpectedEdgesPerElement = 11;
  run_edge_check(numExpectedEdgesPerElement);
}

class NgpMeshRankLimit : public stk::mesh::fixtures::TestHexFixture {};

TEST_F(NgpMeshRankLimit, tooManyRanksThrowWithMessage)
{
  setup_mesh(1,1,1, {"NODE","EDGE","FACE","ELEM","CONSTRAINT","JULIA"});

  try {
    stk::mesh::get_updated_ngp_mesh(get_bulk());
    FAIL()<< "expected throw but didn't throw";
  }
  catch(std::exception& e) {
    std::string expectedMsg("stk::mesh::NgpMesh: too many entity ranks (6). Required to be less-or-equal stk::topology::NUM_RANKS");
    std::string msg(e.what());
    EXPECT_TRUE((msg.find(expectedMsg) != std::string::npos));
  }
}

TEST_F(NgpMeshRankLimit, tooManyRanksThrowWithMessage_custom_NgpMemSpace)
{
  setup_mesh(1,1,1, {"NODE","EDGE","FACE","ELEM","CONSTRAINT","JULIA"});

  try {
    stk::mesh::get_updated_ngp_mesh<NgpMeshDefaultMemSpace>(get_bulk());
    FAIL()<< "expected throw but didn't throw";
  }
  catch(std::exception& e) {
    std::string expectedMsg("stk::mesh::NgpMesh: too many entity ranks (6). Required to be less-or-equal stk::topology::NUM_RANKS");
    std::string msg(e.what());
    EXPECT_TRUE((msg.find(expectedMsg) != std::string::npos));
  }
}

class NgpMeshHostDevice : public stk::mesh::fixtures::TestHexFixture {};

void test_HostMesh_works_on_host_in_any_build(const stk::mesh::BulkData& bulk)
{
  stk::mesh::HostMesh hostMesh(bulk);
  stk::mesh::Selector all = bulk.mesh_meta_data().universal_part();
  stk::NgpVector<unsigned> vec("elem-count", 1, 0);
  stk::mesh::for_each_entity_run(hostMesh, stk::topology::ELEM_RANK, all,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& sideIndex) {
      vec[0] += 1; //we're on host. on device we would use 'vec.device_get(0) += 1;'
    }
  );

  const unsigned numElems = vec[0];
  EXPECT_EQ(numElems, stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK, all));
}

TEST_F(NgpMeshHostDevice, host_mesh_host_space)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  setup_mesh(1,1,1);

  constexpr bool isCPUbuild = std::is_same_v<stk::ngp::MemSpace,stk::ngp::HostMemSpace>;

  if constexpr(isCPUbuild) {
    auto ngpMeshIsHostMesh = stk::mesh::get_updated_ngp_mesh<stk::ngp::HostMemSpace>(get_bulk());
#if defined(STK_USE_DEVICE_MESH) && defined(STK_ENABLE_GPU)
    ASSERT_TRUE((std::is_same_v<decltype(ngpMeshIsHostMesh),stk::mesh::HostMesh>));
#elif defined(STK_USE_DEVICE_MESH)
    ASSERT_TRUE((std::is_same_v<decltype(ngpMeshIsHostMesh),stk::mesh::DeviceMesh>));
#else
    ASSERT_TRUE((std::is_same_v<decltype(ngpMeshIsHostMesh),stk::mesh::HostMesh>));
#endif
  }
  else {
    auto ngpMeshIsDeviceMesh = stk::mesh::get_updated_ngp_mesh<stk::ngp::MemSpace>(get_bulk());
    ASSERT_TRUE((std::is_same_v<decltype(ngpMeshIsDeviceMesh),stk::mesh::DeviceMesh>));
  }

  test_HostMesh_works_on_host_in_any_build(get_bulk());
}

class EntityIndexSpace : public stk::mesh::fixtures::TestHexFixture {};

TEST_F(EntityIndexSpace, accessingLocalData_useLocalOffset)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  setup_mesh(1, 1, 1);
  std::vector<unsigned> entityToLocalOffset(get_bulk().get_size_of_entity_index_space(), 0);

  stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(get_meta().entity_rank_count());
  for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; ++rank)
  {
    unsigned localOffset = 0;
    const stk::mesh::BucketVector &buckets = get_bulk().buckets(stk::topology::NODE_RANK);
    for(const stk::mesh::Bucket *bucket : buckets)
    {
      for(stk::mesh::Entity entity : *bucket)
      {
        entityToLocalOffset[entity.local_offset()] = localOffset;
        localOffset++;
      }
    }
  }

  std::vector<unsigned> gold {0,0,0,1,3,2,4,5,7,6};
  ASSERT_EQ(gold.size(), entityToLocalOffset.size());
  for(size_t i=0; i<gold.size(); i++)
  {
    EXPECT_EQ(gold[i], entityToLocalOffset[i]);
  }
}

void run_vector_gpu_test()
{
  size_t n = 10;
  stk::NgpVector<double> vec("vec", n);
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, n),
                       KOKKOS_LAMBDA(const int i)
                       {
                         vec.device_get(i) = i;
                       });
  vec.copy_device_to_host();
  for(size_t i=0; i<n; i++)
    EXPECT_EQ(i, vec[i]);
}

TEST(StkVectorGpuTest, gpu_runs)
{
  run_vector_gpu_test();
}

void check_volatile_fast_shared_comm_map_values_on_device(
    const stk::mesh::NgpMesh & ngpMesh, int proc,
    const stk::mesh::DeviceCommMapIndices<stk::ngp::MemSpace> & deviceCommMapIndicesGold)
{
  auto test = KOKKOS_LAMBDA(size_t /*i*/)
              {
                stk::mesh::DeviceCommMapIndices<stk::ngp::MemSpace> deviceCommMapIndices =
                    ngpMesh.volatile_fast_shared_comm_map(stk::topology::NODE_RANK, proc);

                for (size_t entry = 0; entry < deviceCommMapIndices.size(); ++entry) {
                  NGP_EXPECT_EQ(deviceCommMapIndicesGold[entry].bucket_id, deviceCommMapIndices[entry].bucket_id);
                  NGP_EXPECT_EQ(deviceCommMapIndicesGold[entry].bucket_ord, deviceCommMapIndices[entry].bucket_ord);
                }
              };

  if constexpr (std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>) {
    test(0);
  } else {
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), test);
  }
}

NGP_TEST_F(NgpMeshTest, volatileFastSharedCommMap)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) { GTEST_SKIP(); }

  setup_mesh(1, 1, 4);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  std::vector<int> comm_procs = get_bulk().all_sharing_procs(stk::topology::NODE_RANK);

  for (int proc : comm_procs) {
    stk::mesh::HostCommMapIndices<stk::ngp::MemSpace> hostNgpMeshIndices =
        get_bulk().template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(stk::topology::NODE_RANK, proc);
    stk::mesh::DeviceCommMapIndices<stk::ngp::MemSpace> deviceNgpMeshIndices("deviceNgpMeshIndices",
                                                                             hostNgpMeshIndices.extent(0));

    Kokkos::deep_copy(deviceNgpMeshIndices, hostNgpMeshIndices);
    check_volatile_fast_shared_comm_map_values_on_device(ngpMesh, proc, deviceNgpMeshIndices);
  }
}

NGP_TEST_F(NgpMeshTest, volatileFastSharedCommMap_custom_NgpMemSpace)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) { GTEST_SKIP(); }

  setup_mesh(1, 1, 4);

  stk::mesh::NgpMeshT<stk::mesh::NgpMeshDefaultMemSpace> & ngpMesh = stk::mesh::get_updated_ngp_mesh<stk::mesh::NgpMeshDefaultMemSpace>(get_bulk());
  std::vector<int> comm_procs = get_bulk().all_sharing_procs(stk::topology::NODE_RANK);

  for (int proc : comm_procs) {
    stk::mesh::HostCommMapIndices<stk::ngp::MemSpace> hostNgpMeshIndices =
        get_bulk().template volatile_fast_shared_comm_map<stk::mesh::NgpMeshDefaultMemSpace>(stk::topology::NODE_RANK, proc);
    stk::mesh::DeviceCommMapIndices<stk::ngp::MemSpace> deviceNgpMeshIndices("deviceNgpMeshIndices",
                                                                             hostNgpMeshIndices.extent(0));

    Kokkos::deep_copy(deviceNgpMeshIndices, hostNgpMeshIndices);
    check_volatile_fast_shared_comm_map_values_on_device(ngpMesh, proc, deviceNgpMeshIndices);
  }
}

void test_ngp_permutations_1side_2perms(const stk::mesh::BulkData& mesh,
                                        const stk::mesh::Part& sidePart)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);

  stk::mesh::EntityRank sideRank = mesh.mesh_meta_data().side_rank();
  stk::mesh::EntityVector sides;
  stk::mesh::get_entities(mesh, sideRank, sidePart, sides);
  EXPECT_EQ(1u, sides.size());
  EXPECT_EQ(2u, mesh.num_connectivity(sides[0], stk::topology::ELEM_RANK));
  const stk::mesh::Permutation* hostPerms = mesh.begin_permutations(sides[0], stk::topology::ELEM_RANK);
  stk::mesh::Permutation expectedPerm1 = hostPerms[0];
  stk::mesh::Permutation expectedPerm2 = hostPerms[1];

  stk::mesh::for_each_entity_run(ngpMesh, sideRank, sidePart,
  KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& sideIndex) {
    stk::mesh::NgpMesh::Permutations perms = ngpMesh.get_permutations(sideRank, sideIndex, stk::topology::ELEM_RANK);
    NGP_EXPECT_EQ(2u, perms.size());
    NGP_EXPECT_EQ(expectedPerm1, perms[0]);
    NGP_EXPECT_EQ(expectedPerm2, perms[1]);
  });
}

NGP_TEST(TestNgpMesh, permutations)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::string meshDesc =
    "0,1,TRI_3_2D,1,2,3,block_1\n"
    "0,2,TRI_3_2D,2,4,3,block_2\n"
    "|dimension:2|sideset:name=surface_1; data=1,2";

  std::shared_ptr<stk::mesh::BulkData> mesh = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                  .set_spatial_dimension(2).create();
  stk::unit_test_util::setup_text_mesh(*mesh, meshDesc);

  stk::mesh::EntityRank sideRank = mesh->mesh_meta_data().side_rank();
  stk::mesh::Part* sidePart = mesh->mesh_meta_data().get_part("surface_1");
  STK_ThrowAssertMsg(sidePart != nullptr, "failed to find part for surface_1");

  stk::mesh::EntityVector sides;
  stk::mesh::get_entities(*mesh, sideRank, *sidePart, sides);
  EXPECT_EQ(1u, sides.size());
  EXPECT_EQ(2u, mesh->num_connectivity(sides[0], stk::topology::ELEM_RANK));

  stk::mesh::Permutation expectedPerm1 = static_cast<stk::mesh::Permutation>(0);
  stk::mesh::Permutation expectedPerm2 = static_cast<stk::mesh::Permutation>(1);
  const stk::mesh::Permutation* permutations = mesh->begin_permutations(sides[0], stk::topology::ELEM_RANK);
  EXPECT_EQ(expectedPerm1, permutations[0]);
  EXPECT_EQ(expectedPerm2, permutations[1]);

  test_ngp_permutations_1side_2perms(*mesh, *sidePart);
}

void test_ngp_local_ids(const stk::mesh::BulkData& mesh)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);

  stk::mesh::Entity elem1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
  unsigned expectedElem1LocalId = mesh.local_id(elem1);
  unsigned expectedNode1LocalId = mesh.local_id(node1);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(const int /*i*/) {
    unsigned elem1LocalId = ngpMesh.local_id(elem1);
    unsigned node1LocalId = ngpMesh.local_id(node1);
    NGP_EXPECT_EQ(expectedElem1LocalId, elem1LocalId);
    NGP_EXPECT_EQ(expectedNode1LocalId, node1LocalId);
  });
}

NGP_TEST(TestNgpMesh, local_ids)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::string meshDesc =
    "0,1,TRI_3_2D,1,2,3,block_1\n"
    "0,2,TRI_3_2D,2,4,3,block_2\n"
    "|dimension:2|sideset:name=surface_1; data=1,2";

  std::shared_ptr<stk::mesh::BulkData> mesh = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                  .set_spatial_dimension(2).create();
  stk::unit_test_util::setup_text_mesh(*mesh, meshDesc);

  test_ngp_local_ids(*mesh);
}

namespace {
double reduce_on_host(stk::mesh::BulkData& bulk)
{
  auto ngp_mesh = stk::mesh::HostMesh(bulk);

  double max_val = 0.0;
  Kokkos::Max<double> max_reduction(max_val);
  stk::mesh::for_each_entity_reduce(
    ngp_mesh,
    stk::topology::NODE_RANK,
    !stk::mesh::Selector(),
    max_reduction,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex /*i*/, double& thread_value) {
      max_reduction.join(thread_value, 1.0);
    });

  return max_val;
}
}

TEST(NgpHostMesh, FieldForEachEntityReduceOnHost_fromTylerVoskuilen)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }
  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, 2, 2, 2);
  fixture.m_meta.commit();
  fixture.generate_mesh();

  auto maxZ = reduce_on_host(fixture.m_bulk_data);
  EXPECT_EQ(1.0, maxZ);
}

TEST(NgpDeviceMesh, dont_let_stacksize_get_out_of_control)
{
  constexpr size_t tol = 64;

#ifdef SIERRA_MIGRATION
  constexpr size_t expectedBulkDataSize = 1320;
#else
  constexpr size_t expectedBulkDataSize = 1256;
#endif
  EXPECT_NEAR(expectedBulkDataSize, sizeof(stk::mesh::BulkData), tol);

  constexpr size_t expectedBucketSize = 976;
  EXPECT_NEAR(expectedBucketSize, sizeof(stk::mesh::Bucket), tol);

  constexpr size_t expectedDeviceMeshSize = 920;
  EXPECT_NEAR(expectedDeviceMeshSize, sizeof(stk::mesh::DeviceMesh), tol);

  constexpr size_t expectedDeviceBucketSize = 264;
  EXPECT_NEAR(expectedDeviceBucketSize, sizeof(stk::mesh::DeviceBucket), tol);
}

void add_elements(std::unique_ptr<stk::mesh::BulkData>& bulk)
{
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  stk::mesh::Part& part_1 = meta.declare_part_with_topology("part_1", stk::topology::HEX_8);

  const int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const stk::mesh::EntityId elemId = rank + 1;
  const stk::mesh::EntityId firstNodeId = rank * 8 + 1;

  stk::mesh::EntityIdVector nodeIds(8, 0);
  for (unsigned i = 0; i < nodeIds.size(); ++i) {
    nodeIds[i] = firstNodeId + i;
  }

  bulk->modification_begin();
  stk::mesh::declare_element(*bulk, part_1, elemId, nodeIds);
  bulk->modification_end();
}

TEST(NgpTeardownOrdering, BulkDataOutlivesNgpMesh)
{
  std::unique_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
  add_elements(bulk);

  [[maybe_unused]] stk::mesh::NgpMesh ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulk);

  // The "expect" for this test is a clean Valgrind run and no seg-faults
}

TEST(NgpTeardownOrdering, NgpMeshOutlivesBulkData)
{
  stk::mesh::NgpMesh ngpMesh;
  std::unique_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
  add_elements(bulk);

  ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulk);

  // The "expect" for this test is a clean Valgrind run and no seg-faults
}

