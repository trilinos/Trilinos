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
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"

#include <limits>

class NgpMeshTest : public stk::mesh::fixtures::TestHexFixture
{
public:
  void run_get_nodes_using_FastMeshIndex_test()
  {
    setup_mesh(1, 1, 4);

    stk::NgpVector<double> numNodesVec("numNodes", 1);

    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
                         KOKKOS_LAMBDA(const int i)
                         {
                           stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, stk::mesh::FastMeshIndex{0,0});
                           numNodesVec.device_get(0) = nodes.size();
                         });

    numNodesVec.copy_device_to_host();
    ASSERT_EQ(8u, numNodesVec[0]);
  }
};
TEST_F(NgpMeshTest, get_nodes_using_FastMeshIndex)
{
  run_get_nodes_using_FastMeshIndex_test();
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

void check_volatile_fast_shared_comm_map_values_on_device(const stk::mesh::NgpMesh & ngpMesh, int proc, const stk::mesh::DeviceCommMapIndices & deviceCommMapIndicesGold)
{
  auto test = KOKKOS_LAMBDA(size_t i)
              {
                stk::mesh::DeviceCommMapIndices deviceCommMapIndices = ngpMesh.volatile_fast_shared_comm_map(stk::topology::NODE_RANK, proc);

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

using HostCommMapIndices = Kokkos::View<stk::mesh::FastMeshIndex*, stk::ngp::HostExecSpace>;

NGP_TEST_F(NgpMeshTest, volatileFastSharedCommMap)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) { GTEST_SKIP(); }

  setup_mesh(1, 1, 4);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  std::vector<int> comm_procs = get_bulk().all_sharing_procs(stk::topology::NODE_RANK);

  for (int proc : comm_procs) {
    stk::mesh::DeviceCommMapIndices::HostMirror hostNgpMeshIndices = get_bulk().volatile_fast_shared_comm_map(stk::topology::NODE_RANK, proc);
    stk::mesh::DeviceCommMapIndices deviceNgpMeshIndices("deviceNgpMeshIndices", hostNgpMeshIndices.extent(0));

    Kokkos::deep_copy(deviceNgpMeshIndices, hostNgpMeshIndices);
    check_volatile_fast_shared_comm_map_values_on_device(ngpMesh, proc, deviceNgpMeshIndices);
  }
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
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex i, double& thread_value) {
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

