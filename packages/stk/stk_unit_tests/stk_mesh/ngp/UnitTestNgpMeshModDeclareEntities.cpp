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
#include <algorithm>
#include <iterator>
#include <map>
#include "Kokkos_Core.hpp"
#include "ngp/NgpUnitTestUtils.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/SkinMesh.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/SkinBoundary.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_util/command_line/CommandLineParser.hpp"
#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_mesh/base/SkinBoundary.hpp"

#include "ngp/UnitTestNgpMeshModificationUtils.hpp"

#ifdef STK_USE_DEVICE_MESH
namespace
{
class NgpBatchDeclareEntities : public NgpBatchDeclareDestroyEntities
{
public:
  NgpBatchDeclareEntities() {};

  std::vector<stk::mesh::PartOrdinal> extract_part_ordinals(const stk::mesh::PartVector& addParts) const {
    std::vector<stk::mesh::PartOrdinal> addPartOrdinals;
    std::ranges::transform(addParts, std::back_inserter(addPartOrdinals), [](const stk::mesh::Part* part) {
      return part->mesh_meta_data_ordinal();
    });
    return addPartOrdinals;
  }
  void declare_on_host(const std::vector<unsigned>& entityIdsToDeclare, const stk::mesh::PartVector& addParts)
  {
    Kokkos::View<unsigned*, Kokkos::HostSpace> entityIds("entityIds", entityIdsToDeclare.size());
    fill_views(entityIds, entityIdsToDeclare);

    HostPartOrdinalsType addPartOrdinals("addPartOrdinals", addParts.size());
    fill_views(addPartOrdinals, extract_part_ordinals(addParts));

    stk::mesh::HostMesh hostMesh(*m_bulk);

    HostEntitiesType createdEntities("createdEntities", entityIdsToDeclare.size());
    hostMesh.batch_declare_entities(stk::topology::NODE_RANK, entityIds, addPartOrdinals, createdEntities);
    confirm_host_mesh_is_not_synchronized_from_device(hostMesh);

    hostMesh.update_bulk_data();
    confirm_host_mesh_is_synchronized_from_device(hostMesh);
  }

  void declare_on_device(const std::vector<unsigned>& entityIdsToDeclare, const stk::mesh::PartVector& addParts)
  {
    Kokkos::View<unsigned*> entityIds("entityIds", entityIdsToDeclare.size());
    fill_views(entityIds, entityIdsToDeclare);

    DevicePartOrdinalsType addPartOrdinals("addPartOrdinals", addParts.size());
    fill_views(addPartOrdinals, extract_part_ordinals(addParts));

    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);

    DeviceEntitiesType createdEntities("createdEntities", entityIdsToDeclare.size());
    ngpMesh.batch_declare_entities(stk::topology::NODE_RANK, entityIds, addPartOrdinals, createdEntities);
    confirm_host_mesh_is_not_synchronized_from_device(ngpMesh);

    ngpMesh.update_bulk_data();
    confirm_host_mesh_is_synchronized_from_device(ngpMesh);
  }
};

class NgpBatchDeclareNodes : public NgpBatchDeclareEntities
{
public:
  NgpBatchDeclareNodes() {};

  void declare_on_host(const std::vector<unsigned>& entityIdsToDeclare, const stk::mesh::PartVector& addParts)
  {
    NgpBatchDeclareEntities::declare_on_host(entityIdsToDeclare, addParts);
  }

  void declare_on_device(const std::vector<unsigned>& entityIdsToDeclare, const stk::mesh::PartVector& addParts)
  {
    NgpBatchDeclareEntities::declare_on_device(entityIdsToDeclare, addParts);
  }

  static constexpr unsigned nodeId1 = 1;
  static constexpr unsigned nodeId2 = 2;
};

NGP_TEST_F(NgpBatchDeclareNodes, CreateOneNode_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_host({nodeId1}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDeclareNodes, CreateOneNode_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_device({nodeId1}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDeclareNodes, CreateTwoNodes_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_host({nodeId1, nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}, {{"part1"}, {nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateTwoNodes_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_device({nodeId1, nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}, {{"part1"}, {nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateTwoNodes_OneBucket_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 2);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_host({nodeId1, nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1, nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateTwoNodes_OneBucket_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 2);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_device({nodeId1, nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1, nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateOneNode_OneExisting_DifferentBucket_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  [[maybe_unused]] const stk::mesh::Entity entity1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  declare_on_host({nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}, {{"part1"}, {nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateOneNode_OneExisting_DifferentBucket_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  [[maybe_unused]] const stk::mesh::Entity entity1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  declare_on_device({nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}, {{"part1"}, {nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateOneNode_OneExisting_SameBucket_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 2);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  [[maybe_unused]] const stk::mesh::Entity entity1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  declare_on_host({nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1, nodeId2}}}, stk::topology::NODE_RANK);
}

TEST_F(NgpBatchDeclareNodes, CreateOneNode_OneExisting_SameBucket_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(2, 2);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  [[maybe_unused]] const stk::mesh::Entity entity1 = create_node(*m_bulk, nodeId1, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1}}}, stk::topology::NODE_RANK);

  declare_on_device({nodeId2}, {&part1});

  check_bucket_layout(*m_bulk, {{{"part1"}, {nodeId1, nodeId2}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDeclareNodes, CreateOneNode_twoParts_ngpHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_host({nodeId1}, {&part1, &part2});

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

NGP_TEST_F(NgpBatchDeclareNodes, CreateOneNode_twoParts_ngpDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part & part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  stk::mesh::Part & part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);

  check_bucket_layout(*m_bulk, {}, stk::topology::NODE_RANK);

  declare_on_device({nodeId1}, {&part1, &part2});

  check_bucket_layout(*m_bulk, {{{"part1", "part2"}, {nodeId1}}}, stk::topology::NODE_RANK);
}

}  // namespace
#endif
