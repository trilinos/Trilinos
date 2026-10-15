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
#include <Kokkos_Core.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace {

// Reads the ghost-inclusive communication-map size on the device.  It is queried inside a device
// kernel so it exercises the same volatileFastSharedCommMap views that a real GPU parallel data
// exchange would index.
std::size_t get_device_ghost_comm_map_size(const stk::mesh::NgpMesh& ngpMesh, int proc)
{
  Kokkos::View<std::size_t, stk::ngp::MemSpace> deviceSize("device_ghost_comm_map_size");
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
                       KOKKOS_LAMBDA(const int /*i*/) {
                         deviceSize() =
                             ngpMesh.volatile_fast_shared_comm_map(stk::topology::NODE_RANK, proc, true).extent(0);
                       });

  auto hostSize = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), deviceSize);
  return hostSize();
}

std::size_t get_host_ghost_comm_map_size(const stk::mesh::BulkData& bulk, int proc)
{
  return bulk.template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(stk::topology::NODE_RANK, proc, true)
      .extent(0);
}

// A sender-side custom ghosting of an already-owned node changes only the sender's communication
// metadata; it does not add, remove, or move any of its local buckets.  get_updated_ngp_mesh() must
// still refresh the device communication map, otherwise NgpParallelDataExchange sizes its work from
// the (updated) host map but indexes a stale device map, silently dropping sender contributions.
TEST(NgpGhostCommMapUpdate, SenderOnlyGhostingRefreshesDeviceCommMap)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP() << "requires exactly two MPI ranks"; }

  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int otherRank = 1 - parallelRank;

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                  .set_spatial_dimension(3)
                                                  .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                                  .set_symmetric_ghost_info(true)
                                                  .create();
  stk::io::fill_mesh("generated:1x1x4", *bulk);

  // Prime the device mesh before ghosting, as an application does when a GPU kernel uses the mesh
  // before a later contact-search modification.
  const stk::mesh::NgpMesh& ngpMeshBefore = stk::mesh::get_updated_ngp_mesh(*bulk);
  const std::size_t hostSizeBefore = get_host_ghost_comm_map_size(*bulk, otherRank);
  EXPECT_EQ(get_device_ghost_comm_map_size(ngpMeshBefore, otherRank), hostSizeBefore);

  bulk->modification_begin();
  stk::mesh::Ghosting& customGhosting = bulk->create_ghosting("sender_only_ghosting");
  std::vector<stk::mesh::EntityProc> entitiesToGhost;
  if (parallelRank == 0) {
    const stk::mesh::Entity node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    ASSERT_TRUE(bulk->is_valid(node1));
    ASSERT_TRUE(bulk->bucket(node1).owned());
    // Rank 0 already owns node 1.  Sending it to rank 1 changes rank 0's communication metadata but
    // does not add, remove, or move a local bucket.
    entitiesToGhost.emplace_back(node1, 1);
  }
  bulk->change_ghosting(customGhosting, entitiesToGhost);
  bulk->modification_end();

  const stk::mesh::NgpMesh& ngpMeshAfter = stk::mesh::get_updated_ngp_mesh(*bulk);
  const std::size_t hostSizeAfter = get_host_ghost_comm_map_size(*bulk, otherRank);
  const std::size_t deviceSizeAfter = get_device_ghost_comm_map_size(ngpMeshAfter, otherRank);

  // Both ranks now have one additional ghost communication entry.  An updated device mesh must
  // expose the same entries as BulkData even on the sender, where bucket membership did not change.
  EXPECT_EQ(hostSizeAfter, hostSizeBefore + 1);
  EXPECT_EQ(deviceSizeAfter, hostSizeAfter);
}

}  // namespace
