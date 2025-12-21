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
#include <stddef.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_topology/topology.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stkMeshTestUtils.hpp>

namespace stk { namespace mesh { class Part; } }

namespace
{
using DeviceEntityViewType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;
using HostEntityViewType = typename DeviceEntityViewType::host_mirror_type;

using DevicePartOrdinalViewType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::MemSpace>;
using HostPartOrdinalViewType = typename DevicePartOrdinalViewType::host_mirror_type;

void check_entity_parts(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& entities, const stk::mesh::PartVector& addedParts, const stk::mesh::PartVector& removedParts)
{
  for(stk::mesh::Entity entity : entities) {
    stk::mesh::Bucket& bucket = bulk.bucket(entity);

    for (auto part : addedParts) {
      EXPECT_TRUE(bucket.member(*part));
    }
    for (auto part : removedParts) {
      EXPECT_FALSE(bucket.member(*part));
    }
  }
}

void check_entity_parts(const stk::mesh::NgpMesh& ngpMesh, const DeviceEntityViewType& deviceEntityView, const DevicePartOrdinalViewType& addedPartOrdinals, const DevicePartOrdinalViewType& removedPartOrdinals)
{
  Kokkos::parallel_for(deviceEntityView.extent(0),
    KOKKOS_LAMBDA(const int entityIdx) {
      auto entity = deviceEntityView(entityIdx);
      auto rank = ngpMesh.entity_rank(entity);
      auto fastMeshIndex = ngpMesh.fast_mesh_index(entity);
      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIndex.bucket_id);

      for (unsigned i = 0; i < addedPartOrdinals.extent(0); ++i) {
        auto partOrdinal = addedPartOrdinals(i);
        NGP_EXPECT_TRUE(bucket.member(partOrdinal));
      }

      for (unsigned i = 0; i < removedPartOrdinals.extent(0); ++i) {
        auto partOrdinal = removedPartOrdinals(i);
        NGP_EXPECT_FALSE(bucket.member(partOrdinal));
      }
    }
  );
}

void populate_device_entity_view(DeviceEntityViewType& deviceEntityView, stk::mesh::EntityVector& entities)
{
  if (deviceEntityView.extent(0) != entities.size()) {
    Kokkos::resize(deviceEntityView, entities.size());
  }

  auto hostEntityView = Kokkos::create_mirror_view(deviceEntityView);

  for (unsigned i = 0; i < hostEntityView.extent(0); ++i) {
    hostEntityView(i) = entities[i];
  }
  Kokkos::deep_copy(deviceEntityView, hostEntityView);
}

//BEGIN
NGP_TEST(StkMeshHowTo, NgpMeshBatchChangeEntityParts)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(communicator) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  unsigned elementCount = 10;
  stk::io::fill_mesh("generated:1x1x" + std::to_string(elementCount), *bulkPtr);

  stk::mesh::Part* block1Part = meta.get_part("block_1");
  stk::mesh::Part* block2Part = &meta.declare_part("block_2", stk::topology::ELEM_RANK);
  stk::mesh::Selector block1Selector(*block1Part);
  stk::mesh::EntityVector entities;

  stk::mesh::get_entities(*bulkPtr, stk::topology::ELEM_RANK,  block1Selector, entities);
  EXPECT_EQ(elementCount, entities.size());
  check_entity_parts(*bulkPtr, entities, stk::mesh::PartVector{block1Part}, stk::mesh::PartVector{});

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);

  DeviceEntityViewType deviceEntityView("entityView", entities.size());
  populate_device_entity_view(deviceEntityView, entities);

  HostPartOrdinalViewType addParts("", 1);
  addParts(0) = block2Part->mesh_meta_data_ordinal();
  auto deviceAddParts = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, addParts);

  HostPartOrdinalViewType removeParts("", 1);
  removeParts(0) = block1Part->mesh_meta_data_ordinal();
  auto deviceRemoveParts = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, removeParts);

  ngpMesh.batch_change_entity_parts(deviceEntityView, deviceAddParts, deviceRemoveParts);

  check_entity_parts(ngpMesh, deviceEntityView, deviceAddParts, deviceRemoveParts);

  ngpMesh.update_bulk_data();

  check_entity_parts(*bulkPtr, entities, stk::mesh::PartVector{block2Part}, stk::mesh::PartVector{block1Part});
}
//END

}

