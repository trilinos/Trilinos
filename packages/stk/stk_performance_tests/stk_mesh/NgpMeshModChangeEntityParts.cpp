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
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/timer.hpp>
#include <algorithm>
#include <string>
#include <vector>

class NgpMeshModChangeEntityPartsPerfTester : public stk::unit_test_util::MeshFixture
{
public:
  using DeviceEntitiesViewType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;
  using DevicePartOrdinalsViewType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::MemSpace>;

  NgpMeshModChangeEntityPartsPerfTester()
    : stk::unit_test_util::MeshFixture(),
      numElemsPerDim(0),
      numIters(0),
      numRuns(0),
      batchTimer(get_comm())
  {}

  void setup_test_params(unsigned defaultNumElemsPerDim = 20, unsigned defaultNumIters = 20, unsigned defaultNumRuns = 5)
  {
    numElemsPerDim = stk::unit_test_util::get_command_line_option("-e", defaultNumElemsPerDim);
    numIters = stk::unit_test_util::get_command_line_option("-i", defaultNumIters);
    numRuns = stk::unit_test_util::get_command_line_option("-r", defaultNumRuns);
    batchTimer.initialize_batch_timer();
  }

  void setup_host_mesh(unsigned initialBucketCapacity = stk::mesh::get_default_initial_bucket_capacity(),
                       unsigned maximumBucketCapacity = stk::mesh::get_default_maximum_bucket_capacity())
  {
    const std::string elemsPerDim = std::to_string(numElemsPerDim);
    setup_mesh("generated:" + elemsPerDim + "x" + elemsPerDim + "x" + elemsPerDim,
               stk::mesh::BulkData::NO_AUTO_AURA, initialBucketCapacity, maximumBucketCapacity);
  }

  stk::mesh::Part& declare_part(std::string const& name, stk::mesh::EntityRank rank = stk::topology::INVALID_RANK)
  {
    if (rank == stk::topology::INVALID_RANK) {
      return get_meta().declare_part(name);
    }
    return get_meta().declare_part(name, rank);
  }

  stk::mesh::EntityVector get_entities(stk::mesh::EntityRank rank, stk::mesh::Selector const& selector)
  {
    stk::mesh::EntityVector entities;
    stk::mesh::get_entities(get_bulk(), rank, selector, entities);
    return entities;
  }

  stk::mesh::EntityVector get_every_nth_entity(stk::mesh::EntityVector const& entities, unsigned stride)
  {
    stk::mesh::EntityVector strided;
    for (unsigned i = 0; i < entities.size(); i += stride) {
      strided.push_back(entities[i]);
    }
    return strided;
  }

  DeviceEntitiesViewType get_device_entities_view(stk::mesh::EntityVector const& entities)
  {
    DeviceEntitiesViewType deviceEntities("deviceEntities", entities.size());
    auto hostEntities = Kokkos::create_mirror_view(deviceEntities);
    for (unsigned i = 0; i < entities.size(); ++i) {
      hostEntities(i) = entities[i];
    }
    Kokkos::deep_copy(deviceEntities, hostEntities);
    return deviceEntities;
  }

  DevicePartOrdinalsViewType get_device_part_ordinals_view(stk::mesh::PartVector const& parts)
  {
    DevicePartOrdinalsViewType devicePartOrdinals("devicePartOrdinals", parts.size());
    auto hostPartOrdinals = Kokkos::create_mirror_view(devicePartOrdinals);
    for (unsigned i = 0; i < parts.size(); ++i) {
      hostPartOrdinals(i) = parts[i]->mesh_meta_data_ordinal();
    }
    Kokkos::deep_copy(devicePartOrdinals, hostPartOrdinals);
    return devicePartOrdinals;
  }

  void time_device_change_parts(stk::mesh::EntityVector const& entities,
                                stk::mesh::PartVector const& addParts,
                                stk::mesh::PartVector const& removeParts,
                                unsigned iters = 1)
  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    DeviceEntitiesViewType deviceEntities = get_device_entities_view(entities);
    DevicePartOrdinalsViewType deviceAddParts = get_device_part_ordinals_view(addParts);
    DevicePartOrdinalsViewType deviceRemoveParts = get_device_part_ordinals_view(removeParts);

    batchTimer.start_batch_timer();
    Kokkos::Profiling::pushRegion("batch_device_change_entity_parts_perf_test");

    for (unsigned iter = 0; iter < iters; ++iter) {
      if (iter % 2 == 0) {
        ngpMesh.batch_change_entity_parts(deviceEntities, deviceAddParts, deviceRemoveParts);
      }
      else {
        ngpMesh.batch_change_entity_parts(deviceEntities, deviceRemoveParts, deviceAddParts);
      }
    }
    Kokkos::fence();

    Kokkos::Profiling::popRegion();
    batchTimer.stop_batch_timer();
  }

  unsigned numElemsPerDim;
  unsigned numIters;
  unsigned numRuns;
  stk::unit_test_util::BatchTimer batchTimer;
};

TEST_F(NgpMeshModChangeEntityPartsPerfTester, AddRankedNodeParts)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  setup_test_params(80, 1, 3);

  for (unsigned run = 0; run < numRuns; ++run) {
    setup_host_mesh();
    auto& nodePart = declare_part("node_part", stk::topology::NODE_RANK);
    auto elements = get_entities(stk::topology::NODE_RANK, get_meta().universal_part());

    time_device_change_parts(elements, {&nodePart}, {});
    reset_mesh();
  }
  batchTimer.print_batch_timing(1);
}

TEST_F(NgpMeshModChangeEntityPartsPerfTester, MoveUnrankedElementParts)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  setup_test_params(80, 1, 3);

  unsigned bucketCapacity = stk::mesh::get_default_maximum_bucket_capacity();
  for (unsigned run = 0; run < numRuns; ++run) {
    setup_host_mesh(bucketCapacity, bucketCapacity);
    auto& unrankedPart = declare_part("unranked_part");

    time_device_change_parts(get_entities(stk::topology::ELEM_RANK, get_meta().universal_part()), {&unrankedPart}, {});
    reset_mesh();
  }
  batchTimer.print_batch_timing(1);
}

TEST_F(NgpMeshModChangeEntityPartsPerfTester, PartialMoveRankedElementParts)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  setup_test_params(70, 1, 3);

  auto fraction = stk::unit_test_util::get_command_line_option("-f", 10u);
  auto stride = std::max(1u, 100u / std::max(1u, fraction));

  for (unsigned run = 0; run < numRuns; ++run) {
    setup_host_mesh();
    auto& block1 = *get_meta().get_part("block_1");
    auto& block2 = declare_part("block_2", stk::topology::ELEM_RANK);
    auto elements = get_entities(stk::topology::ELEM_RANK, block1);

    time_device_change_parts(get_every_nth_entity(elements, stride), {&block2}, {&block1});
    reset_mesh();
  }
  batchTimer.print_batch_timing(1);
}

TEST_F(NgpMeshModChangeEntityPartsPerfTester, MoveRankedElementPartsInSequence)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  setup_test_params(30, 6, 3);

  for (unsigned run = 0; run < numRuns; ++run) {
    setup_host_mesh();
    auto& block1 = *get_meta().get_part("block_1");
    auto& block2 = declare_part("block_2", stk::topology::ELEM_RANK);
    auto elements = get_entities(stk::topology::ELEM_RANK, block1);

    time_device_change_parts(elements, {&block2}, {&block1}, numIters);
    reset_mesh();
  }
  batchTimer.print_batch_timing(numIters);
}
