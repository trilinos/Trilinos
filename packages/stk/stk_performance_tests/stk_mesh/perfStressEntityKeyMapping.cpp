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

#include <string>
#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <cstdlib>

using EntityIdPair = std::pair<stk::mesh::EntityId,stk::mesh::EntityId>;

class StressEntityKeyMapping : public stk::unit_test_util::MeshFixture
{
public:
  StressEntityKeyMapping()
  { }

  void setup_host_mesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    stk::io::fill_mesh("generated:500x500x10", get_bulk());
  }

  EntityIdPair get_min_max_elem_ids_on_local_proc()
  {
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stk::topology::ELEM_RANK, get_meta().locally_owned_part());
    stk::mesh::EntityId minElemId = std::numeric_limits<stk::mesh::EntityId>::max();
    stk::mesh::EntityId maxElemId = 0;

    for(const stk::mesh::Bucket* bucket : buckets) {
      minElemId = std::min(minElemId, get_bulk().identifier((*bucket)[0]));
      stk::mesh::Entity lastElem = *std::prev(bucket->end());
      maxElemId = std::max(maxElemId, get_bulk().identifier(lastElem));
    }
    return std::make_pair(minElemId, maxElemId);
  }

  stk::mesh::Entity get_owned_element(stk::mesh::Entity node)
  {
    const unsigned numElems = get_bulk().num_elements(node);
    const stk::mesh::Entity* elems = get_bulk().begin_elements(node);
    for(unsigned i=0; i<numElems; ++i) {
      if (get_bulk().bucket(elems[i]).owned()) {
        return elems[i];
      }
    }
    return stk::mesh::Entity();
  }

  void add_1_element_to_block_2()
  {
    stk::mesh::Part& block2 = *get_meta().get_part("block_2");
    stk::mesh::Selector sharedButNotBlock2 = get_meta().globally_shared_part() & !block2;
    const stk::mesh::BucketVector& nodeBuckets = get_bulk().get_buckets(stk::topology::NODE_RANK, sharedButNotBlock2);
    STK_ThrowRequire(!nodeBuckets.empty());
    STK_ThrowRequire(nodeBuckets[0]->size() > 0);

    stk::mesh::Entity node = (*nodeBuckets[0])[0];
    stk::mesh::Entity elem = get_owned_element(node);
    STK_ThrowRequire(get_bulk().is_valid(elem));
    get_bulk().change_entity_parts(elem, stk::mesh::ConstPartVector{&block2});
  }
};

stk::mesh::EntityId get_random_id(const EntityIdPair& minMaxElemIds)
{
  unsigned numIds = minMaxElemIds.second - minMaxElemIds.first;
  return minMaxElemIds.first + std::rand()%numIds;
}

TEST_F( StressEntityKeyMapping, Timing )
{
  if (get_parallel_size() != 2) return;

  unsigned arbitrarySeed = 1919;
  std::srand(arbitrarySeed);

  const unsigned NUM_RUNS = 5;
  const int numElementsToChange = 10;
  const int numGetEntityQueries = 1000000;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_host_mesh();
    batchTimer.start_batch_timer();
  
    for(int i=0; i<numElementsToChange; ++i) {
      get_bulk().modification_begin();
      add_1_element_to_block_2();
      get_bulk().modification_end();
    }
  
    std::pair<stk::mesh::EntityId,stk::mesh::EntityId> minMaxElemIds = get_min_max_elem_ids_on_local_proc();

    for (int i=0; i<numGetEntityQueries; i++) {
      stk::mesh::EntityId id = get_random_id(minMaxElemIds);
      get_bulk().get_entity(stk::topology::ELEM_RANK, id);
    }

    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(numGetEntityQueries);
}
