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
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/timer.hpp>

class NgpMeshChangeElementPartMembership : public stk::unit_test_util::MeshFixture
{
public:
  NgpMeshChangeElementPartMembership()
    : stk::unit_test_util::MeshFixture(),
      newPartName("block2")
  { }

  void setup_host_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
#ifdef NDEBUG
    numElements = 750000;
    setup_mesh("generated:300x250x10", auraOption);
#else
    numElements = 10000;
    setup_mesh("generated:10x10x100", auraOption);
#endif
    get_meta().declare_part(newPartName);
  }

  void change_element_part_membership(int cycle)
  {
    get_bulk().modification_begin();
    get_bulk().change_entity_parts<stk::mesh::ConstPartVector>(get_element(cycle), {get_part()});
    get_bulk().modification_end();
    stk::mesh::get_updated_ngp_mesh(get_bulk());
  }

  void batch_change_element_part_membership(int cycle)
  {
    get_bulk().batch_change_entity_parts(stk::mesh::EntityVector{get_element(cycle)},
                                         stk::mesh::PartVector{get_part()}, {});
    stk::mesh::get_updated_ngp_mesh(get_bulk());
  }

private:
  stk::mesh::Entity get_element(int cycle)
  {
    stk::mesh::EntityId firstLocalElemId = get_parallel_rank()*numElements/2 + 1;
    stk::mesh::EntityId elemId = firstLocalElemId + cycle;
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    return elem;
  }

  stk::mesh::Part* get_part()
  {
    return get_meta().get_part(newPartName);
  }

  std::string newPartName;
  unsigned numElements;
};

class NgpMeshCreateEntity : public stk::unit_test_util::MeshFixture
{
public:
  NgpMeshCreateEntity()
    : stk::unit_test_util::MeshFixture(),
      numElements(1000000)
  { }

  void setup_host_mesh()
  {
    setup_mesh("generated:100x100x100", stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void create_entity(int cycle)
  {
    get_bulk().modification_begin();
    get_bulk().declare_element(get_new_entity_id(cycle));
    get_bulk().modification_end();
    stk::mesh::get_updated_ngp_mesh(get_bulk());
  }

private:
  stk::mesh::EntityId get_new_entity_id(int cycle)
  {
    return numElements + cycle + 1;
  }

  int numElements;
};

class NgpMeshGhosting : public stk::unit_test_util::MeshFixture
{
public:
  NgpMeshGhosting()
    : stk::unit_test_util::MeshFixture(),
      ghostingName("testGhosting")
  { }

protected:
  void setup_host_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
#ifdef NDEBUG
    numElements = 1000000;
    setup_mesh("generated:400x250x10", auraOption);
#else
    numElements = 10000;
    setup_mesh("generated:10x10x100", auraOption);
#endif
    get_bulk().modification_begin();
    ghosting = &get_bulk().create_ghosting(ghostingName);
    get_bulk().modification_end();
  }

  void ghost_element(int cycle)
  {
    get_bulk().modification_begin();
    get_bulk().change_ghosting(*ghosting, element_to_ghost(cycle));
    get_bulk().modification_end();
    stk::mesh::get_updated_ngp_mesh(get_bulk());
  }

private:
  stk::mesh::EntityProcVec element_to_ghost(int cycle)
  {
    stk::mesh::EntityId firstLocalElemId = get_parallel_rank()*numElements/2 + 1;
    stk::mesh::EntityId elemId = firstLocalElemId + cycle;
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    int ghostRank = 1 - get_parallel_rank();
    return {{elem, ghostRank}};
  }

  int numElements;
  std::string ghostingName;
  stk::mesh::Ghosting* ghosting;
};

TEST_F( NgpMeshChangeElementPartMembership, Timing )
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 200;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_host_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    batchTimer.start_batch_timer();
  
    for (int i = 0; i < NUM_ITERS; i++) {
      change_element_part_membership(i);
    }
    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F( NgpMeshChangeElementPartMembership, TimingWithAura )
{
  if (get_parallel_size() != 2) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 50;
    
  stk::parallel_machine_barrier(get_comm());

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_host_mesh(stk::mesh::BulkData::AUTO_AURA);
    batchTimer.start_batch_timer();
  
    for (int i = 0; i < NUM_ITERS; i++) {
      change_element_part_membership(i);
    }
    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F( NgpMeshChangeElementPartMembership, TimingBatch )
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 200;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    setup_host_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    for (int i = 0; i < NUM_ITERS; i++) {
      batch_change_element_part_membership(i);
    }
    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F( NgpMeshCreateEntity, Timing )
{
  if (get_parallel_size() != 1) return;

  const unsigned NUM_RUNS = 5;
  #ifdef STK_ENABLE_GPU
  const int NUM_ITERS = 100;
  #else
  const int NUM_ITERS = 5000;
  #endif
  const int NUM_FAKE_ITERS = 5000;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_host_mesh();
    batchTimer.start_batch_timer();

    for (int i=0; i<NUM_ITERS; i++) {
      create_entity(i);
    }
    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_FAKE_ITERS);
}

TEST_F( NgpMeshGhosting, Timing )
{
  if (get_parallel_size() != 2) return;

  std::string perfCheck = stk::unit_test_util::get_option("-perf_check", "PERF_CHECK");
#ifdef NDEBUG
  const int NUM_INNER_ITERS = (perfCheck=="NO_PERF_CHECK" ? 1 : 100);
#else
  const int NUM_INNER_ITERS = 1;
#endif

  const unsigned NUM_RUNS = 5;
  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    setup_host_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  
    for (int i = 0; i < NUM_INNER_ITERS; i++) {
      ghost_element(i);
    }

    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_INNER_ITERS);
}

TEST_F( NgpMeshGhosting, TimingWithAura )
{
  if (get_parallel_size() != 2) return;

#ifdef NDEBUG
  const int NUM_INNER_ITERS = 50;
#else
  const int NUM_INNER_ITERS = 1;
#endif

  const unsigned NUM_RUNS = 5;
  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_host_mesh(stk::mesh::BulkData::AUTO_AURA);
    batchTimer.start_batch_timer();

    for (int i = 0; i < NUM_INNER_ITERS; i++) {
      ghost_element(i);
    }

    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_INNER_ITERS);
}
