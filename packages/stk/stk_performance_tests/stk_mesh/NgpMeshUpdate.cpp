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
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

class Timer
{
public:
  Timer(MPI_Comm comm)
    : communicator(comm),
      iterationStartTime(0.0),
      cumulativeTime(0.0),
      iterationStartHwm(0),
      meshOperationHwm(1)
  { }

  void start_timing()
  {
    iterationStartTime = stk::wall_time();
    iterationStartHwm = stk::get_max_hwm_across_procs(communicator);
  }

  void update_timing()
  {
    cumulativeTime += stk::wall_dtime(iterationStartTime);
    size_t currentHwm = stk::get_max_hwm_across_procs(communicator) - iterationStartHwm;
    meshOperationHwm = std::max(currentHwm, meshOperationHwm);
  }

  void print_timing()
  {
    double timeAll = stk::get_global_sum(communicator, cumulativeTime);
    stk::print_stats_for_performance_compare(std::cout, timeAll, meshOperationHwm, 1, communicator);
  }

private:
  MPI_Comm communicator;
  double iterationStartTime;
  double cumulativeTime;
  size_t iterationStartHwm;
  size_t meshOperationHwm;
};

class NgpMeshChangeElementPartMembership : public stk::unit_test_util::MeshFixture
{
public:
  NgpMeshChangeElementPartMembership()
    : stk::unit_test_util::MeshFixture(),
      newPartName("block2")
  { }

  void setup_host_mesh()
  {
    setup_mesh("generated:1x1x1000000", stk::mesh::BulkData::NO_AUTO_AURA);
    get_meta().declare_part(newPartName);
  }

  void change_element_part_membership(int cycle)
  {
    get_bulk().modification_begin();
    get_bulk().change_entity_parts<stk::mesh::ConstPartVector>(get_element(cycle), {get_part()});
    get_bulk().modification_end();
    get_bulk().get_ngp_mesh();
  }

private:
  stk::mesh::Entity get_element(int cycle)
  {
    stk::mesh::EntityId elemId = cycle+1;
    return get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
  }

  const stk::mesh::Part* get_part()
  {
    return get_meta().get_part(newPartName);
  }

  std::string newPartName;
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
    setup_mesh("generated:1x1x1000000", stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void create_entity(int cycle)
  {
    get_bulk().modification_begin();
    get_bulk().declare_element(get_new_entity_id(cycle));
    get_bulk().modification_end();
    get_bulk().get_ngp_mesh();
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
      numElements(1000000),
      ghostingName("testGhosting")
  { }

protected:
  void setup_host_mesh()
  {
    setup_mesh("generated:1x1x1000000", stk::mesh::BulkData::NO_AUTO_AURA);
    get_bulk().modification_begin();
    ghosting = &get_bulk().create_ghosting(ghostingName);
    get_bulk().modification_end();
  }

  void ghost_element(int cycle)
  {
    get_bulk().modification_begin();
    get_bulk().change_ghosting(*ghosting, element_to_ghost(cycle));
    get_bulk().modification_end();
    get_bulk().get_ngp_mesh();
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
  if (get_parallel_size() != 1) return;

  Timer timer(get_comm());
  setup_host_mesh();

  for (int i=0; i<50; i++) {
    timer.start_timing();
    change_element_part_membership(i);
    timer.update_timing();
  }
  timer.print_timing();
}

TEST_F( NgpMeshCreateEntity, Timing )
{
  if (get_parallel_size() != 1) return;

  Timer timer(get_comm());
  setup_host_mesh();

  for (int i=0; i<50; i++) {
    timer.start_timing();
    create_entity(i);
    timer.update_timing();
  }
  timer.print_timing();
}

TEST_F( NgpMeshGhosting, Timing )
{
  if (get_parallel_size() != 2) return;

  Timer timer(get_comm());
  setup_host_mesh();

  for (int i=0; i<50; i++) {
    timer.start_timing();
    ghost_element(i);
    timer.update_timing();
  }
  timer.print_timing();
}
