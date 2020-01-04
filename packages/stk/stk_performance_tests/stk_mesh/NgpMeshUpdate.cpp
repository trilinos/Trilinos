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
#include <stk_ngp/NgpMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

class NgpMeshUpdate : public stk::unit_test_util::MeshFixture {
public:
  NgpMeshUpdate()
    : stk::unit_test_util::MeshFixture(),
      cumulativeTime(0.0),
      ngpMeshHwm(0),
      newPartName("block2"),
      ghostingName("testGhosting")
  { }

protected:
  void setup_host_mesh()
  {
    setup_mesh("generated:1x1x1000000", stk::mesh::BulkData::NO_AUTO_AURA);
    numElements = 1000000;
    get_meta().declare_part(newPartName);

    get_bulk().modification_begin();
    ghosting = &get_bulk().create_ghosting(ghostingName);
    get_bulk().modification_end();
  }

  void change_element_part_membership(int cycle)
  {
    stk::mesh::EntityId elemId = cycle+1;
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    const stk::mesh::Part* part = get_meta().get_part(newPartName);

    get_bulk().modification_begin();
    get_bulk().change_entity_parts<stk::mesh::ConstPartVector>(elem, {part});
    get_bulk().modification_end();
  }

  void create_entity(int cycle)
  {
    stk::mesh::EntityId elemId = numElements + cycle + 1;

    get_bulk().modification_begin();
    get_bulk().declare_element(elemId);
    get_bulk().modification_end();
  }

  void ghost_element(int cycle)
  {
    stk::mesh::EntityId firstLocalElemId = get_parallel_rank()*numElements/2 + 1;
    stk::mesh::EntityId elemId = firstLocalElemId + cycle;
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    int ghostRank = 1 - get_parallel_rank();
    stk::mesh::EntityProcVec toGhost = {{elem, ghostRank}};

    get_bulk().modification_begin();
    get_bulk().change_ghosting(*ghosting, toGhost);
    get_bulk().modification_end();
  }

  void time_get_ngp_mesh()
  {
    iterationStartTime = stk::wall_time();
    bulkDataHwm = stk::get_max_hwm_across_procs(get_comm());

    ngp::StaticMesh ngpMesh(get_bulk());

    cumulativeTime += stk::wall_dtime(iterationStartTime);
    size_t currentNgpMeshHwm = stk::get_max_hwm_across_procs(get_comm()) - bulkDataHwm;
    ngpMeshHwm = std::max(currentNgpMeshHwm, ngpMeshHwm);
  }

  void print_timing()
  {
    double timeAll = stk::get_global_sum(get_comm(), cumulativeTime);
    stk::print_stats_for_performance_compare(std::cout, timeAll, ngpMeshHwm, 1, get_comm());
  }

private:
  double iterationStartTime;
  double cumulativeTime;
  size_t bulkDataHwm;
  size_t ngpMeshHwm;

  std::string newPartName;
  int numElements;
  std::string ghostingName;
  stk::mesh::Ghosting* ghosting;
};

TEST_F( NgpMeshUpdate, PartMembership )
{
  if (get_parallel_size() != 1) return;

  setup_host_mesh();
  for (int i=0; i<400; i++) {
    change_element_part_membership(i);
    time_get_ngp_mesh();
  }
  print_timing();
}

TEST_F( NgpMeshUpdate, EntityCreation )
{
  if (get_parallel_size() != 1) return;

  setup_host_mesh();
  for (int i=0; i<400; i++) {
    create_entity(i);
    time_get_ngp_mesh();
  }
  print_timing();
}

TEST_F( NgpMeshUpdate, Ghosting )
{
  if (get_parallel_size() != 2) return;

  setup_host_mesh();
  for (int i=0; i<400; i++) {
    ghost_element(i);
    time_get_ngp_mesh();
  }
  print_timing();
}
