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

#ifndef __IBMCPP__
#include <gtest/gtest.h>
#include <unordered_map>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/perf_util.hpp>

#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>
#include <stk_performance_tests/stk_mesh/hex_refine_info.hpp>

#include <iostream>

using namespace stk::mesh;

namespace stk {
namespace performance_tests {

namespace  {

void create_entities( BulkData & bulk,
                      PartVector& node_parts,
                      PartVector& elem_parts,
                      // FIXME Part& active_elements_part,
                      HexRefineInfo& refine_info)
{
  HexRefineInfo refine_info_half(refine_info.m_level-1, refine_info.m_nx, refine_info.m_ny, refine_info.m_nz);

  unsigned eid_start = 1 + refine_info.elem_id_offset();
  unsigned eid_end = eid_start + refine_info.num_elems();

  unsigned nid_start = 1 + refine_info.node_id_offset();
  unsigned nid_end = nid_start + refine_info.num_nodes();

  std::unordered_map<unsigned, Entity> node_map;

  for(unsigned nid=nid_start; nid<nid_end; ++nid) {
    Entity e = bulk.declare_node(nid, node_parts);
    node_map[nid] = e;
  }

  EntityIdVector elem_node(8);
  for (unsigned entity_id = eid_start; entity_id < eid_end; ++entity_id)  {
    unsigned ix = 0, iy = 0, iz = 0;
    refine_info.elem_x_y_z(entity_id, ix, iy, iz);
    EntityId ie_check = refine_info.elem_id(ix, iy, iz);
    EXPECT_EQ(ie_check, entity_id);

    elem_node[0] = refine_info.node_id( ix   , iy   , iz   );
    elem_node[1] = refine_info.node_id( ix+1 , iy   , iz   );
    elem_node[2] = refine_info.node_id( ix+1 , iy   , iz+1 );
    elem_node[3] = refine_info.node_id( ix   , iy   , iz+1 );
    elem_node[4] = refine_info.node_id( ix   , iy+1 , iz   );
    elem_node[5] = refine_info.node_id( ix+1 , iy+1 , iz   );
    elem_node[6] = refine_info.node_id( ix+1 , iy+1 , iz+1 );
    elem_node[7] = refine_info.node_id( ix   , iy+1 , iz+1 );

    // check if a parent node
    for (unsigned i = 0; i<8; ++i) {
      unsigned ixn = 0, iyn = 0, izn = 0;
      refine_info.node_x_y_z(elem_node[i], ixn, iyn, izn);
      EntityId in_check = refine_info.node_id(ixn, iyn, izn);
      EXPECT_EQ(in_check, elem_node[i]);

      if (
        ((ixn - 1) % 2 == 0) &&
        ((iyn - 1) % 2 == 0) &&
        ((izn - 1) % 2 == 0))
      {
        elem_node[i] = refine_info_half.node_id(ixn/2, iyn/2, izn/2);
      }
    }

    stk::mesh::declare_element( bulk, elem_parts, entity_id , elem_node);
  }
}

} // empty namespace

TEST( hex_refine, hex_refine)
{
  double start_time = stk::cpu_time();
  const int max_levels = 2;
  unsigned nn = 100/(1 << (max_levels-1));
  unsigned ex=nn, ey=nn, ez=nn;

  fixtures::HexFixture fixture(MPI_COMM_WORLD, ex, ey, ez);
  fixture.m_meta.commit();
  fixture.generate_mesh();
  double create_time = stk::cpu_time() - start_time;

  std::vector<double> avg_centroid(3, 0.0);

  start_time = stk::cpu_time();

  for (int level=1; level <= max_levels; ++level) {

    HexRefineInfo refine_info(level, ex, ey, ez);

    //Selector hex_elem_selector(fixture.m_hex_part & !fixture.m_node_part);

    unsigned num_elems_new = refine_info.num_elems();
    std::cout << "num_elems_new for level = " << level << " = " << num_elems_new << std::endl;
    fixture.m_bulk_data.modification_begin();
    create_entities(fixture.m_bulk_data, fixture.m_node_parts, fixture.m_elem_parts, refine_info);
    fixture.m_bulk_data.modification_end();
  }

  double refine_time = stk::cpu_time() - start_time;
  double total_time = refine_time + create_time;

  static const int NUM_TIMERS = 3;
  const double timers[NUM_TIMERS] = {create_time, refine_time, total_time};
  const char* timer_names[NUM_TIMERS] = {"Create mesh", "Refine", "Total time"};

  stk::print_timers_and_memory(timer_names, timers, NUM_TIMERS);

  stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
}

} // namespace performance_tests
} // namespace std

#endif
