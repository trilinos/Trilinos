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

#include <stk_util/environment/CPUTime.hpp>
#include <gtest/gtest.h>
#include <stk_util/environment/perf_util.hpp>

#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>

#include "stk_unit_test_utils/stk_mesh_fixtures/Gear.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/GearsFixture.hpp"
#include <sstream>

namespace stk {
namespace performance_tests {

namespace {

//This very simple test will visit all local elements, gather coordinates,
//compute element-centroid (simple average of nodal coords) for each, and
//store the sum of the centroids in sum_centroid.
void do_stk_gather_gears_test(stk::mesh::BulkData& bulk, std::vector<double>& sum_centroid)
{
  using namespace stk::mesh;
  typedef Field<double> VectorField;

  const MetaData& meta = bulk.mesh_meta_data();
  const unsigned spatial_dim = meta.spatial_dimension();
  for(unsigned d=0; d<spatial_dim; ++d) sum_centroid[d] = 0;

  std::vector<double> elem_centroid(spatial_dim, 0);

  const VectorField * coord_field = meta.get_field<double>(stk::topology::NODE_RANK, "coordinates");
  STK_ThrowAssert(coord_field != nullptr);

  Selector local = meta.locally_owned_part();

  BucketVector const& buckets = bulk.get_buckets(stk::topology::ELEMENT_RANK, local);

  std::vector<double> elem_node_coords;

  for(size_t ib=0; ib<buckets.size(); ++ib) {
    const Bucket& b = *buckets[ib];

    for(size_t i=0; i<b.size(); ++i) {

      Entity const *node_rels = b.begin_nodes(i);
      size_t num_nodes = b.num_nodes(i);
      size_t len = num_nodes*spatial_dim;
      if (elem_node_coords.size() != len) elem_node_coords.resize(len);

      //here's the gather:

      unsigned offset = 0;
      for(size_t n=0; n<num_nodes; ++n) {
        Entity node = node_rels[n];
        double* node_coords = stk::mesh::field_data(*coord_field, node);
        elem_node_coords[offset++] = node_coords[0];
        elem_node_coords[offset++] = node_coords[1];
        elem_node_coords[offset++] = node_coords[2];
      }

      stk::performance_tests::calculate_centroid_3d(num_nodes, &elem_node_coords[0], &elem_centroid[0]);

      //add this element-centroid to the sum_centroid vector, and
      //re-zero the element-centroid vector:
      sum_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
      sum_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
      sum_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
    }
  }
}

} // empty namespace

TEST(gather_gears, gather_gears)
{
  stk::mesh::fixtures::GearsFixture fixture(MPI_COMM_WORLD, 1,
      stk::mesh::fixtures::GearParams(0.01, 0.4, 1.5, -0.4, 0.4));
  fixture.meta_data.commit();

  double start_time = stk::cpu_time();

  fixture.generate_mesh();

  double mesh_create_time = stk::cpu_time() - start_time;

  const size_t spatial_dim = fixture.meta_data.spatial_dimension();

  std::vector<double> avg_centroid(spatial_dim, 0);

  start_time = stk::cpu_time();

  const int num_iters = 150;
  for(int t=0; t<num_iters; ++t) {

    do_stk_gather_gears_test(fixture.bulk_data, avg_centroid);

    double tolerance = 1.e-6;

    for(size_t d=0; d<spatial_dim; ++d) {
      double expected = 0;
      EXPECT_LT(std::abs(avg_centroid[d] - expected), tolerance);
    }
  }

  double gather_time = stk::cpu_time() - start_time;

  std::cout << "Gear: ";
  std::cout << "\tNum Nodes: " << fixture.get_gear(0).num_nodes;
  std::cout << "\tNum Elements: " << fixture.get_gear(0).num_elements << std::endl;

  double total_time = mesh_create_time + gather_time;

  static const int NUM_TIMERS = 3;
  const double timers[NUM_TIMERS] = {mesh_create_time, gather_time, total_time};
  const char* timer_names[NUM_TIMERS] = {"Create mesh", "Gather", "Total time"};

  stk::print_timers_and_memory(timer_names, timers, NUM_TIMERS);

  stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
}

} // namespace performance_tests
} // namespace stk
