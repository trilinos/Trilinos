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

#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>

#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <sstream>

//This very simple test will visit all local elements and traverse the
//element's node relations. It will return the number of nodes visited.
//The purpose of this test is to stress the relation-traversal for a
//performance test.
size_t do_stk_node_rel_test(stk::mesh::BulkData& bulk)
{
  using namespace stk::mesh;

  const MetaData& meta = bulk.mesh_meta_data();

  Selector local = meta.locally_owned_part();

  BucketVector const& buckets = bulk.get_buckets(stk::topology::ELEMENT_RANK, local);

  size_t nodes_visited = 0;

  for(size_t ib=0; ib<buckets.size(); ++ib) {
    const Bucket& b = *buckets[ib];

    for(size_t i=0; i<b.size(); ++i) {
      Entity const *node_itr = b.begin_nodes(i);
      Entity const *nodes_end = b.end_nodes(i);

      for (; node_itr != nodes_end; ++node_itr)
      {
        Entity node = *node_itr;
        EXPECT_TRUE(bulk.is_valid(node));
        ++nodes_visited;
      }
    }
  }
  return nodes_visited;
}

TEST(node_rels, node_rels)
{
  //vector of mesh-dimensions holds the number of elements in each dimension.
  //Hard-wired to 3. This test can run with spatial-dimension less than 3,
  //(if generated-mesh can do that) but not greater than 3.
  std::vector<int> mesh_dims(3);
#ifndef NDEBUG
  mesh_dims[0]=50; //num_elems_x
  mesh_dims[1]=50; //num_elems_y
  mesh_dims[2]=50; //num_elems_z
#else
  mesh_dims[0]=150; //num_elems_x
  mesh_dims[1]=150; //num_elems_y
  mesh_dims[2]=150; //num_elems_z
#endif

  std::ostringstream oss;
  oss << mesh_dims[0] << "x" << mesh_dims[1] << "x" << mesh_dims[2] << "|sideset:xXyYzZ";

  double start_time = stk::cpu_time();

  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, oss.str());
  fixture.commit();

  double mesh_create_time = stk::cpu_time() - start_time;

  //compute total number of elements in the mesh:
  const size_t spatial_dim = fixture.getMetaData().spatial_dimension();
  size_t num_elems = mesh_dims[0];
  for(size_t d=1; d<spatial_dim; ++d) {
    num_elems *= mesh_dims[d];
  }

  std::cout << "num_elems: " << num_elems << std::endl;
  std::cout << "sizeof(stk::mesh::Relation): " << sizeof(stk::mesh::Relation) << std::endl;
  std::cout << "sizeof(stk::mesh::Entity): " << sizeof(stk::mesh::Entity) << std::endl;

  start_time = stk::cpu_time();

  const int num_iters = 50;
  for(int t=0; t<num_iters; ++t) {

    size_t nodes_visited = do_stk_node_rel_test(fixture.getBulkData());
    size_t expected_nodes_visited = num_elems*8;

    EXPECT_EQ(nodes_visited, expected_nodes_visited);
  }

  double traverse_time = stk::cpu_time() - start_time;
  double total_time    = mesh_create_time + traverse_time;

  static const int NUM_TIMERS = 3;
  const double timers[NUM_TIMERS] = {mesh_create_time, traverse_time, total_time};
  const char* timer_names[NUM_TIMERS] = {"Create mesh", "Traverse", "Total time"};

  stk::print_timers_and_memory(timer_names, timers, NUM_TIMERS);

  stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
}
