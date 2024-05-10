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

#include <sstream>
#include <gtest/gtest.h>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>
#include <stk_util/environment/perf_util.hpp>

#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>
#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>

#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/Comm.hpp>

namespace stk {
namespace performance_tests {

TEST(hex_edges, hex_edges)
{
  //vector of mesh-dimensions holds the number of elements in each dimension.
  //Hard-wired to 3. This test can run with spatial-dimension less than 3,
  //(if generated-mesh can do that) but not greater than 3.
  std::vector<int> mesh_dims(3);
  int proc = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  int numprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

#ifndef NDEBUG
  mesh_dims[0]=50; //num_elems_x
  mesh_dims[1]=50; //num_elems_y
  mesh_dims[2]=50; //num_elems_z
#else
  mesh_dims[0]=200; //num_elems_x
  mesh_dims[1]=200; //num_elems_y
  mesh_dims[2]=2*numprocs; //num_elems_z
#endif

  std::ostringstream oss;
  oss << mesh_dims[0] << "x" << mesh_dims[1] << "x" << mesh_dims[2];

  double start_time = stk::cpu_time();

  size_t num_nodes = mesh_dims[0]+1;
  num_nodes *= mesh_dims[1]+1;
  num_nodes *= mesh_dims[2]+1;

  //if num_nodes is greater than the int limit (2.1B) then we need to
  //use the 64-bit api in the IO (which lies beneath the gmesh fixture).
  bool use_64bit_IO_api = num_nodes > 2150000000;

  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, oss.str(), use_64bit_IO_api);
  fixture.commit();

  double mesh_create_time = stk::cpu_time() - start_time;

  start_time = stk::cpu_time();

  stk::mesh::create_edges(fixture.getBulkData());

  double create_edges_time = stk::cpu_time() - start_time;
  double total_time = mesh_create_time + create_edges_time;

  //quick check to make sure all nodes got induced into the edge part.
  const stk::mesh::Part& edgePart = fixture.getBulkData().mesh_meta_data().get_topology_root_part(stk::topology::LINE_2);
  const stk::mesh::BucketVector& nodeBuckets = fixture.getBulkData().get_buckets(stk::topology::NODE_RANK, edgePart);
  const stk::mesh::BucketVector& allNodeBuckets = fixture.getBulkData().buckets(stk::topology::NODE_RANK);
  EXPECT_EQ(allNodeBuckets.size(), nodeBuckets.size());

  std::vector<size_t> mesh_counts;
  stk::mesh::comm_mesh_counts(fixture.getBulkData(), mesh_counts);

  if (proc == 0) {
    std::cout<< "\nnum nodes: "<<mesh_counts[stk::topology::NODE_RANK]<<std::endl;
    std::cout<< "num edges: "<<mesh_counts[stk::topology::EDGE_RANK]<<std::endl;
    std::cout<< "num faces: "<<mesh_counts[stk::topology::FACE_RANK]<<std::endl;
    std::cout<< "num elems: "<<mesh_counts[stk::topology::ELEMENT_RANK]<<std::endl;
  }

  size_t now=0, hwm=0;
  stk::get_memory_usage(now,hwm);
  size_t global_hwm=0;
  MPI_Allreduce(&hwm, &global_hwm, 1, MPI_LONG_LONG, MPI_MAX, fixture.getBulkData().parallel());

  if (proc == 0) {
    static const int NUM_TIMERS = 3;
    const double timers[NUM_TIMERS] = {mesh_create_time, create_edges_time, total_time};
    const char* timer_names[NUM_TIMERS] = {"Create mesh", "Create edges", "Total time"};

    stk::print_timers_and_memory(timer_names, timers, NUM_TIMERS);
    std::cout<<"Global HWM: "<<stk::human_bytes(global_hwm)<<std::endl;
  }

  stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
}

TEST(hex_edges, minimal_hex_edges)
{
  //vector of mesh-dimensions holds the number of elements in each dimension.
  //Hard-wired to 3. This test can run with spatial-dimension less than 3,
  //(if generated-mesh can do that) but not greater than 3.
  std::vector<int> mesh_dims(3);
  int proc = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  int numprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

#ifndef NDEBUG
  mesh_dims[0]=50; //num_elems_x
  mesh_dims[1]=50; //num_elems_y
  mesh_dims[2]=50; //num_elems_z
#else
  mesh_dims[0]=100; //num_elems_x
  mesh_dims[1]=100; //num_elems_y
  mesh_dims[2]=100*numprocs; //num_elems_z
#endif

  std::ostringstream oss;
  oss << mesh_dims[0] << "x" << mesh_dims[1] << "x" << mesh_dims[2];

  double start_time = stk::cpu_time();

  size_t num_nodes = mesh_dims[0]+1;
  num_nodes *= mesh_dims[1]+1;
  num_nodes *= mesh_dims[2]+1;

  //if num_nodes is greater than the int limit (2.1B) then we need to
  //use the 64-bit api in the IO (which lies beneath the gmesh fixture).
  bool use_64bit_IO_api = num_nodes > 2150000000;

  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, oss.str(), use_64bit_IO_api);
  fixture.commit();

  double mesh_create_time = stk::cpu_time() - start_time;

  start_time = stk::cpu_time();

  stk::mesh::create_edges(fixture.getBulkData());

  double create_edges_time = stk::cpu_time() - start_time;
  double total_time = mesh_create_time + create_edges_time;

  std::vector<size_t> mesh_counts;
  stk::mesh::comm_mesh_counts(fixture.getBulkData(), mesh_counts);

  if (proc == 0) {
    std::cout<< "\nnum nodes: "<<mesh_counts[stk::topology::NODE_RANK]<<std::endl;
    std::cout<< "num edges: "<<mesh_counts[stk::topology::EDGE_RANK]<<std::endl;
    std::cout<< "num faces: "<<mesh_counts[stk::topology::FACE_RANK]<<std::endl;
    std::cout<< "num elems: "<<mesh_counts[stk::topology::ELEMENT_RANK]<<std::endl;
  }

  size_t now=0, hwm=0;
  stk::get_memory_usage(now,hwm);
  size_t global_hwm=0;
  MPI_Allreduce(&hwm, &global_hwm, 1, MPI_LONG_LONG, MPI_MAX, fixture.getBulkData().parallel());

  if (proc == 0) {
    static const int NUM_TIMERS = 3;
    const double timers[NUM_TIMERS] = {mesh_create_time, create_edges_time, total_time};
    const char* timer_names[NUM_TIMERS] = {"Create mesh", "Create edges", "Total time"};

    stk::print_timers_and_memory(timer_names, timers, NUM_TIMERS);
    std::cout<<"Global HWM: "<<stk::human_bytes(global_hwm)<<std::endl;
  }

  stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
}

} //namespace performance_tests
} //namespace stk
