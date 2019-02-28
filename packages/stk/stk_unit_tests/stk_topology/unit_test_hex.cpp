// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/environment/CPUTime.hpp>  // for cpu_time

TEST(stk_topology, DISABLED_hex8_defined_on_spatial_dimension_perf)
{
  stk::topology hex8 = stk::topology::HEX_8;

  const size_t numIterations = 10000000000;
  size_t result = 0;

  double startTime = stk::cpu_time();
  for (size_t i = 0; i < numIterations; ++i) {
    result += hex8.defined_on_spatial_dimension(i % 4) ? 1 : 0;
  }
  double stopTime = stk::cpu_time();
  std::cout << "result = " << result << ", time = " << stopTime - startTime << "s" << std::endl;
}

TEST(stk_topology, hex8)
{
  stk::topology hex8 = stk::topology::HEX_8;

  EXPECT_TRUE(hex8.is_valid());
  EXPECT_TRUE(hex8.has_homogeneous_faces());
  EXPECT_FALSE(hex8.is_shell());

  EXPECT_EQ(hex8.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(hex8.side_rank(),stk::topology::FACE_RANK);

  const unsigned num_nodes = 8;
  const unsigned num_edges = 12;
  const unsigned num_faces = 6;

  EXPECT_EQ(hex8.num_nodes(),num_nodes);
  EXPECT_EQ(hex8.num_vertices(),num_nodes);
  EXPECT_EQ(hex8.num_edges(),num_edges);
  EXPECT_EQ(hex8.num_faces(),num_faces);

  EXPECT_FALSE(hex8.defined_on_spatial_dimension(1));
  EXPECT_FALSE(hex8.defined_on_spatial_dimension(2));
  EXPECT_TRUE(hex8.defined_on_spatial_dimension(3));

  EXPECT_EQ(hex8.base(),stk::topology::HEX_8);

}

void check_hex_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology hex8 = stk::topology::HEX_8;

    NGP_EXPECT_TRUE(hex8.is_valid());
    NGP_EXPECT_TRUE(hex8.has_homogeneous_faces());
    NGP_EXPECT_FALSE(hex8.is_shell());

    NGP_EXPECT_EQ(hex8.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(hex8.side_rank(),stk::topology::FACE_RANK);

    const unsigned num_nodes = 8;
    const unsigned num_edges = 12;
    const unsigned num_faces = 6;

    NGP_EXPECT_EQ(hex8.num_nodes(),num_nodes);
    NGP_EXPECT_EQ(hex8.num_vertices(),num_nodes);
    NGP_EXPECT_EQ(hex8.num_edges(),num_edges);
    NGP_EXPECT_EQ(hex8.num_faces(),num_faces);

    NGP_EXPECT_FALSE(hex8.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(hex8.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(hex8.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(hex8.base(),stk::topology::HEX_8);

  });
}

NGP_TEST(stk_topology_ngp, hex8)
{
  check_hex_on_device();
}

