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
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_topology/topology.hpp>
#include "topology_test_utils.hpp"

TEST(stk_topology, wedge_6)
{
  stk::topology t = stk::topology::WEDGE_6;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),6u);
  EXPECT_EQ(t.num_vertices(),6u);
  EXPECT_EQ(t.num_edges(),9u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::WEDGE_6);

  EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_4);
  EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_4);
  EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_4);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(4), stk::topology::TRI_3);

  std::vector<std::vector<unsigned>> gold_edge_node_ordinals = { {0, 1},
                                                                 {1, 2},
                                                                 {2, 0},
                                                                 {3, 4},
                                                                 {4, 5},
                                                                 {5, 3},
                                                                 {0, 3},
                                                                 {1, 4},
                                                                 {2, 5} };
  check_edge_node_ordinals(t, gold_edge_node_ordinals);
  check_edge_nodes(t, gold_edge_node_ordinals);

  std::vector<std::vector<unsigned>> gold_face_node_ordinals = { {0, 1, 4, 3},
                                                                 {1, 2, 5, 4},
                                                                 {0, 3, 5, 2},
                                                                 {0, 2, 1},
                                                                 {3, 4, 5} };
  check_side_node_ordinals(t, gold_face_node_ordinals);
  check_face_node_ordinals(t, gold_face_node_ordinals);
  check_side_nodes(t, gold_face_node_ordinals);
  check_face_nodes(t, gold_face_node_ordinals);

  std::vector<std::vector<unsigned>> gold_permutation_node_ordinals = { {0, 1, 2, 3, 4, 5},
                                                                        {1, 2, 0, 4, 5, 3},
                                                                        {2, 0, 1, 5, 3, 4},
                                                                        {3, 5, 4, 0, 2, 1},
                                                                        {5, 4, 3, 2, 1, 0},
                                                                        {4, 3, 5, 1, 0, 2} };
  check_permutation_node_ordinals(t, gold_permutation_node_ordinals);
  check_permutation_nodes(t, gold_permutation_node_ordinals);

  check_equivalent(t, gold_permutation_node_ordinals);
  check_lexicographical_smallest_permutation(t, gold_permutation_node_ordinals);
}

void check_wedge_6_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::WEDGE_6;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),6u);
    NGP_EXPECT_EQ(t.num_vertices(),6u);
    NGP_EXPECT_EQ(t.num_edges(),9u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::WEDGE_6);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_4);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_4);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_4);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::TRI_3);

    unsigned gold_edge_node_ordinals[9][2] = { {0, 1},
                                               {1, 2},
                                               {2, 0},
                                               {3, 4},
                                               {4, 5},
                                               {5, 3},
                                               {0, 3},
                                               {1, 4},
                                               {2, 5} };
    check_edge_node_ordinals_ngp(t, gold_edge_node_ordinals);
    check_edge_nodes_ngp(t, gold_edge_node_ordinals);

    unsigned gold_face_node_ordinals[5][4] = { {0, 1, 4, 3},
                                               {1, 2, 5, 4},
                                               {0, 3, 5, 2},
                                               {0, 2, 1, INVALID},
                                               {3, 4, 5, INVALID} };
    check_side_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_face_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_side_nodes_ngp(t, gold_face_node_ordinals);
    check_face_nodes_ngp(t, gold_face_node_ordinals);

    unsigned gold_permutation_node_ordinals[6][6] = { {0, 1, 2, 3, 4, 5},
                                                      {1, 2, 0, 4, 5, 3},
                                                      {2, 0, 1, 5, 3, 4},
                                                      {3, 5, 4, 0, 2, 1},
                                                      {5, 4, 3, 2, 1, 0},
                                                      {4, 3, 5, 1, 0, 2} };
    check_permutation_node_ordinals_ngp(t, gold_permutation_node_ordinals);
    check_permutation_nodes_ngp(t, gold_permutation_node_ordinals);

    check_equivalent_ngp(t, gold_permutation_node_ordinals);
    check_lexicographical_smallest_permutation_ngp(t, gold_permutation_node_ordinals);
  });
}

NGP_TEST(stk_topology_ngp, wedge_6)
{
  check_wedge_6_on_device();
}

TEST(stk_topology, wedge_12)
{
  stk::topology t = stk::topology::WEDGE_12;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),12u);
  EXPECT_EQ(t.num_vertices(),6u);
  EXPECT_EQ(t.num_edges(),9u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::WEDGE_6);

  EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_6);
  EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_6);
  EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_6);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(4), stk::topology::TRI_6);

  std::vector<std::vector<unsigned>> gold_edge_node_ordinals = { {0, 1,  6},
                                                                 {1, 2,  7},
                                                                 {2, 0,  8},
                                                                 {3, 4,  9},
                                                                 {4, 5, 10},
                                                                 {5, 3, 11},
                                                                 {0, 3},
                                                                 {1, 4},
                                                                 {2, 5} };
  check_edge_node_ordinals(t, gold_edge_node_ordinals);
  check_edge_nodes(t, gold_edge_node_ordinals);

  std::vector<std::vector<unsigned>> gold_face_node_ordinals = { {0, 1, 4, 3,  6, 9},
                                                                 {1, 2, 5, 4,  7, 10},
                                                                 {0, 3, 5, 2,  8, 11},
                                                                 {0, 2, 1,   8, 7, 6},
                                                                 {3, 4, 5,  9, 10, 11} };

  check_side_node_ordinals(t, gold_face_node_ordinals);
  check_face_node_ordinals(t, gold_face_node_ordinals);
  check_side_nodes(t, gold_face_node_ordinals);
  check_face_nodes(t, gold_face_node_ordinals);

  std::vector<std::vector<unsigned>> gold_permutation_node_ordinals = {
    {0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11},
    {1, 2, 0, 4, 5, 3,  7,  8,  6, 10, 11,  9},
    {2, 0, 1, 5, 3, 4,  8,  6,  7, 11,  9, 10},
    {3, 5, 4, 0, 2, 1,  9, 11, 10,  8,  7,  6},
    {5, 4, 3, 2, 1, 0, 11, 10,  9,  7,  6,  8},
    {4, 3, 5, 1, 0, 2, 10,  9, 11,  6,  8,  7}
  };
  std::cout<<"Reminder: we still need to enable permutation for wedge_12"<<std::endl;
  const bool enabledPermutation = false;
  if (enabledPermutation) {
    check_permutation_node_ordinals(t, gold_permutation_node_ordinals);
  }
  check_permutation_nodes(t, gold_permutation_node_ordinals);

  check_equivalent(t, gold_permutation_node_ordinals);
  check_lexicographical_smallest_permutation(t, gold_permutation_node_ordinals);
}

void check_wedge_12_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::WEDGE_12;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),12u);
    NGP_EXPECT_EQ(t.num_vertices(),6u);
    NGP_EXPECT_EQ(t.num_edges(),9u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::WEDGE_6);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_6);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_6);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_6);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::TRI_6);

    unsigned gold_edge_node_ordinals[9][3] = { {0, 1, 6},
                                               {1, 2, 7},
                                               {2, 0, 8},
                                               {3, 4, 9},
                                               {4, 5, 10},
                                               {5, 3, 11},
                                               {0, 3, INVALID},
                                               {1, 4, INVALID},
                                               {2, 5, INVALID} };
    check_edge_node_ordinals_ngp(t, gold_edge_node_ordinals);
    check_edge_nodes_ngp(t, gold_edge_node_ordinals);

    unsigned gold_face_node_ordinals[5][6] = { {0, 1, 4,  3,  6, 9},
                                               {1, 2, 5,  4,  7, 10},
                                               {0, 3, 5,  2,  8, 11},
                                               {0, 2, 1,  8,  7,  6},
                                               {3, 4, 5,  9, 10, 11} };
    check_side_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_face_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_side_nodes_ngp(t, gold_face_node_ordinals);
    check_face_nodes_ngp(t, gold_face_node_ordinals);

    unsigned gold_permutation_node_ordinals[6][12] = {
      {0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11},
      {1, 2, 0, 4, 5, 3,  7,  8,  6, 10, 11,  9},
      {2, 0, 1, 5, 3, 4,  8,  6,  7, 11,  9, 10},
      {3, 5, 4, 0, 2, 1,  9, 11, 10,  8,  7,  6},
      {5, 4, 3, 2, 1, 0, 11, 10,  9,  7,  6,  8},
      {4, 3, 5, 1, 0, 2, 10,  9, 11,  6,  8,  7}
    };
    printf("Reminder: we still need to enable permutation for wedge_12\n");
    const bool enabledPermutation = false;
    if (enabledPermutation) {
      check_permutation_node_ordinals_ngp(t, gold_permutation_node_ordinals);
    }
    check_permutation_nodes_ngp(t, gold_permutation_node_ordinals);

    check_equivalent_ngp(t, gold_permutation_node_ordinals);
    check_lexicographical_smallest_permutation_ngp(t, gold_permutation_node_ordinals);
  });
}

NGP_TEST(stk_topology_ngp, wedge_12)
{
  check_wedge_12_on_device();
}

TEST(stk_topology, wedge_15)
{
  stk::topology t = stk::topology::WEDGE_15;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),15u);
  EXPECT_EQ(t.num_vertices(),6u);
  EXPECT_EQ(t.num_edges(),9u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::WEDGE_6);

  EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_8);
  EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_8);
  EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_8);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(4), stk::topology::TRI_6);

  std::vector<std::vector<unsigned>> gold_edge_node_ordinals = { {0, 1, 6},
                                                                 {1, 2, 7},
                                                                 {2, 0, 8},
                                                                 {3, 4, 12},
                                                                 {4, 5, 13},
                                                                 {5, 3, 14},
                                                                 {0, 3, 9},
                                                                 {1, 4, 10},
                                                                 {2, 5, 11} };
  check_edge_node_ordinals(t, gold_edge_node_ordinals);
  check_edge_nodes(t, gold_edge_node_ordinals);

  std::vector<std::vector<unsigned>> gold_face_node_ordinals = { {0, 1, 4,  3,  6, 10, 12,  9},
                                                                 {1, 2, 5,  4,  7, 11, 13, 10},
                                                                 {0, 3, 5,  2,  9, 14, 11,  8},
                                                                 {0, 2, 1,  8,  7,  6},
                                                                 {3, 4, 5, 12, 13, 14} };
  check_side_node_ordinals(t, gold_face_node_ordinals);
  check_face_node_ordinals(t, gold_face_node_ordinals);
  check_side_nodes(t, gold_face_node_ordinals);
  check_face_nodes(t, gold_face_node_ordinals);

  std::vector<std::vector<unsigned>> gold_permutation_node_ordinals = {
    {0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14},
    {1, 2, 0, 4, 5, 3,  7,  8,  6, 10, 11,  9, 13, 14, 12},
    {2, 0, 1, 5, 3, 4,  8,  6,  7, 11,  9, 10, 14, 12, 13},
    {3, 5, 4, 0, 2, 1, 14, 13, 12,  9, 11, 10,  8,  7,  6},
    {5, 4, 3, 2, 1, 0, 13, 12, 14, 11, 10,  9,  7,  6,  8},
    {4, 3, 5, 1, 0, 2, 12, 14, 13, 10,  9, 11,  6,  8,  7}
  };
  check_permutation_node_ordinals(t, gold_permutation_node_ordinals);
  check_permutation_nodes(t, gold_permutation_node_ordinals);

  check_equivalent(t, gold_permutation_node_ordinals);
  check_lexicographical_smallest_permutation(t, gold_permutation_node_ordinals);
}

void check_wedge_15_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::WEDGE_15;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),15u);
    NGP_EXPECT_EQ(t.num_vertices(),6u);
    NGP_EXPECT_EQ(t.num_edges(),9u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::WEDGE_6);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_8);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_8);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_8);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::TRI_6);

    unsigned gold_edge_node_ordinals[9][3] = { {0, 1, 6},
                                               {1, 2, 7},
                                               {2, 0, 8},
                                               {3, 4, 12},
                                               {4, 5, 13},
                                               {5, 3, 14},
                                               {0, 3, 9},
                                               {1, 4, 10},
                                               {2, 5, 11} };
    check_edge_node_ordinals_ngp(t, gold_edge_node_ordinals);
    check_edge_nodes_ngp(t, gold_edge_node_ordinals);

    unsigned gold_face_node_ordinals[5][8] = { {0, 1, 4,  3,  6, 10, 12,  9},
                                               {1, 2, 5,  4,  7, 11, 13, 10},
                                               {0, 3, 5,  2,  9, 14, 11,  8},
                                               {0, 2, 1,  8,  7,  6, INVALID, INVALID},
                                               {3, 4, 5, 12, 13, 14, INVALID, INVALID} };
    check_side_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_face_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_side_nodes_ngp(t, gold_face_node_ordinals);
    check_face_nodes_ngp(t, gold_face_node_ordinals);

    unsigned gold_permutation_node_ordinals[6][15] = {
      {0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14},
      {1, 2, 0, 4, 5, 3,  7,  8,  6, 10, 11,  9, 13, 14, 12},
      {2, 0, 1, 5, 3, 4,  8,  6,  7, 11,  9, 10, 14, 12, 13},
      {3, 5, 4, 0, 2, 1, 14, 13, 12,  9, 11, 10,  8,  7,  6},
      {5, 4, 3, 2, 1, 0, 13, 12, 14, 11, 10,  9,  7,  6,  8},
      {4, 3, 5, 1, 0, 2, 12, 14, 13, 10,  9, 11,  6,  8,  7}
    };
    check_permutation_node_ordinals_ngp(t, gold_permutation_node_ordinals);
    check_permutation_nodes_ngp(t, gold_permutation_node_ordinals);

    check_equivalent_ngp(t, gold_permutation_node_ordinals);
    check_lexicographical_smallest_permutation_ngp(t, gold_permutation_node_ordinals);
  });
}

NGP_TEST(stk_topology_ngp, wedge_15)
{
  check_wedge_15_on_device();
}


TEST(stk_topology, wedge_18)
{
  stk::topology t = stk::topology::WEDGE_18;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),18u);
  EXPECT_EQ(t.num_vertices(),6u);
  EXPECT_EQ(t.num_edges(),9u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::WEDGE_6);

  EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_9);
  EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_9);
  EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_9);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(4), stk::topology::TRI_6);

  std::vector<std::vector<unsigned>> gold_edge_node_ordinals = { {0, 1, 6},
                                                                 {1, 2, 7},
                                                                 {2, 0, 8},
                                                                 {3, 4, 12},
                                                                 {4, 5, 13},
                                                                 {5, 3, 14},
                                                                 {0, 3, 9},
                                                                 {1, 4, 10},
                                                                 {2, 5, 11} };
  check_edge_node_ordinals(t, gold_edge_node_ordinals);
  check_edge_nodes(t, gold_edge_node_ordinals);

  std::vector<std::vector<unsigned>> gold_face_node_ordinals = { {0, 1, 4,  3,  6, 10, 12,  9, 15},
                                                                 {1, 2, 5,  4,  7, 11, 13, 10, 16},
                                                                 {0, 3, 5,  2,  9, 14, 11,  8, 17},
                                                                 {0, 2, 1,  8,  7,  6},
                                                                 {3, 4, 5, 12, 13, 14} };
  check_side_node_ordinals(t, gold_face_node_ordinals);
  check_face_node_ordinals(t, gold_face_node_ordinals);
  check_side_nodes(t, gold_face_node_ordinals);
  check_face_nodes(t, gold_face_node_ordinals);

  std::vector<std::vector<unsigned>> gold_permutation_node_ordinals = {
    {0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17},
    {1, 2, 0, 4, 5, 3,  7,  8,  6, 10, 11,  9, 13, 14, 12, 16, 17, 15},
    {2, 0, 1, 5, 3, 4,  8,  6,  7, 11,  9, 10, 14, 12, 13, 17, 15, 16},
    {3, 5, 4, 0, 2, 1, 14, 13, 12,  9, 11, 10,  8,  7,  6, 17, 16, 15},
    {5, 4, 3, 2, 1, 0, 13, 12, 14, 11, 10,  9,  7,  6,  8, 16, 15, 17},
    {4, 3, 5, 1, 0, 2, 12, 14, 13, 10,  9, 11,  6,  8,  7, 15, 17, 16}
  };
  check_permutation_node_ordinals(t, gold_permutation_node_ordinals);
  check_permutation_nodes(t, gold_permutation_node_ordinals);

  check_equivalent(t, gold_permutation_node_ordinals);
  check_lexicographical_smallest_permutation(t, gold_permutation_node_ordinals);
}

void check_wedge_18_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::WEDGE_18;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),18u);
    NGP_EXPECT_EQ(t.num_vertices(),6u);
    NGP_EXPECT_EQ(t.num_edges(),9u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::WEDGE_6);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_9);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_9);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_9);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::TRI_6);

    unsigned gold_edge_node_ordinals[9][3] = { {0, 1, 6},
                                               {1, 2, 7},
                                               {2, 0, 8},
                                               {3, 4, 12},
                                               {4, 5, 13},
                                               {5, 3, 14},
                                               {0, 3, 9},
                                               {1, 4, 10},
                                               {2, 5, 11} };
    check_edge_node_ordinals_ngp(t, gold_edge_node_ordinals);
    check_edge_nodes_ngp(t, gold_edge_node_ordinals);

    unsigned gold_face_node_ordinals[5][9] = { {0, 1, 4,  3,  6, 10, 12,  9, 15},
                                               {1, 2, 5,  4,  7, 11, 13, 10, 16},
                                               {0, 3, 5,  2,  9, 14, 11,  8, 17},
                                               {0, 2, 1,  8,  7,  6, INVALID, INVALID, INVALID},
                                               {3, 4, 5, 12, 13, 14, INVALID, INVALID, INVALID} };
    check_side_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_face_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_side_nodes_ngp(t, gold_face_node_ordinals);
    check_face_nodes_ngp(t, gold_face_node_ordinals);

    unsigned gold_permutation_node_ordinals[6][18] = {
      {0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17},
      {1, 2, 0, 4, 5, 3,  7,  8,  6, 10, 11,  9, 13, 14, 12, 16, 17, 15},
      {2, 0, 1, 5, 3, 4,  8,  6,  7, 11,  9, 10, 14, 12, 13, 17, 15, 16},
      {3, 5, 4, 0, 2, 1, 14, 13, 12,  9, 11, 10,  8,  7,  6, 17, 16, 15},
      {5, 4, 3, 2, 1, 0, 13, 12, 14, 11, 10,  9,  7,  6,  8, 16, 15, 17},
      {4, 3, 5, 1, 0, 2, 12, 14, 13, 10,  9, 11,  6,  8,  7, 15, 17, 16}
    };
    check_permutation_node_ordinals_ngp(t, gold_permutation_node_ordinals);
    check_permutation_nodes_ngp(t, gold_permutation_node_ordinals);

    check_equivalent_ngp(t, gold_permutation_node_ordinals);
    check_lexicographical_smallest_permutation_ngp(t, gold_permutation_node_ordinals);
  });
}

NGP_TEST(stk_topology_ngp, wedge_18)
{
  check_wedge_18_on_device();
}
