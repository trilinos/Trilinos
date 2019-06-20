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

TEST(stk_topology, pyramid_5)
{
  stk::topology t = stk::topology::PYRAMID_5;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),5u);
  EXPECT_EQ(t.num_vertices(),5u);
  EXPECT_EQ(t.num_edges(),8u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_4);

  std::vector<std::vector<unsigned>> gold_edge_node_ordinals = { {0, 1},
                                                                 {1, 2},
                                                                 {2, 3},
                                                                 {3, 0},
                                                                 {0, 4},
                                                                 {1, 4},
                                                                 {2, 4},
                                                                 {3, 4} };
  check_edge_node_ordinals(t, gold_edge_node_ordinals);
  check_edge_nodes(t, gold_edge_node_ordinals);

  std::vector<std::vector<unsigned>> gold_face_node_ordinals = { {0, 1, 4},
                                                                 {1, 2, 4},
                                                                 {2, 3, 4},
                                                                 {0, 4, 3},
                                                                 {0, 3, 2, 1} };
  check_side_node_ordinals(t, gold_face_node_ordinals);
  check_face_node_ordinals(t, gold_face_node_ordinals);
  check_side_nodes(t, gold_face_node_ordinals);
  check_face_nodes(t, gold_face_node_ordinals);

  std::vector<std::vector<unsigned>> gold_permutation_node_ordinals = { {0, 1, 2, 3, 4},
                                                                        {1, 2, 3, 0, 4},
                                                                        {2, 3, 0, 1, 4},
                                                                        {3, 0, 1, 2, 4} };
  check_permutation_node_ordinals(t, gold_permutation_node_ordinals);
  check_permutation_nodes(t, gold_permutation_node_ordinals);

  check_equivalent(t, gold_permutation_node_ordinals);
  check_lexicographical_smallest_permutation(t, gold_permutation_node_ordinals);
}

void check_pyramid_5_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::PYRAMID_5;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),5u);
    NGP_EXPECT_EQ(t.num_vertices(),5u);
    NGP_EXPECT_EQ(t.num_edges(),8u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_4);

    unsigned gold_edge_node_ordinals[8][2] = { {0, 1},
                                               {1, 2},
                                               {2, 3},
                                               {3, 0},
                                               {0, 4},
                                               {1, 4},
                                               {2, 4},
                                               {3, 4} };
    check_edge_node_ordinals_ngp(t, gold_edge_node_ordinals);
    check_edge_nodes_ngp(t, gold_edge_node_ordinals);

    unsigned gold_face_node_ordinals[5][4] = { {0, 1, 4, INVALID},
                                               {1, 2, 4, INVALID},
                                               {2, 3, 4, INVALID},
                                               {0, 4, 3, INVALID},
                                               {0, 3, 2, 1} };
    check_side_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_face_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_side_nodes_ngp(t, gold_face_node_ordinals);
    check_face_nodes_ngp(t, gold_face_node_ordinals);

    unsigned gold_permutation_node_ordinals[4][5] = { {0, 1, 2, 3, 4},
                                                      {1, 2, 3, 0, 4},
                                                      {2, 3, 0, 1, 4},
                                                      {3, 0, 1, 2, 4} };
    check_permutation_node_ordinals_ngp(t, gold_permutation_node_ordinals);
    check_permutation_nodes_ngp(t, gold_permutation_node_ordinals);

    check_equivalent_ngp(t, gold_permutation_node_ordinals);
    check_lexicographical_smallest_permutation_ngp(t, gold_permutation_node_ordinals);
  });
}

NGP_TEST(stk_topology_ngp, pyramid_5)
{
  check_pyramid_5_on_device();
}


TEST(stk_topology, pyramid_13)
{
  stk::topology t = stk::topology::PYRAMID_13;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),13u);
  EXPECT_EQ(t.num_vertices(),5u);
  EXPECT_EQ(t.num_edges(),8u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_8);

  std::vector<std::vector<unsigned>> gold_edge_node_ordinals = { {0, 1, 5},
                                                                 {1, 2, 6},
                                                                 {2, 3, 7},
                                                                 {3, 0, 8},
                                                                 {0, 4, 9},
                                                                 {1, 4, 10},
                                                                 {2, 4, 11},
                                                                 {3, 4, 12} };
  check_edge_node_ordinals(t, gold_edge_node_ordinals);
  check_edge_nodes(t, gold_edge_node_ordinals);

  std::vector<std::vector<unsigned>> gold_face_node_ordinals = { {0, 1, 4, 5, 10,  9},
                                                                 {1, 2, 4, 6, 11, 10},
                                                                 {2, 3, 4, 7, 12, 11},
                                                                 {3, 0, 4, 8,  9, 12},
                                                                 {0, 3, 2, 1,  8,  7, 6, 5} };
  check_side_node_ordinals(t, gold_face_node_ordinals);
  check_face_node_ordinals(t, gold_face_node_ordinals);
  check_side_nodes(t, gold_face_node_ordinals);
  check_face_nodes(t, gold_face_node_ordinals);

  std::vector<std::vector<unsigned>> gold_permutation_node_ordinals = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 11, 12},
    {1, 2, 3, 0, 4, 6, 7, 8, 5, 10, 11, 12,  9},
    {2, 3, 0, 1, 4, 7, 8, 5, 6, 11, 12,  9, 10},
    {3, 0, 1, 2, 4, 8, 5, 6, 7, 12,  9, 10, 11}
  };
  check_permutation_node_ordinals(t, gold_permutation_node_ordinals);
  check_permutation_nodes(t, gold_permutation_node_ordinals);

  check_equivalent(t, gold_permutation_node_ordinals);
  check_lexicographical_smallest_permutation(t, gold_permutation_node_ordinals);
}

void check_pyramid_13_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::PYRAMID_13;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),13u);
    NGP_EXPECT_EQ(t.num_vertices(),5u);
    NGP_EXPECT_EQ(t.num_edges(),8u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_8);

    unsigned gold_edge_node_ordinals[8][3] = { {0, 1, 5},
                                               {1, 2, 6},
                                               {2, 3, 7},
                                               {3, 0, 8},
                                               {0, 4, 9},
                                               {1, 4, 10},
                                               {2, 4, 11},
                                               {3, 4, 12} };
    check_edge_node_ordinals_ngp(t, gold_edge_node_ordinals);
    check_edge_nodes_ngp(t, gold_edge_node_ordinals);

    unsigned gold_face_node_ordinals[5][8] = { {0, 1, 4, 5, 10,  9, INVALID, INVALID},
                                               {1, 2, 4, 6, 11, 10, INVALID, INVALID},
                                               {2, 3, 4, 7, 12, 11, INVALID, INVALID},
                                               {3, 0, 4, 8,  9, 12, INVALID, INVALID},
                                               {0, 3, 2, 1,  8,  7, 6, 5} };
    check_side_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_face_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_side_nodes_ngp(t, gold_face_node_ordinals);
    check_face_nodes_ngp(t, gold_face_node_ordinals);

    unsigned gold_permutation_node_ordinals[4][13] = {
      {0, 1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 11, 12},
      {1, 2, 3, 0, 4, 6, 7, 8, 5, 10, 11, 12,  9},
      {2, 3, 0, 1, 4, 7, 8, 5, 6, 11, 12,  9, 10},
      {3, 0, 1, 2, 4, 8, 5, 6, 7, 12,  9, 10, 11}
    };
    check_permutation_node_ordinals_ngp(t, gold_permutation_node_ordinals);
    check_permutation_nodes_ngp(t, gold_permutation_node_ordinals);

    check_equivalent_ngp(t, gold_permutation_node_ordinals);
    check_lexicographical_smallest_permutation_ngp(t, gold_permutation_node_ordinals);
  });
}

NGP_TEST(stk_topology_ngp, pyramid_13)
{
  check_pyramid_13_on_device();
}


TEST(stk_topology, pyramid_14)
{
  stk::topology t = stk::topology::PYRAMID_14;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),14u);
  EXPECT_EQ(t.num_vertices(),5u);
  EXPECT_EQ(t.num_edges(),8u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_9);

  std::vector<std::vector<unsigned>> gold_edge_node_ordinals = { {0, 1, 5},
                                                                 {1, 2, 6},
                                                                 {2, 3, 7},
                                                                 {3, 0, 8},
                                                                 {0, 4, 9},
                                                                 {1, 4, 10},
                                                                 {2, 4, 11},
                                                                 {3, 4, 12} };
  check_edge_node_ordinals(t, gold_edge_node_ordinals);
  check_edge_nodes(t, gold_edge_node_ordinals);

  std::vector<std::vector<unsigned>> gold_face_node_ordinals = { { 0, 1, 4, 5, 10,  9},
                                                                 { 1, 2, 4, 6, 11, 10},
                                                                 { 2, 3, 4, 7, 12, 11},
                                                                 { 3, 0, 4, 8,  9, 12},
                                                                 { 0, 3, 2, 1,  8,  7, 6, 5, 13} };
  check_side_node_ordinals(t, gold_face_node_ordinals);
  check_face_node_ordinals(t, gold_face_node_ordinals);
  check_side_nodes(t, gold_face_node_ordinals);
  check_face_nodes(t, gold_face_node_ordinals);

  std::vector<std::vector<unsigned>> gold_permutation_node_ordinals = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 11, 12, 13},
    {1, 2, 3, 0, 4, 6, 7, 8, 5, 10, 11, 12,  9, 13},
    {2, 3, 0, 1, 4, 7, 8, 5, 6, 11, 12,  9, 10, 13},
    {3, 0, 1, 2, 4, 8, 5, 6, 7, 12,  9, 10, 11, 13}
  };
  check_permutation_node_ordinals(t, gold_permutation_node_ordinals);
  check_permutation_nodes(t, gold_permutation_node_ordinals);

  check_equivalent(t, gold_permutation_node_ordinals);
  check_lexicographical_smallest_permutation(t, gold_permutation_node_ordinals);
}

void check_pyramid_14_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::PYRAMID_14;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),14u);
    NGP_EXPECT_EQ(t.num_vertices(),5u);
    NGP_EXPECT_EQ(t.num_edges(),8u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_9);

    unsigned gold_edge_node_ordinals[8][3] = { {0, 1, 5},
                                               {1, 2, 6},
                                               {2, 3, 7},
                                               {3, 0, 8},
                                               {0, 4, 9},
                                               {1, 4, 10},
                                               {2, 4, 11},
                                               {3, 4, 12} };
    check_edge_node_ordinals_ngp(t, gold_edge_node_ordinals);
    check_edge_nodes_ngp(t, gold_edge_node_ordinals);

    unsigned gold_face_node_ordinals[5][9] = { { 0, 1, 4, 5, 10,  9, INVALID, INVALID, INVALID},
                                               { 1, 2, 4, 6, 11, 10, INVALID, INVALID, INVALID},
                                               { 2, 3, 4, 7, 12, 11, INVALID, INVALID, INVALID},
                                               { 3, 0, 4, 8,  9, 12, INVALID, INVALID, INVALID},
                                               { 0, 3, 2, 1,  8,  7, 6, 5, 13} };
    check_side_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_face_node_ordinals_ngp(t, gold_face_node_ordinals);
    check_side_nodes_ngp(t, gold_face_node_ordinals);
    check_face_nodes_ngp(t, gold_face_node_ordinals);

    unsigned gold_permutation_node_ordinals[4][14] = {
      {0, 1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 11, 12, 13},
      {1, 2, 3, 0, 4, 6, 7, 8, 5, 10, 11, 12,  9, 13},
      {2, 3, 0, 1, 4, 7, 8, 5, 6, 11, 12,  9, 10, 13},
      {3, 0, 1, 2, 4, 8, 5, 6, 7, 12,  9, 10, 11, 13}
    };
    check_permutation_node_ordinals_ngp(t, gold_permutation_node_ordinals);
    check_permutation_nodes_ngp(t, gold_permutation_node_ordinals);

    check_equivalent_ngp(t, gold_permutation_node_ordinals);
    check_lexicographical_smallest_permutation_ngp(t, gold_permutation_node_ordinals);
  });
}

NGP_TEST(stk_topology_ngp, pyramid_14)
{
  check_pyramid_14_on_device();
}
