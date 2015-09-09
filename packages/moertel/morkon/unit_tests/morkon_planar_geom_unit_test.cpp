/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2015) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */


#include <gtest/gtest.h>

#include <mrk_planar_geom.hpp>
#include <morkon_unit_test_utils.hpp>

TEST(morkon_planar_geom, geom_primitives)
{
  using namespace morkon_exp;

  double seg_head[2] = {0, 1};
  double seg_tail[2] = {4, 0};

  double seg_midpoint[2];
  add_2d(seg_tail, seg_head, seg_midpoint);
  EXPECT_DOUBLE_EQ(4, seg_midpoint[0]);
  EXPECT_DOUBLE_EQ(1, seg_midpoint[1]);

  scale_2d(0.5, seg_midpoint);
  EXPECT_DOUBLE_EQ(2, seg_midpoint[0]);
  EXPECT_DOUBLE_EQ(0.5, seg_midpoint[1]);

  double outward_nml[2];
  outward_nml[0] = seg_head[1] - seg_tail[1];
  outward_nml[1] = seg_tail[0] - seg_head[0];
  EXPECT_DOUBLE_EQ(1, outward_nml[0]);
  EXPECT_DOUBLE_EQ(4, outward_nml[1]);

  normalize(outward_nml);
  EXPECT_DOUBLE_EQ(1, len_2d(outward_nml));

  double seg_dir[2];
  subtract_2d(seg_head, seg_tail, seg_dir);
  normalize(seg_dir);
  EXPECT_DOUBLE_EQ(1, len_2d(seg_dir));
  EXPECT_DOUBLE_EQ(0, dot_prod_2d(outward_nml, seg_dir));

  double point_above[2] = {2, .6};
  double point_above_rel[2];
  subtract_2d(point_above, seg_midpoint, point_above_rel);
  EXPECT_DOUBLE_EQ(0, point_above_rel[0]);
  EXPECT_DOUBLE_EQ(.1, point_above_rel[1]);

  EXPECT_TRUE(!is_zero(0.1));
  EXPECT_TRUE(!is_zero(0.1, 1e-100));
  EXPECT_TRUE(is_zero(0.1, 1e20));
}

TEST(morkon_planar_geom, classify_point_2d)
{
  using namespace morkon_exp;

  double seg_head[2] = {0, 1};
  double seg_tail[2] = {4, 0};

  double line_witness[2], outward_normal[2];
  compute_supporting_line(seg_tail, seg_head, line_witness, outward_normal);

  double seg_dir[2];
  subtract_2d(seg_head, seg_tail, seg_dir);
  normalize(seg_dir);
  EXPECT_DOUBLE_EQ(1, len_2d(seg_dir));
  EXPECT_DOUBLE_EQ(1, len_2d(outward_normal));
  EXPECT_DOUBLE_EQ(0, dot_prod_2d(outward_normal, seg_dir));

  double point_above[2] = {2, .6};
  double point_below[2] = {2, .4};
  double point_on[2]    = {2, .5};

  EXPECT_EQ(POINT_ABOVE, classify_point_2d(point_above, line_witness, outward_normal));
  EXPECT_EQ(POINT_BELOW, classify_point_2d(point_below, line_witness, outward_normal));
  EXPECT_EQ(POINT_ON, classify_point_2d(point_on, line_witness, outward_normal));

  point_on[1] += 5 * std::numeric_limits<double>::epsilon();
  EXPECT_EQ(POINT_ON, classify_point_2d(point_on, line_witness, outward_normal));

  point_on[1] += 500 * std::numeric_limits<double>::epsilon();
  EXPECT_EQ(POINT_ABOVE, classify_point_2d(point_on, line_witness, outward_normal));
}

template <typename ScalarType, unsigned DIM>
void complain_not_epsilon_equal(int vec_len,
                                const ScalarType vecListA[][DIM], const ScalarType vecListB[][DIM])
{
  std::cerr << "expected";
  streamit(std::cerr, vec_len, vecListA) << " and ";
  streamit(std::cerr, vec_len, vecListB) << " to be epsilon-equal" << std::endl;
}

TEST(morkon_planar_geom, clip_convex_polygon_with_halfspace)
{
  using namespace morkon_exp;

  const double triangle_0[3][2] = {{1,1}, {5,1}, {3,5}};
  double gon_out_buff[10][2];

  const double nml_0[2]      = {0,-1};
  const double witness_0[2] = {3, 0.5};

  int num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_0, nml_0, gon_out_buff);
  EXPECT_EQ(3, num_gon_out_verts);
  bool ok = are_lists_epsilon_equal(3, triangle_0, gon_out_buff);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(3, triangle_0, gon_out_buff);

  const double witness_1[2] = {3, 1 - 1e-15};
  num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_1, nml_0, gon_out_buff);
  EXPECT_EQ(3, num_gon_out_verts);
  ok = are_lists_epsilon_equal(3, triangle_0, gon_out_buff);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(3, triangle_0, gon_out_buff);

  const double witness_2[2] = {3, 1 + 1e-15};
  num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_2, nml_0, gon_out_buff);
  EXPECT_EQ(3, num_gon_out_verts);
  ok = are_lists_epsilon_equal(3, triangle_0, gon_out_buff);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(3, triangle_0, gon_out_buff);

  const double witness_3[2] = {3, 2};
  const double answer_3[3][2] =  { {4.5,2}, {3,5}, {1.5,2} };
  num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_3, nml_0, gon_out_buff);
  EXPECT_EQ(3, num_gon_out_verts);
  ok = are_lists_epsilon_equal(3, answer_3, gon_out_buff);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(3, answer_3, gon_out_buff);

  const double witness_4[2] = {3, 1 + 1e-5};
  const double answer_4[3][2] =  { {4.99999,1.00001}, {3,5}, {1,1.00001} };
  num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_4, nml_0, gon_out_buff);
  EXPECT_EQ(3, num_gon_out_verts);
  ok = are_lists_epsilon_equal(3, answer_4, gon_out_buff, 0.00001);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(3, answer_4, gon_out_buff);

  const double nml_1[2]       = {0,1};
  const double witness_5[2]   = {1,4};
  const double answer_5[4][2] = {{1,1}, {5,1}, {3.5,4}, {2.5,4}};
  num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_5, nml_1, gon_out_buff);
  EXPECT_EQ(4, num_gon_out_verts);
  ok = are_lists_epsilon_equal(4, answer_5, gon_out_buff);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(4, answer_5, gon_out_buff);

  const double dir_3[2] = {2,-1};
  double nml_3[2];
  normalize(dir_3, nml_3);
  const double witness_6[2]   = {3,1};
  const double answer_6[4][2] = {{1,1}, {3,1}, {4,3}, {3,5}};
  num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_6, nml_3, gon_out_buff);
  EXPECT_EQ(4, num_gon_out_verts);
  ok = are_lists_epsilon_equal(4, answer_6, gon_out_buff);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(4, answer_6, gon_out_buff);

  const double dir_4[2] = {-2,-1};
  double nml_4[2];
  normalize(dir_4, nml_4);
  const double witness_7[2]   = {3,1};
  const double answer_7[4][2] = {{3,1}, {5,1}, {3,5}, {2,3}};
  num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_7, nml_4, gon_out_buff);
  EXPECT_EQ(4, num_gon_out_verts);
  ok = are_lists_epsilon_equal(4, answer_7, gon_out_buff);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(4, answer_7, gon_out_buff);

  const double dir_5[2] = {2,-3};
  double nml_5[2];
  normalize(dir_5, nml_5);
  const double witness_8[2]   = {1,1};
  const double answer_8[3][2] = {{1,1}, {4,3}, {3,5}};
  num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_8, nml_5, gon_out_buff);
  EXPECT_EQ(3, num_gon_out_verts);
  ok = are_lists_epsilon_equal(3, answer_8, gon_out_buff);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(3, answer_8, gon_out_buff);

  const double nml_6[2]       = {1, 0};
  const double witness_9[2]   = {3, 3};
  const double answer_9[3][2] = {{1,1}, {3,1}, {3,5}};
  num_gon_out_verts = clip_convex_polygon(3, triangle_0, witness_9, nml_6, gon_out_buff);
  EXPECT_EQ(3, num_gon_out_verts);
  ok = are_lists_epsilon_equal(3, answer_9, gon_out_buff);
  EXPECT_TRUE(ok);
  if (!ok)
    complain_not_epsilon_equal(3, answer_9, gon_out_buff);

}


