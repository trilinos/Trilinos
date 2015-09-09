/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
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

#ifndef MORKON_EXP_PLANAR_GEOM_H
#define MORKON_EXP_PLANAR_GEOM_H

#include <cmath>
#include <mrk_data_types.hpp>

namespace morkon_exp {

enum PointLineRelationship {POINT_ON = 0, POINT_ABOVE, POINT_BELOW, NUM_POINTLINE_RELS};

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
bool is_zero(const ScalarType val, const ScalarType epsilon_coeff_arg = 10.0)
{
  ScalarType epsilon_coefficient = (epsilon_coeff_arg >= 10.0 ? epsilon_coeff_arg : 10.0);
  return std::abs(val) <= epsilon_coefficient * std::numeric_limits<ScalarType>::epsilon();
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
void scale_2d(const ScalarType factor, ScalarType point[2])
{
  point[0] = factor * point[0];
  point[1] = factor * point[1];
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
void scale_2d(const ScalarType factor, const ScalarType pointIn[2], ScalarType pointOut[2])
{
  pointOut[0] = factor * pointIn[0];
  pointOut[1] = factor * pointIn[1];
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
ScalarType dot_prod_2d(const ScalarType vecA[2], const ScalarType vecB[2])
{
  return vecA[0] * vecB[0] + vecA[1] * vecB[1];
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
ScalarType len_2d(const ScalarType vec[2])
{
  return sqrt(dot_prod_2d(vec, vec));
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
void normalize(ScalarType vec[2])
{
  ScalarType inv_len = 1.0 / len_2d(vec);
  vec[0] *= inv_len;
  vec[1] *= inv_len;
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
void normalize(const ScalarType vec_in[2], ScalarType vec_out[2])
{
  ScalarType inv_len = 1.0 / len_2d(vec_in);
  vec_out[0] = inv_len * vec_in[0];
  vec_out[1] = inv_len * vec_in[1];
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
void add_2d(const ScalarType pointA[2], const ScalarType pointB[2], ScalarType result[2])
{
  result[0] = pointA[0] + pointB[0];
  result[1] = pointA[1] + pointB[1];
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
void subtract_2d(const ScalarType pointA[2], const ScalarType pointB[2], ScalarType result[2])
{
  result[0] = pointA[0] - pointB[0];
  result[1] = pointA[1] - pointB[1];
}


template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
void compute_supporting_line(const ScalarType segment_tail[2], const ScalarType segment_head[2],
                             ScalarType line_witness[2], ScalarType outward_normal[2])
{
  add_2d(segment_tail, segment_tail, line_witness);
  scale_2d(0.5, line_witness);
  outward_normal[0] = segment_head[1] - segment_tail[1];
  outward_normal[1] = segment_tail[0] - segment_head[0];

  normalize(outward_normal);
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
PointLineRelationship
classify_point_2d(const ScalarType point[2],
                  const ScalarType line_witness[2], const ScalarType outward_normal[2],
                  ScalarType epsilon_arg = 1e-10)
{
  ScalarType local_point[2];
  subtract_2d(point, line_witness, local_point);
  ScalarType val = dot_prod_2d(outward_normal, local_point);

  const ScalarType min_epsilon = 100.0 * len_2d(local_point);
  ScalarType epsilon = (epsilon_arg > min_epsilon ? epsilon_arg : min_epsilon);

  if (is_zero(val, epsilon))
  {
    return POINT_ON;
  }
  else
  {
    PointLineRelationship retval = (val < 0 ? POINT_BELOW : POINT_ABOVE);
    return retval;
  }
}


template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
void find_crossing_point_2d(const ScalarType pointA[2], const ScalarType pointB[2],
                            const ScalarType line_witness[2], const ScalarType outward_normal[2],
                            ScalarType crossing_point[2])
{
  ScalarType A_B[2], W_B[2];
  subtract_2d(pointA, pointB, A_B);
  subtract_2d(line_witness, pointB, W_B);

  ScalarType A_proj  = dot_prod_2d(A_B, outward_normal);
  ScalarType W_proj = dot_prod_2d(W_B, outward_normal);
  ScalarType frac = W_proj / A_proj;

  ScalarType scaled_A_B[2];
  scale_2d(frac, A_B, scaled_A_B);
  add_2d(pointB, scaled_A_B, crossing_point);
}


// Assumes that gon_in is convex!
template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
int clip_convex_polygon(int num_pts, const ScalarType gon_in[][2],
                        const ScalarType line_witness[2], const ScalarType outward_normal[2],
                        ScalarType gon_out[][2] // Assume big enough!
                                      )
{
  int out_idx = 0;
  int curr_tail = 0;

  PointLineRelationship curr_tail_rel = classify_point_2d(gon_in[curr_tail], line_witness, outward_normal);

  for (; curr_tail < num_pts; ++curr_tail)
  {
    if (curr_tail_rel != POINT_ABOVE)
    {
      gon_out[out_idx][0] = gon_in[curr_tail][0];
      gon_out[out_idx][1] = gon_in[curr_tail][1];
      ++out_idx;
    }
    int curr_head = (curr_tail + 1) % num_pts;
    PointLineRelationship curr_head_rel = classify_point_2d(gon_in[curr_head], line_witness, outward_normal);

    if ((curr_tail_rel == POINT_ABOVE && curr_head_rel == POINT_BELOW)
        || (curr_tail_rel == POINT_BELOW && curr_head_rel == POINT_ABOVE))
    {
      ScalarType crossing_pt[2];
      find_crossing_point_2d(gon_in[curr_tail], gon_in[curr_head], line_witness, outward_normal,
                             crossing_pt);
      gon_out[out_idx][0] = crossing_pt[0];
      gon_out[out_idx][1] = crossing_pt[1];
      ++out_idx;
    }

    curr_tail_rel = curr_head_rel;
  }

  return out_idx;
}


template <typename ScalarType, int MAX_VERTS_IN = 3>
KOKKOS_INLINE_FUNCTION
void clip_convex_polygon(int num_subject_gon_verts, ScalarType subject_gon[][2],
                         int num_clip_gon_verts, ScalarType clip_gon[][2],
                         ScalarType gon_out[][2])
{
  const int BuffGonCapacity= 2 * MAX_VERTS_IN;
  ScalarType scratch_gons[2][BuffGonCapacity][2];

  // YOU ARE HERE!
}


}

#endif
