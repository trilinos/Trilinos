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

#ifndef STK_STK_SEARCH_UTIL_PARMETRICDISTANCE_HPP_
#define STK_STK_SEARCH_UTIL_PARMETRICDISTANCE_HPP_

#include <stk_util/stk_config.h>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_topology/topology.hpp>
#include <array>
#include <algorithm>
#include <cmath>

namespace stk {
namespace search {

namespace impl {

inline
double param_dist_util_wedge(const double X, const double Y)
{
  const double dist0 = -3 * X;
  const double dist1 = -3 * Y;
  const double dist2 = 3 * (X + Y);
  const double dist = std::max(std::max(dist0, dist1), dist2);
  return dist;
}

inline
double vector_norm(const double* vect, int N)
{
  double norm_sq = 0.0;
  for (int i = 0; i < N; ++i) {
    norm_sq += vect[i] * vect[i];
  }
  return norm_sq;
}

} // namespace impl

template<typename ArrayType>
double parametric_distance_line2(const ArrayType& x)
{
  STK_ThrowAssertMsg(x.size() == 2, "parametric_distance_line2 requires 2 coords");
  constexpr double ELEMENT_THICKNESS = 0.01;
  double dist = std::fabs(x[0]);
  if(ELEMENT_THICKNESS < x[1] && dist < 1 + x[1]) dist = 1 + x[1];
  return dist;
}

template<typename ArrayType>
double parametric_distance_linear_quad(const ArrayType& x)
{
  STK_ThrowAssertMsg(x.size() == 3, "parametric_distance_linear_quad requires 3 coords");
  constexpr double ELEM_THICK = 0.01;
  std::array<double, 3> y = { { std::fabs(x[0]), std::fabs(x[1]), std::fabs(x[2]) } };
  double dist = y[0];
  if(dist < y[1]) dist = y[1];
  if(ELEM_THICK < y[2] && dist < 1 + y[2]) dist = 1 + y[2];
  return dist;
}

template<typename ArrayType>
double parametric_distance_quad_2d(const ArrayType& x) {

  STK_ThrowAssertMsg(x.size() == 2, "parametric_distance_quad_2d requires 2 coords");

  ArrayType y = {{ std::fabs(x[0]), std::fabs(x[1]) }};
  const double d = *std::max_element(y.begin(), y.end());
  return d;  

}

inline
double compute_linear_quad_cpp_distance(const double* normal,
                                        const double* solcur,
                                        const double* deltasol,
                                        const int& npar_coord,
                                        double* par_coor)
{
  // Rescale the distance vector by the length of the (non-unit) normal vector,
  // which was used above in the NR iteration.
  const double area = std::sqrt(impl::vector_norm(normal, 3));
  const double length = std::sqrt(area);
  const double par_coor_2 = (solcur[2] + deltasol[2]) * length;
  if(npar_coord == 3) par_coor[2] = par_coor_2;

  std::array<double, 3> xtmp = { { par_coor[0], par_coor[1], par_coor_2 } };
  return parametric_distance_linear_quad(xtmp);
}

template<typename ArrayType>
double parametric_distance_hex(const ArrayType& x)
{
  STK_ThrowAssertMsg(x.size() == 3, "parametric_distance_hex requires 3 coords");
  std::array<double, 3> y = {{ std::fabs(x[0]), std::fabs(x[1]), std::fabs(x[2]) }};

  double dist = 0.0;
  for(int i = 0; i < 3; ++i) {
    if(dist < y[i]) {
      dist = y[i];
    }
  }
  return dist;
}

template<typename ArrayType>
double parametric_distance_pyramid(const ArrayType& x)
{
  STK_ThrowAssertMsg(x.size() == 3, "parametric_distance_pyramid requires 3 coords");
  const double X = x[0];
  const double Y = x[1];
  const double Z = x[2] - 1. / 3.;
  const double dist0 = (3. / 2.) * (Z + std::max(std::fabs(X), std::fabs(Y)));
  const double dist1 = -3 * Z;
  const double dist = std::max(dist0, dist1);
  return dist;
}

template<typename ArrayType>
double parametric_distance_tri_2d(const ArrayType& x) {

  STK_ThrowAssertMsg(x.size() == 2, "parametric_distance_tri_2d requires 2 coords");

  const double X = x[0] - 1. / 3.;
  const double Y = x[1] - 1. / 3.;
  const double dist0 = -3 * X;
  const double dist1 = -3 * Y;
  const double dist2 = 3 * (X + Y);
  const double dist = std::max(std::max(dist0, dist1), dist2);
  return dist;
}

template<typename ArrayType>
double parametric_distance_tet(const ArrayType& x)
{
  STK_ThrowAssertMsg(x.size() == 3, "parametric_distance_tet requires 3 coords");
  const double X = x[0] - 1. / 4.;
  const double Y = x[1] - 1. / 4.;
  const double Z = x[2] - 1. / 4.;
  const double dist0 = -4 * X;
  const double dist1 = -4 * Y;
  const double dist2 = -4 * Z;
  const double dist3 = 4 * (X + Y + Z);
  const double dist = std::max(std::max(dist0, dist1), std::max(dist2, dist3));
  return dist;
}


template<typename ArrayType>
double parametric_distance_wedge(const ArrayType& x)
{
  STK_ThrowAssertMsg(x.size() == 3, "parametric_distance_wedge requires 3 coords");
  const double X = x[0] - 1. / 3.;
  const double Y = x[1] - 1. / 3.;
  const double Z = x[2];
  const double dist_t = impl::param_dist_util_wedge(X, Y);
  const double dist_z = std::fabs(Z);
  const double dist = std::max(dist_z, dist_t);
  return dist;
}

} // namespace search
} // namespace stk

#endif /* STK_STK_SEARCH_UTIL_PARMETRICDISTANCE_HPP_ */
