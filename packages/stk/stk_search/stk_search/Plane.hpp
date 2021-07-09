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

#ifndef STK_SEARCH_PLANE_HPP
#define STK_SEARCH_PLANE_HPP

#include <cmath>
#include <Kokkos_ArithTraits.hpp>
#include <stk_search/Point.hpp>

namespace stk { namespace search {

template <typename T>
class Plane
{
public:
  typedef T value_type;
  typedef Point<value_type> point_type;
  static const int Dim = 3;

  static KOKKOS_FUNCTION constexpr value_type max() { return Kokkos::Details::ArithTraits<T>::max() ;}

  KOKKOS_FUNCTION point_type cross(const point_type& a, const point_type& b) const{
    return point_type(a[1]*b[2] - a[2]*b[1],
                      a[2]*b[0] - a[0]*b[2],
                      a[0]*b[1] - a[1]*b[0]);
  }

  KOKKOS_FUNCTION void normalize(point_type& a) const{
    const value_type vec_len = Kokkos::Details::ArithTraits<T>::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    const value_type denom = vec_len == 0.0 ? 1.0 : vec_len;
    const value_type vec_len_inv = 1.0 / denom;
    a[0] *= vec_len_inv;
    a[1] *= vec_len_inv;
    a[2] *= vec_len_inv;
  }

  KOKKOS_FUNCTION value_type dot(const point_type& a, const point_type& b) const {
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
  }

  KOKKOS_FUNCTION value_type SignedDistance(const point_type& P) const { return dot(normal, P) - constant; }

  KOKKOS_FUNCTION
  int WhichSide(const point_type& P) const {
    value_type Distance = SignedDistance(P);
    if (Distance < 0.0) {
      return -1;
    } else if (Distance > 0.0) {
      return 1;
    } else {  // Distance == 0.0
      return 0;
    }
  }

  KOKKOS_FUNCTION Plane() :
    normal(max()),
    constant(max())
  {
  }


  KOKKOS_FUNCTION Plane( const point_type normal_, const value_type constant_) :
    normal(normal_),
    constant(constant_)
  {
  }

  KOKKOS_FUNCTION Plane( const point_type pointOnPlane, const point_type normal_) :
    normal(normal_),
    constant()
  {
    constant = dot(normal, pointOnPlane);
  }

  KOKKOS_FUNCTION Plane( const point_type P0, const point_type P1, const point_type P2) :
    normal(),
    constant()
  {
    point_type Edge1 = P1 - P0;
    point_type Edge2 = P2 - P0;
    normal = cross(Edge1, Edge2);
    normalize(normal);
    constant = dot(normal, P0);
  }

  KOKKOS_FUNCTION bool operator==(Plane<value_type> const& b) const
  { return normal == b.normal && constant == b.constant; }

  KOKKOS_FUNCTION bool operator!=(Plane<value_type> const& b) const
  { return !(*this == b); }

  KOKKOS_DEFAULTED_FUNCTION ~Plane() = default;

private:
  point_type normal;
  value_type constant;
};

}} //namespace stk::search

#endif //STK_SEARCH_Plane_HPP

