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

#ifndef STK_SEARCH_SPHERE_HPP
#define STK_SEARCH_SPHERE_HPP

#include <stk_search/Point.hpp>
#include <Kokkos_Core.hpp>

namespace stk { namespace search {

template <typename T>
class Sphere
{
public:
  typedef T value_type;
  typedef Point<value_type> point_type;
  static const int Dim = 3;

  KOKKOS_FORCEINLINE_FUNCTION Sphere( point_type const& x_center = point_type(), value_type const& x_radius = static_cast<value_type>(-1))
    : m_center(x_center)
    , m_radius(x_radius)
  {}

  KOKKOS_FORCEINLINE_FUNCTION void set_box(const Sphere &val) {
    m_center = val.m_center;
    m_radius = val.m_radius;
  }

  KOKKOS_FORCEINLINE_FUNCTION point_type const& center() const { return m_center; }
  KOKKOS_FORCEINLINE_FUNCTION point_type      & center()       { return m_center; }

  KOKKOS_FORCEINLINE_FUNCTION value_type const& radius() const { return m_radius; }
  KOKKOS_FORCEINLINE_FUNCTION value_type      & radius()       { return m_radius; }

  void set_center(point_type const& c) { m_center = c; }
  void set_radius(value_type const& r) { m_radius = r; }

  bool operator==(Sphere<value_type> const& s) const
  { return m_radius == s.m_radius && m_center == s.m_center; }

  bool operator!=(Sphere<value_type> const& s) const
  { return !(*this == s); }

  KOKKOS_FORCEINLINE_FUNCTION value_type get_x_min() const { return m_center[0] - m_radius; }
  KOKKOS_FORCEINLINE_FUNCTION value_type get_y_min() const { return m_center[1] - m_radius; }
  KOKKOS_FORCEINLINE_FUNCTION value_type get_z_min() const { return m_center[2] - m_radius; }
  KOKKOS_FORCEINLINE_FUNCTION value_type get_x_max() const { return m_center[0] + m_radius; }
  KOKKOS_FORCEINLINE_FUNCTION value_type get_y_max() const { return m_center[1] + m_radius; }
  KOKKOS_FORCEINLINE_FUNCTION value_type get_z_max() const { return m_center[2] + m_radius; }

  KOKKOS_FORCEINLINE_FUNCTION
  float get_expanded_radius() const { 
    return Kokkos::nextafter(static_cast<float>(m_radius), Kokkos::Experimental::finite_max_v<float>); 
  }

  KOKKOS_FORCEINLINE_FUNCTION float get_expanded_x_min() const { return static_cast<float>(m_center[0]) - this->get_expanded_radius(); }
  KOKKOS_FORCEINLINE_FUNCTION float get_expanded_y_min() const { return static_cast<float>(m_center[1]) - this->get_expanded_radius(); }
  KOKKOS_FORCEINLINE_FUNCTION float get_expanded_z_min() const { return static_cast<float>(m_center[2]) - this->get_expanded_radius(); }
  KOKKOS_FORCEINLINE_FUNCTION float get_expanded_x_max() const { return static_cast<float>(m_center[0]) + this->get_expanded_radius(); }
  KOKKOS_FORCEINLINE_FUNCTION float get_expanded_y_max() const { return static_cast<float>(m_center[1]) + this->get_expanded_radius(); }
  KOKKOS_FORCEINLINE_FUNCTION float get_expanded_z_max() const { return static_cast<float>(m_center[2]) + this->get_expanded_radius(); }

  friend std::ostream& operator<<(std::ostream & out, Sphere<value_type> const& s)
  {
    out << "{" << s.center() << ":" << s.radius() << "}";
    return out;
  }

private:
  point_type m_center;
  value_type m_radius;
};

}} //namespace stk::search



#endif //STK_SEARCH_SPHERE_HPP
