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

#ifndef STK_SEARCH_BOX_HPP
#define STK_SEARCH_BOX_HPP

#include <Kokkos_Core.hpp>
#include <stk_search/Point.hpp>

namespace stk { namespace search {

template <typename T>
class Box
{
public:
  typedef T value_type;
  typedef Point<value_type> point_type;
  static const int Dim = 3;

  static KOKKOS_FUNCTION constexpr value_type max() { return Kokkos::Experimental::finite_max_v<T>;}
  static KOKKOS_FUNCTION constexpr value_type min() {
    // finite_min_v<T> returns the most negative real value (equivalent to numeric_limits<T>::lowest).
    // it is the 'lowest' value that we want here.  
    return Kokkos::Experimental::finite_min_v<T>;
  }

  KOKKOS_FUNCTION Box()
  : m_min_corner(max(), max(), max()), m_max_corner(min(), min(), min()) {}

  KOKKOS_FUNCTION Box( point_type const& minCorner,
                       point_type const& maxCorner)
    : m_min_corner(minCorner)
    , m_max_corner(maxCorner)
  {
  }

  KOKKOS_FUNCTION Box( const T x_min, const T y_min, const T z_min,
       const T x_max, const T y_max, const T z_max) :
    m_min_corner(x_min, y_min, z_min),
    m_max_corner(x_max, y_max, z_max) 
  {}

  KOKKOS_FUNCTION point_type const& min_corner() const { return m_min_corner; }
  KOKKOS_FUNCTION point_type      & min_corner()       { return m_min_corner; }
  KOKKOS_FUNCTION point_type const& max_corner() const { return m_max_corner; }
  KOKKOS_FUNCTION point_type      & max_corner()       { return m_max_corner; }

  KOKKOS_FUNCTION value_type get_x_min() const { return m_min_corner[0]; }
  KOKKOS_FUNCTION value_type get_y_min() const { return m_min_corner[1]; }
  KOKKOS_FUNCTION value_type get_z_min() const { return m_min_corner[2]; }
  KOKKOS_FUNCTION value_type get_x_max() const { return m_max_corner[0]; }
  KOKKOS_FUNCTION value_type get_y_max() const { return m_max_corner[1]; }
  KOKKOS_FUNCTION value_type get_z_max() const { return m_max_corner[2]; }

  KOKKOS_FORCEINLINE_FUNCTION
  float get_expanded_x_min() const { return Kokkos::nextafter(static_cast<float>(m_min_corner[0]),
                                                              Kokkos::Experimental::finite_min_v<float>); }

  KOKKOS_FORCEINLINE_FUNCTION
  float get_expanded_y_min() const { return Kokkos::nextafter(static_cast<float>(m_min_corner[1]),
                                                              Kokkos::Experimental::finite_min_v<float>); }

  KOKKOS_FORCEINLINE_FUNCTION
  float get_expanded_z_min() const { return Kokkos::nextafter(static_cast<float>(m_min_corner[2]),  
                                                              Kokkos::Experimental::finite_min_v<float>);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  float get_expanded_x_max() const { return Kokkos::nextafter(static_cast<float>(m_max_corner[0]),
                                                              Kokkos::Experimental::finite_max_v<float>); }

  KOKKOS_FORCEINLINE_FUNCTION
  float get_expanded_y_max() const { return Kokkos::nextafter(static_cast<float>(m_max_corner[1]),
                                                              Kokkos::Experimental::finite_max_v<float>); }

  KOKKOS_FORCEINLINE_FUNCTION
  float get_expanded_z_max() const { return Kokkos::nextafter(static_cast<float>(m_max_corner[2]),
                                                              Kokkos::Experimental::finite_max_v<float>); }

  KOKKOS_FUNCTION void set_min_corner(point_type const& x_min_corner) { m_min_corner = x_min_corner; }
  KOKKOS_FUNCTION void set_max_corner(point_type const& x_max_corner) { m_max_corner = x_max_corner; }

  KOKKOS_FUNCTION void set_box(const value_type x1, const value_type y1, const value_type z1,
               const value_type x2, const value_type y2, const value_type z2)
  {
    set_min_corner(point_type(x1, y1, z1));
    set_max_corner(point_type(x2, y2, z2));
  }

  KOKKOS_FUNCTION bool operator==(Box<value_type> const& b) const
  { return m_min_corner == b.m_min_corner && m_max_corner == b.m_max_corner; }

  KOKKOS_FUNCTION bool operator!=(Box<value_type> const& b) const
  { return !(*this == b); }

  KOKKOS_FUNCTION void set_box(const Box& b)
  {
    set_min_corner(b.min_corner());
    set_max_corner(b.max_corner());
  }

private:
  point_type m_min_corner;
  point_type m_max_corner;
};

template <typename T>
std::ostream& operator<<(std::ostream & out, Box<T> const& b)
{
  out << "{" << b.min_corner() << "->" << b.max_corner() << "}";
  return out;
}

template<typename T>
std::istream& operator>>(std::istream& in, Box<T>& b) {
  char c;
  in >> c >> b.min_corner() >> c >> c >> b.max_corner() >> c;
  return in;
}

}} //namespace stk::search

#endif //STK_SEARCH_BOX_HPP

