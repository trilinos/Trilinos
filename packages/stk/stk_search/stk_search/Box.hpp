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

#ifndef STK_SEARCH_BOX_HPP
#define STK_SEARCH_BOX_HPP

#include <stk_search/Point.hpp>

namespace stk { namespace search {

template <typename T>
class Box
{
public:
  typedef T value_type;
  typedef Point<value_type> point_type;
  static const int Dim = 3;

  Box( point_type const& x_min_corner = point_type(std::numeric_limits<T>::max(),std::numeric_limits<T>::max(),std::numeric_limits<T>::max()), 
       point_type const& x_max_corner = point_type(std::numeric_limits<T>::lowest(),std::numeric_limits<T>::lowest(),std::numeric_limits<T>::lowest()))
    : m_min_corner(x_min_corner)
    , m_max_corner(x_max_corner)
  {}

  point_type const& min_corner() const { return m_min_corner; }
  point_type      & min_corner()       { return m_min_corner; }
  point_type const& max_corner() const { return m_max_corner; }
  point_type      & max_corner()       { return m_max_corner; }

  value_type get_x_min() const { return m_min_corner[0]; }
  value_type get_y_min() const { return m_min_corner[1]; }
  value_type get_z_min() const { return m_min_corner[2]; }
  value_type get_x_max() const { return m_max_corner[0]; }
  value_type get_y_max() const { return m_max_corner[1]; }
  value_type get_z_max() const { return m_max_corner[2]; }

  void set_min_corner(point_type const& x_min_corner) { m_min_corner = x_min_corner; }
  void set_max_corner(point_type const& x_max_corner) { m_max_corner = x_max_corner; }

  void set_box(const value_type x1, const value_type y1, const value_type z1,
               const value_type x2, const value_type y2, const value_type z2)
  {
    set_min_corner(point_type(x1, y1, z1));
    set_max_corner(point_type(x2, y2, z2));
  }

  bool operator==(Box<value_type> const& b) const
  { return m_min_corner == b.m_min_corner && m_max_corner == b.m_max_corner; }

  bool operator!=(Box<value_type> const& b) const
  { return !(*this == b); }

  //  NKC, keeping track of things may need to put here to support usage of this box in GTK search
  /*
  Box<value_type>& operator = (const Box<value_type>& box1) {
    m_min_corner[0] = box1.m_min_corner[0];
    m_min_corner[1] = box1.m_min_corner[1];
    m_min_corner[2] = box1.m_min_corner[2];
    m_max_corner[0] = box1.m_max_corner[0];
    m_max_corner[1] = box1.m_max_corner[1];
    m_max_corner[2] = box1.m_max_corner[2];
    return *this;
  }
  */

  friend std::ostream& operator<<(std::ostream & out, Box<value_type> const& b)
  {
    out << "{" << b.min_corner() << "->" << b.max_corner() << "}";
    return out;
  }

private:
  point_type m_min_corner;
  point_type m_max_corner;
};

}} //namespace stk::search


namespace boost {
namespace geometry {
namespace traits {

// traits for stk::search::Box<T>
template <typename T> struct tag< stk::search::Box<T> > { typedef box_tag type; };
template <typename T> struct point_type< stk::search::Box<T> > { typedef stk::search::Point<T> type; };

template <typename T, size_t Index>
struct indexed_access< stk::search::Box<T>, min_corner, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::Box<T> const& s) { return s.min_corner()[Index]; }
};




template <typename T, size_t Index>
struct indexed_access< stk::search::Box<T>, max_corner, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::Box<T> const& s) { return s.max_corner()[Index]; }
};


}}} // namespace boost::geometry::traits

#endif //STK_SEARCH_BOX_HPP

