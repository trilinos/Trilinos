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

#ifndef STK_SEARCH_POINT_HPP
#define STK_SEARCH_POINT_HPP

#include <stk_util/environment/ReportHandler.hpp>
#include <boost/geometry/geometry.hpp>
#include <iosfwd>

namespace stk { namespace search {

template <typename T>
class Point
{
public:
  static const unsigned Dim = 3;
  typedef T value_type;

  Point(value_type x = value_type(), value_type y = value_type(), value_type z = value_type())
    : m_value()
  {
    m_value[0] = x; m_value[1] = y; m_value[2] = z;
  }

  value_type const& operator[](size_t index) const
  {
    ThrowAssert(index < Dim);
    return m_value[index];
  }

  value_type & operator[](size_t index)
  {
    ThrowAssert(index < Dim);
    return m_value[index];
  }

  bool operator==(Point<value_type> const& p) const
  {
    return  m_value[0] == p.m_value[0]
         && m_value[1] == p.m_value[1]
         && m_value[2] == p.m_value[2];
  }

  bool operator!=(Point<value_type> const& p) const
  { return !(*this == p); }

  value_type get_x_min() const { return m_value[0]; }
  value_type get_y_min() const { return m_value[1]; }
  value_type get_z_min() const { return m_value[2]; }
  value_type get_x_max() const { return m_value[0]; }
  value_type get_y_max() const { return m_value[1]; }
  value_type get_z_max() const { return m_value[2]; }


  friend std::ostream& operator<<(std::ostream & out, Point<value_type> const& p)
  {
    out << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
    return out;
  }

private:
  value_type m_value[Dim];
};

}} // stk::search

namespace boost { namespace geometry { namespace traits {

// traits for stk::search::Point<T>
template <typename T> struct tag< stk::search::Point<T> > { typedef point_tag type; };
template <typename T> struct coordinate_type< stk::search::Point<T> > { typedef T type; };
template <typename T> struct coordinate_system< stk::search::Point<T> > { typedef cs::cartesian type; };
template <typename T> struct dimension< stk::search::Point<T> > : public boost::mpl::int_<3> {};

template <typename T, size_t Index>
struct access< stk::search::Point<T>, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::Point<T> const& p) { return p[Index]; }
  static inline void set( stk::search::Point<T> const& p, T const& v) { p[Index] = v; }
};

}}} // namespace boost::geometry::traits

#endif //STK_SEARCH_POINT_HPP
