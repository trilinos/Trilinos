/*
 * Copyright (c) 2013, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/*
 * gtkTraitsForSearch.h
 *
 *  Created on: Dec 18, 2013
 *      Author: mkbhard
 */

#ifndef GTKTRAITSFORSEARCH_H_
#define GTKTRAITSFORSEARCH_H_

#include <Geom_AxisAlignedBB.h>

namespace boost {
namespace geometry {
namespace traits {

// traits for ::gtk::AxisAlignedBB
template <> struct tag< ::gtk::AxisAlignedBB > { typedef box_tag type; };
template <> struct point_type< ::gtk::AxisAlignedBB > { typedef stk::search::Point<Real> type; };

template <>
struct indexed_access< ::gtk::AxisAlignedBB, min_corner, 0 >
{  static inline Real get( ::gtk::AxisAlignedBB const& s) { return s.get_x_min(); } };

template <>
struct indexed_access< ::gtk::AxisAlignedBB, min_corner, 1 >
{  static inline Real get( ::gtk::AxisAlignedBB const& s) { return s.get_y_min(); } };

template <>
struct indexed_access< ::gtk::AxisAlignedBB, min_corner, 2 >
{  static inline Real get( ::gtk::AxisAlignedBB const& s) { return s.get_z_min(); } };

template <>
struct indexed_access< ::gtk::AxisAlignedBB, max_corner, 0 >
{  static inline Real get( ::gtk::AxisAlignedBB const& s) { return s.get_x_max(); } };

template <>
struct indexed_access< ::gtk::AxisAlignedBB, max_corner, 1 >
{  static inline Real get( ::gtk::AxisAlignedBB const& s) { return s.get_y_max(); } };

template <>
struct indexed_access< ::gtk::AxisAlignedBB, max_corner, 2 >
{  static inline Real get( ::gtk::AxisAlignedBB const& s) { return s.get_z_max(); } };

}}} // namespace boost::geometry::traits

namespace stk {
namespace search {

inline bool intersects(const gtk::AxisAlignedBB& a, const gtk::AxisAlignedBB& b)
{
  return a.overlap(b);
}

// intersects: Sphere,Box
template <typename T>
inline bool intersects(Sphere<T> const& a, gtk::AxisAlignedBB const& b)
{
  Point<T> const& ac = a.center();
  Point<T> const& bmin = Point<T>(b.get_x_min(), b.get_y_min(), b.get_z_min());
  Point<T> const& bmax = Point<T>(b.get_x_max(), b.get_y_max(), b.get_z_max());

  const T r2 = a.radius() * a.radius();

  // check that the nearest point in the bounding box is within the sphere
  T dmin = 0;
  for( int i = 0; i < 3; ++i ) {
    if( ac[i] < bmin[i] ) dmin += (ac[i]-bmin[i])*(ac[i]-bmin[i]);
    else if( ac[i] > bmax[i] ) dmin += (ac[i]-bmax[i])*(ac[i]-bmax[i]);
  }
  return dmin <= r2;
}

// intersects: Box,Sphere
template <typename T>
inline bool intersects(gtk::AxisAlignedBB const& a, Sphere<T> const& b)
{ return intersects(b,a); }

}} // namespace stk::search

#endif /* GTKTRAITSFORSEARCH_H_ */
