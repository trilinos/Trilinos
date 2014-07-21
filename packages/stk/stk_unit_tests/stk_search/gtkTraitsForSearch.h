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
