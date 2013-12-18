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

// traits for ::geometry::AxisAlignedBB
template <> struct tag< ::geometry::AxisAlignedBB > { typedef box_tag type; };
template <> struct point_type< ::geometry::AxisAlignedBB > { typedef stk::search::Point<Real> type; };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, min_corner, 0 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_x_min(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, min_corner, 1 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_y_min(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, min_corner, 2 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_z_min(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, max_corner, 0 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_x_max(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, max_corner, 1 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_y_max(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, max_corner, 2 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_z_max(); } };

}}} // namespace boost::geometry::traits

namespace stk {
namespace search {

inline bool intersects(const geometry::AxisAlignedBB& a, const geometry::AxisAlignedBB& b)
{
  return a.overlap(b);
}

}} // namespace stk::search

#endif /* GTKTRAITSFORSEARCH_H_ */
