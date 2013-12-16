#ifndef UNITTESTUTILS_HPP
#define UNITTESTUTILS_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <Geom_AxisAlignedBB.h>

typedef stk::search::IdentProc<int,int> Ident;
typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef stk::search::Box<double> Box;

template<class VolumeType>
VolumeType generateBoundingVolume(double x, double y, double z, double radius);

template<>
inline
Point generateBoundingVolume<Point>(double x, double y, double z, double /*radius*/)
{
  return Point(x,y,z);
}

template<>
inline
Sphere generateBoundingVolume<Sphere>(double x, double y, double z, double radius)
{
  return Sphere(Point(x,y,z),radius);
}

//       ------------
//      |            |
//      |      radius|
//      |      ------|
//      |            |
//      |            |
//       ------------
// width = 2*radius
template<>
inline
Box generateBoundingVolume<Box>(double x, double y, double z, double radius)
{
  Point min_corner(x-radius,y-radius,z-radius);
  Point max_corner(x+radius,y+radius,z+radius);
  return Box(min_corner,max_corner);
}

template <typename VolumeType>
std::pair<VolumeType, Ident> generateBoundingVolume(double x, double y, double z, double radius, int id, int proc)
{
  return std::make_pair(generateBoundingVolume<VolumeType>(x,y,z,radius), Ident(id,proc));
}

namespace stk {
namespace search {

inline bool intersects(const geometry::AxisAlignedBB& a, const geometry::AxisAlignedBB& b)
{
  return a.overlap(b);
}

}} // namespace stk::search

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


#endif
