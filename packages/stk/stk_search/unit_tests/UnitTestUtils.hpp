#ifndef UNITTESTUTILS_HPP
#define UNITTESTUTILS_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

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


#endif
