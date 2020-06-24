
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

#ifndef STK_SEARCH_BOUNDINGBOX_HPP
#define STK_SEARCH_BOUNDINGBOX_HPP

#include <stk_math/StkVector.hpp>
#include <stk_search/Box.hpp>
#include <stk_search/Plane.hpp>
#include <stk_search/Point.hpp>
#include <stk_search/Sphere.hpp>
#include <stk_math/StkMath.hpp>          // for stk::math::max, stk::math::min

namespace stk { namespace search {

// is_valid: Point
template <typename T>
inline bool is_valid(Point<T> const&)
{ return true; }

// is_valid: Sphere
template <typename T>
inline bool is_valid(Sphere<T> const& s)
{ return static_cast<T>(0) <= s.radius(); }

// is_valid: Box
template <typename T>
inline bool is_valid(Box<T> const& b)
{
  return   b.min_corner()[0] <= b.max_corner()[0]
        && b.min_corner()[1] <= b.max_corner()[1]
        && b.min_corner()[2] <= b.max_corner()[2] ;
}

//Point
template <typename T>
inline Point<T> const& min_corner(Point<T> const& p)
{ return p; }

template <typename T>
inline Point<T> const& max_corner(Point<T> const& p)
{ return p; }

template <typename T>
inline Point<T> const& center(Point<T> const& p)
{ return p; }

//Point,index
template <typename T>
inline T min_corner(Point<T> const& p, int index)
{ return p[index]; }

template <typename T>
inline T max_corner(Point<T> const& p, int index)
{ return p[index]; }

template <typename T>
inline T center(Point<T> const& p, int index)
{ return p[index]; }


//Sphere
template <typename T> inline Point<T> const min_corner(Sphere<T> const& s)
{
  const Point<T> & c = s.center();
  const T r = s.radius();
  return Point<T>(c[0]-r,c[1]-r,c[2]-r);
}

template <typename T> inline Point<T> const max_corner(Sphere<T> const& s)
{
  const Point<T> & c = s.center();
  const T r = s.radius();
  return Point<T>(c[0]+r,c[1]+r,c[2]+r);
}

template <typename T>
inline Point<T> const& center(Sphere<T> const& s)
{ return s.center(); }

//Sphere,index
template <typename T>
inline T min_corner(Sphere<T> const& s, int index)
{ return s.center()[index] - s.radius(); }

template <typename T>
inline T max_corner(Sphere<T> const& s, int index)
{ return s.center()[index] + s.radius(); }

template <typename T>
inline T center(Sphere<T> const& s, int index)
{ return s.center()[index]; }


//Box
template <typename T>
inline Point<T> const& min_corner(Box<T> const& b)
{ return b.min_corner(); }

template <typename T>
inline Point<T> const& max_corner(Box<T> const& b)
{ return b.max_corner(); }

template <typename T>
inline Point<T> const center(Box<T> const& b)
{
  const Point<T> & l = b.min_corner();
  const Point<T> & u = b.max_corner();
  return Point<T>( (l[0]+u[0])/2, (l[1]+u[1])/2, (l[2]+u[2])/2);
}

//Box,index
template <typename T>
inline T min_corner(Box<T> const& b, int index)
{ return b.min_corner()[index]; }

template <typename T>
inline T max_corner(Box<T> const& b, int index)
{ return b.max_corner()[index]; }

template <typename T>
inline T center(Box<T> const& b, int index)
{ return (b.min_corner()[index] + b.max_corner()[index])/2; }


// intersects: Point,Point
template <typename T>
inline bool intersects(Point<T> const& a, Point<T> const& b)
{ return (a==b); }


// intersects: Point,Sphere
template <typename T>
inline bool intersects(Point<T> const& a, Sphere<T> const& b)
{
  const T dist2 =   (a[0]-b.center()[0])*(a[0]-b.center()[0])
                  + (a[1]-b.center()[1])*(a[1]-b.center()[1])
                  + (a[2]-b.center()[2])*(a[2]-b.center()[2]);
  return (dist2 <= b.radius()*b.radius());
}

// intersects: Sphere,Point
template <typename T>
inline bool intersects(Sphere<T> const& a, Point<T> const& b)
{ return intersects(b,a); }

// intersects: Point,Box
template <typename T>
inline bool intersects(Point<T> const& a, Box<T> const& b)
{
  return b.min_corner()[0] <= a[0] && a[0] <= b.max_corner()[0]
      && b.min_corner()[1] <= a[1] && a[1] <= b.max_corner()[1]
      && b.min_corner()[2] <= a[2] && a[2] <= b.max_corner()[2];
}

// intersects: Box,Point
template <typename T>
inline bool intersects(Box<T> const& a, Point<T> const& b)
{ return intersects(b,a); }

// intersects: Sphere,Sphere
template <typename T>
inline bool intersects(Sphere<T> const& a, Sphere<T> const& b)
{
  const Point<T> & ac = a.center();
  const Point<T> & bc = b.center();
  const T r2 = (a.radius()+b.radius())*(a.radius()+b.radius());
  const T dist2 =  (ac[0]-bc[0])*(ac[0]-bc[0])
                 + (ac[1]-bc[1])*(ac[1]-bc[1])
                 + (ac[2]-bc[2])*(ac[2]-bc[2]);
  return dist2 < r2;
}

// intersects: Sphere,Box
template <typename T1, typename T2>
inline bool intersects(Sphere<T1> const& a, Box<T2> const& b)
{
  Point<T1> const& ac   = a.center();
  const T1         ar   = a.radius(); 
  T1 distToBoxSquared = 0.0;
  for(unsigned idir=0; idir<3; ++idir) {
      const T2 boxMin = b.min_corner()[idir];
      const T2 boxMax = b.max_corner()[idir];
      const T1 sphereCenter = ac[idir];

      if(       sphereCenter + ar < boxMin) {
          return false;
      } else if(sphereCenter - ar > boxMax) {
          return false;
      } else if(sphereCenter      < boxMin) {
          T1 dist = boxMin-sphereCenter;
          distToBoxSquared += dist*dist;
      } else if(sphereCenter      > boxMax) {
          T1 dist = sphereCenter-boxMax;
          distToBoxSquared += dist*dist;
      }
  }
   return (distToBoxSquared <= ar*ar);
}

// intersects: Box,Sphere
template <typename T1, typename T2>
inline bool intersects(Box<T1> const& a, Sphere<T2> const& b)
{ return intersects(b,a); }

// intersects: Box,Box
template <typename T1, typename T2>
inline bool intersects(Box<T1> const& a, Box<T2> const& b)
{
  Point<T1> const& amin = a.min_corner();
  Point<T1> const& amax = a.max_corner();

  Point<T2> const& bmin = b.min_corner();
  Point<T2> const& bmax = b.max_corner();

  // check that the boxes are not disjoint
  return !((amax[0] < bmin[0]) || (bmax[0] < amin[0])
        || (amax[1] < bmin[1]) || (bmax[1] < amin[1])
        || (amax[2] < bmin[2]) || (bmax[2] < amin[2]));

}

// intersects: Plane,Box
template <typename T1, typename T2>
inline bool intersects(Plane<T1> const& p, Box<T2> const& b)
{

  std::array<T2,6> points = {b.get_x_min(), b.get_x_max(), b.get_y_min(), b.get_y_max(), b.get_z_min(), b.get_z_max()};
  int previousWhichSide = p.WhichSide(Point<T2>(points[0], points[0], points[0]));
  for (int x=0; x<2; ++x){
    for (int y=0; y<2; ++y){
      for (int z=0; z<2; ++z){
        int whichSide = p.WhichSide(Point<T2>(points[x], points[y], points[z]));
        if (whichSide != previousWhichSide){
          return true;
        }
        else if(whichSide == 0){
          return true;
        }
      }
    }
  }
  return false;
}

template <typename T1, typename T2>
inline bool intersects(Box<T1> const& b, Plane<T2> const& p)
{
  return intersects(p,b);
}

template <typename T, typename U>
inline void scale_by(Point<T> &p, U const& mult_fact, U const& add_fact = 0)
{
}

template <typename T, typename U>
inline void scale_by(Sphere<T> &s, U const& mult_fact, U const& add_fact = 0)
{
  s.set_radius(s.radius()*mult_fact + add_fact);
}

template <typename T, typename U>
inline void scale_by(Box<T> &b, U const& mult_fact, U const& add_fact = 0)
{
  Point<T> & min_corner = b.min_corner();
  Point<T> & max_corner = b.max_corner();
  const U factor = (mult_fact-1)/2;
  for (int i=0; i<3; ++i) {
    const T d = factor*(max_corner[i] - min_corner[i]) + add_fact;
    min_corner[i] -= d;
    max_corner[i] += d;
  }
}

template <typename T1, typename T2>
KOKKOS_INLINE_FUNCTION void add_to_box(Box<T1> &box, const Box<T2>& addBox) {
  box.set_box(stk::math::min(box.get_x_min(), static_cast<T1>(addBox.get_x_min())),
              stk::math::min(box.get_y_min(), static_cast<T1>(addBox.get_y_min())),
              stk::math::min(box.get_z_min(), static_cast<T1>(addBox.get_z_min())),
              stk::math::max(box.get_x_max(), static_cast<T1>(addBox.get_x_max())),
              stk::math::max(box.get_y_max(), static_cast<T1>(addBox.get_y_max())),
              stk::math::max(box.get_z_max(), static_cast<T1>(addBox.get_z_max())));
}

template <typename T1, typename T2>
KOKKOS_INLINE_FUNCTION void add_to_box(Box<T1> &box, const Sphere<T2>& addBox) {
  box.set_box(stk::math::min(box.get_x_min(), addBox.get_x_min()), 
              stk::math::min(box.get_y_min(), addBox.get_y_min()), 
              stk::math::min(box.get_z_min(), addBox.get_z_min()), 
              stk::math::max(box.get_x_max(), addBox.get_x_max()),
              stk::math::max(box.get_y_max(), addBox.get_y_max()), 
              stk::math::max(box.get_z_max(), addBox.get_z_max())); 
}

template <typename T1, typename T2>
KOKKOS_INLINE_FUNCTION void add_to_box(Box<T1> &box, const Point<T2>& addBox) {
  box.set_box(stk::math::min(box.get_x_min(), addBox.get_x_min()), 
              stk::math::min(box.get_y_min(), addBox.get_y_min()), 
              stk::math::min(box.get_z_min(), addBox.get_z_min()), 
              stk::math::max(box.get_x_max(), addBox.get_x_max()),
              stk::math::max(box.get_y_max(), addBox.get_y_max()), 
              stk::math::max(box.get_z_max(), addBox.get_z_max())); 
}

// This algorithm is based off the minimum circle for a triangle blog post
// by Christer Ericson at http://realtimecollisiondetection.net/blog/?p=20
template <typename NumT>
Sphere<NumT> minimumBoundingSphere(const Point<NumT>& ptA, const Point<NumT>& ptB, const Point<NumT>& ptC)
{
    typedef stk::math::Vector3d Vec;

    Vec a = Vec(ptA[0], ptA[1], ptA[2]);
    Vec b = Vec(ptB[0], ptB[1], ptB[2]);
    Vec c = Vec(ptC[0], ptC[1], ptC[2]);

    Vec AB = b - a;
    Vec AC = c - a;

    NumT dotABAB = Dot(AB, AB);
    NumT dotABAC = Dot(AB, AC);
    NumT dotACAC = Dot(AC, AC);

    NumT d = 2.0*(dotABAB*dotACAC - dotABAC*dotABAC);
    Vec referencePt = a;

    Vec center;
    if (std::abs(d) <= 100 * std::numeric_limits<NumT>::epsilon()) {
        // a, b, and c lie on a line. Circle center is center of AABB of the
        // points, and radius is distance from circle center to AABB corner
        Box<NumT> aabb = Box<NumT>(Point<NumT>(a[0],a[1],a[2]), Point<NumT>(b[0],b[1],b[2]));
        add_to_box(aabb, Point<NumT>(c[0],c[1],c[2]));

        Point<NumT> minCornerPt = aabb.min_corner();
        Point<NumT> maxCornerPt = aabb.max_corner();
        Vec minCorner = Vec(minCornerPt[0], minCornerPt[1], minCornerPt[2]);
        Vec maxCorner = Vec(maxCornerPt[0], maxCornerPt[1], maxCornerPt[2]);
        center = 0.5 * (minCorner + maxCorner);
        referencePt = minCorner;
    } else {
        NumT s = (dotABAB*dotACAC - dotACAC*dotABAC) / d;
        NumT t = (dotACAC*dotABAB - dotABAB*dotABAC) / d;
        // s controls height over AC, t over AB, (1-s-t) over BC
        if (s <= 0.0f) {
            center = 0.5 * (a + c);
        } else if (t <= 0.0f) {
            center = 0.5 * (a + b);
        } else if (s + t >= 1.0) {
            center = 0.5 * (b + c);
            referencePt = b;
        } else {
            center = a + s*(b - a) + t*(c - a);
        }
    }
    NumT radius = std::sqrt(Dot(center - referencePt, center - referencePt));

    return Sphere<NumT>(Point<NumT> (center[0],center[1],center[2]), radius);
}


}} //namespace stk::search


#endif //STK_SEARCH_BOUNDINGBOX_HPP
