// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#ifndef STK_SEARCH_BOUNDINGBOX_HPP
#define STK_SEARCH_BOUNDINGBOX_HPP

#include <stk_search/Point.hpp>
#include <stk_search/Sphere.hpp>
#include <stk_search/Box.hpp>

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
template <typename T>
inline bool intersects(Sphere<T> const& a, Box<T> const& b)
{
  Point<T> const& ac = a.center();
  Point<T> const& bmin = b.min_corner();
  Point<T> const& bmax = b.max_corner();

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
inline bool intersects(Box<T> const& a, Sphere<T> const& b)
{ return intersects(b,a); }

// intersects: Box,Box
template <typename T>
inline bool intersects(Box<T> const& a, Box<T> const& b)
{
  Point<T> const& amin = a.min_corner();
  Point<T> const& amax = a.max_corner();

  Point<T> const& bmin = b.min_corner();
  Point<T> const& bmax = b.max_corner();

  // check that the boxes are not disjoint
  return !((amax[0] < bmin[0]) || (bmax[0] < amin[0])
        || (amax[1] < bmin[1]) || (bmax[1] < amin[1])
        || (amax[2] < bmin[2]) || (bmax[2] < amin[2]));

}

template <typename T, typename U>
inline void scale_by(Sphere<T> &s, U const& c)
{
  s.set_radius(s.radius()*c);
}

template <typename T, typename U>
inline void scale_by(Box<T> &b, U const& c)
{
  Point<T> & min_corner = b.min_corner();
  Point<T> & max_corner = b.max_corner();
  const U factor = (c-1)/2;
  for (int i=0; i<3; ++i) {
    const T d = factor*(max_corner[i] - min_corner[i]);
    min_corner[i] -= d;
    max_corner[i] += d;
  }
}

}} //namespace stk::search


#endif //STK_SEARCH_BOUNDINGBOX_HPP
