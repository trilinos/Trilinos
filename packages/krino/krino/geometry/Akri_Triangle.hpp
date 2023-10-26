// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Triangle_h
#define Akri_Triangle_h

#include <iostream>
#include <assert.h>
#include <cmath>
#include <array>

#include <stk_math/StkVector.hpp>

namespace krino {

enum class ProjectionType{NODE, EDGE, FACE, NUM_PROJ_TYPE};

template<typename REAL>
class CalcTriangle3 {
 public:
  using Vec3d = stk::math::Vec<REAL,3>;
  using Vec2d = stk::math::Vec<REAL,2>;
  static Vec3d normal(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2) { return normal_dir(p0,p1,p2).unit_vector(); } /// Unit vector
  static Vec3d normal_dir(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2) { return Cross(p1-p0,p2-p0); } /// Non-unit normal (faster)
  static REAL area(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2) { return 0.5*normal_dir(p0,p1,p2).length(); }

  static Vec3d parametric_to_real_coords(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec2d & paramPt)
    {
      // Avoids temporary constructions
      return Vec3d(
          p0[0]*(1.0-paramPt[0]-paramPt[1]) + p1[0]*paramPt[0] + p2[0]*paramPt[1],
          p0[1]*(1.0-paramPt[0]-paramPt[1]) + p1[1]*paramPt[0] + p2[1]*paramPt[1],
          p0[2]*(1.0-paramPt[0]-paramPt[1]) + p1[2]*paramPt[0] + p2[2]*paramPt[1]
          );
    }

  ///  Compute the closest point on the triangle to an input point.  Return the type of projection, i.e.,
  ///  is the projected point on a node, edge, surface of the triangle.  Optionally calculate parametric coordinates.
  template<class T = std::nullptr_t>
  static ProjectionType closest_point(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec3d& queryPt, Vec3d& ClosestPt, T paramPt = nullptr);
  static REAL distance_squared(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec3d& queryPt);


  static Vec3d normal(const std::array<Vec3d,3> & triCoords) { return normal(triCoords[0], triCoords[1], triCoords[2]); }
  static Vec3d normal_dir(const std::array<Vec3d,3> & triCoords) { return normal_dir(triCoords[0], triCoords[1], triCoords[2]); }
  static REAL area(const std::array<Vec3d,3> & triCoords) { return area(triCoords[0], triCoords[1], triCoords[2]); }
  static Vec3d parametric_to_real_coords(const std::array<Vec3d,3> & triCoords, const Vec2d & paramPt) { return parametric_to_real_coords(triCoords[0], triCoords[1], triCoords[2], paramPt); }
  template<class T = std::nullptr_t>
  static ProjectionType closest_point(const std::array<Vec3d,3> & triCoords, const Vec3d& queryPt, Vec3d& closestPt, T paramPt = nullptr) { return closest_point(triCoords[0], triCoords[1], triCoords[2], queryPt, closestPt, paramPt); }
  static REAL distance_squared(const std::array<Vec3d,3> & triCoords, const Vec3d& queryPt)  { return distance_squared(triCoords[0], triCoords[1], triCoords[2], queryPt); }
};

namespace detail {
template <class REAL, class T>
struct assign_parametric_coords {
  static void apply(const REAL & x, const REAL & y, T paramPt)
  {
    static_assert(std::is_pointer<T>::value, "Expecting pointer");
    paramPt[0] = x;
    paramPt[1] = y;
  }
};
template <class REAL>
struct assign_parametric_coords<REAL, std::nullptr_t>  {
  static void apply(const REAL, const REAL, std::nullptr_t) {}
};
}

template<typename REAL>
inline REAL CalcTriangle3<REAL>::distance_squared(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec3d& queryPt)
{
  Vec3d closestPt(stk::math::MemberInit::NONE);
  closest_point(p0,p1,p2, queryPt, closestPt);
  return (queryPt-closestPt).length_squared();
}

//  Adapted from closest face projection from "Real time Collision Detection" text, highly optimized
template<typename REAL>
template<class T>
ProjectionType CalcTriangle3<REAL>::closest_point(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d& p, Vec3d& closestPt, T paramPt)
{
  const Vec3d ab(b-a);
  const Vec3d ac(c-a);
  const Vec3d ap(p-a);

  REAL d1(Dot(ab,ap));
  REAL d2(Dot(ac,ap));

  if(d1 <= 0.0 && d2 <= 0.0) {
    closestPt = a;
    detail::assign_parametric_coords<REAL,T>::apply(0.0, 0.0, paramPt);
    return ProjectionType::NODE;
   }

  Vec3d bp(p-b);
  REAL d3(Dot(ab,bp));
  REAL d4(Dot(ac,bp));

  if(d3 >= 0.0 && d4  <= d3) {
    closestPt = b;
    detail::assign_parametric_coords<REAL,T>::apply(1.0, 0.0, paramPt);
    return ProjectionType::NODE;
  }

  REAL vc(d1*d4-d3*d2);
  if(vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    REAL v = d1/(d1-d3);
    closestPt = a + v * ab;
    detail::assign_parametric_coords<REAL,T>::apply(v, 0.0, paramPt);
    return ProjectionType::EDGE;
  }

  Vec3d cp(p-c);
  REAL d5(Dot(ab,cp));
  REAL d6(Dot(ac,cp));
  if(d6 >= 0.0 && d5 <= d6) {
    closestPt = c;
    detail::assign_parametric_coords<REAL,T>::apply(0.0, 1.0, paramPt);
    return ProjectionType::NODE;
  }

  REAL vb(d5*d2 - d1*d6);
  if(vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    REAL w(d2/(d2-d6));
    closestPt = a + w * ac;
    detail::assign_parametric_coords<REAL,T>::apply(0.0, w, paramPt);
    return ProjectionType::EDGE;
  }

  REAL va(d3*d6 - d5*d4);
  if(va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0) {
    REAL w((d4-d3)/((d4-d3) + (d5-d6)));
    closestPt =  b + w * (c-b);
    detail::assign_parametric_coords<REAL,T>::apply(1.0-w, w, paramPt);
    return ProjectionType::EDGE;
  }

  REAL denom(1.0/(va+vb+vc));
  REAL v(vb*denom);
  REAL w(vc*denom);
  closestPt = a + ab*v + ac*w;
  detail::assign_parametric_coords<REAL,T>::apply(v, w, paramPt);
  return ProjectionType::FACE;
}

}
#endif // Akri_Triangle_h
