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

namespace detail {

template <class REAL, class T>
static void assign_location_and_parametric_coords(const stk::math::Vec<REAL,3> & closestPt, const REAL closestParamX, const REAL closestParamY, T & result);

template <class REAL>
static void assign_location_and_parametric_coords(const stk::math::Vec<REAL,3> & closestPt, const REAL closestParamX, const REAL closestParamY, stk::math::Vec<REAL,3> & result)
{
  result = closestPt;
}

template <class REAL>
static void assign_location_and_parametric_coords(const stk::math::Vec<REAL,3> & closestPt, const REAL closestParamX, const REAL closestParamY, const std::tuple<stk::math::Vec<REAL,3> &, stk::math::Vec<REAL,2> &> & result)
{
  std::get<0>(result) = closestPt;
  std::get<1>(result)[0] = closestParamX;
  std::get<1>(result)[1] = closestParamY;
}

template<class REAL, class T>
void assign_vertex(const stk::math::Vec<REAL,3> & pt, const REAL x, const REAL y, T & result)
{
  assign_location_and_parametric_coords(pt, x, y, result);
}

template<class REAL, class T>
void assign_edge_AB_location(const stk::math::Vec<REAL,3> & a, const stk::math::Vec<REAL,3> & ab, const REAL d1, const REAL d3, T & result)
{
  REAL v = d1/(d1-d3);
  assign_location_and_parametric_coords(a + v * ab, v, 0.0, result);
}

template<class REAL, class T>
void assign_edge_AC_location(const stk::math::Vec<REAL,3> & a, const stk::math::Vec<REAL,3> & ac, const REAL d2, const REAL d6, T & result)
{
  REAL w(d2/(d2-d6));
  assign_location_and_parametric_coords(a + w * ac, 0.0, w, result);
}

template<class REAL, class T>
void assign_edge_BC_location(const stk::math::Vec<REAL,3> & b, const stk::math::Vec<REAL,3> & c, const REAL d3, const REAL d4, const REAL d5, const REAL d6, T & result)
{
  REAL w((d4-d3)/((d4-d3) + (d5-d6)));
  assign_location_and_parametric_coords(b + w * (c-b), 1.0-w, w, result);
}

template<class REAL, class T>
void assign_face_location(const stk::math::Vec<REAL,3> & a, const stk::math::Vec<REAL,3> & ab, const stk::math::Vec<REAL,3> & ac, const REAL va, const REAL vb, const REAL vc, T & result)
{
  REAL denom(1.0/(va+vb+vc));
  REAL v(vb*denom);
  REAL w(vc*denom);
  assign_location_and_parametric_coords(a + ab*v + ac*w, v, w, result);
}

}

template<typename REAL>
class CalcTriangle3 {
 public:
  using Vec3d = stk::math::Vec<REAL,3>;
  using Vec2d = stk::math::Vec<REAL,2>;

  struct ClosestPointIsInside
  {
    // The inequalities in closest_point_projection() are such that the edges and vertices are not considered to be inside.
    // If we want these points to be inside, we need to check for ==0.
    typedef bool ResultType;
    static void vertex_projection(const stk::math::Vec<REAL,3> &pt, const REAL x, const REAL y, ResultType & result)
      { result = false; }
    static void edge_AB_projection(const stk::math::Vec<REAL,3> &a, const stk::math::Vec<REAL,3> &ab, const REAL d1, const REAL d3, ResultType & result)
      { result = false; }
    static void edge_AC_projection(const stk::math::Vec<REAL,3> &a, const stk::math::Vec<REAL,3> &ac, const REAL d2, const REAL d6, ResultType & result)
      { result = false; }
    static void edge_BC_projection(const stk::math::Vec<REAL,3> &b, const stk::math::Vec<REAL,3> &c, const REAL d3, const REAL d4, const REAL d5, const REAL d6, ResultType & result)
      { result = false; }
    static void face_projection(const stk::math::Vec<REAL,3> &a, const stk::math::Vec<REAL,3> &ab, const stk::math::Vec<REAL,3> &ac, const REAL va, const REAL vb, const REAL vc, ResultType & result)
      { result = true; }
  };

  template<class T>
  struct ClosestPoint
  {
    typedef T ResultType;
    static void vertex_projection(const stk::math::Vec<REAL,3> &pt, const REAL x, const REAL y, ResultType & result)
      { detail::assign_vertex(pt, x, y, result); }
    static void edge_AB_projection(const stk::math::Vec<REAL,3> &a, const stk::math::Vec<REAL,3> &ab, const REAL d1, const REAL d3, ResultType & result)
      { detail::assign_edge_AB_location(a,ab,d1,d3, result); }
    static void edge_AC_projection(const stk::math::Vec<REAL,3> &a, const stk::math::Vec<REAL,3> &ac, const REAL d2, const REAL d6, ResultType & result)
      { detail::assign_edge_AC_location(a,ac,d2,d6, result); }
    static void edge_BC_projection(const stk::math::Vec<REAL,3> &b, const stk::math::Vec<REAL,3> &c, const REAL d3, const REAL d4, const REAL d5, const REAL d6, ResultType & result)
      { detail::assign_edge_BC_location(b,c,d3,d4,d5,d6, result); }
    static void face_projection(const stk::math::Vec<REAL,3> &a, const stk::math::Vec<REAL,3> &ab, const stk::math::Vec<REAL,3> &ac, const REAL va, const REAL vb, const REAL vc, ResultType & result)
      { detail::assign_face_location(a,ab,ac,va,vb,vc, result); }
  };

  using ClosestPointLocation = ClosestPoint<stk::math::Vec<REAL,3>>;
  using ClosestPointLocationAndParametricCoords = ClosestPoint<const std::tuple<stk::math::Vec<REAL,3> &, stk::math::Vec<REAL,2> &>>;

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

  ///  Compute the closest point on the triangle to an input point.
  static REAL distance_squared(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec3d& queryPt);
  static REAL distance_squared(const std::array<Vec3d,3> & triCoords, const Vec3d& queryPt)  { return distance_squared(triCoords[0], triCoords[1], triCoords[2], queryPt); }

  static Vec3d normal(const std::array<Vec3d,3> & triCoords) { return normal(triCoords[0], triCoords[1], triCoords[2]); }
  static Vec3d normal_dir(const std::array<Vec3d,3> & triCoords) { return normal_dir(triCoords[0], triCoords[1], triCoords[2]); }
  static REAL area(const std::array<Vec3d,3> & triCoords) { return area(triCoords[0], triCoords[1], triCoords[2]); }
  static Vec3d parametric_to_real_coords(const std::array<Vec3d,3> & triCoords, const Vec2d & paramPt) { return parametric_to_real_coords(triCoords[0], triCoords[1], triCoords[2], paramPt); }

  template<class ProjectionType>
  static void closest_point_projection(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d& p, typename ProjectionType::ResultType & result);

  static void closest_point(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec3d& queryPt, Vec3d& closestPt)
  {
    closest_point_projection<ClosestPointLocation>(p0, p1, p2, queryPt, closestPt);
  }

  static Vec3d closest_point(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec3d& queryPt)
  {
    Vec3d closestPt(stk::math::MemberInit::NONE);
    closest_point(p0, p1, p2, queryPt, closestPt);
    return closestPt;
  }

  static void closest_point(const std::array<Vec3d,3> & triCoords, const Vec3d& queryPt, Vec3d& closestPt)
  {
    closest_point(triCoords[0], triCoords[1], triCoords[2], queryPt, closestPt);
  }

  static Vec3d closest_point(const std::array<Vec3d,3> & triCoords, const Vec3d& queryPt)
  {
    return closest_point(triCoords[0], triCoords[1], triCoords[2], queryPt);
  }

  static void closest_point_and_parametric_coords(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec3d& queryPt, Vec3d& closestPt, Vec2d & paramPt)
  {
    closest_point_projection<ClosestPointLocationAndParametricCoords>(p0, p1, p2, queryPt, std::tie(closestPt, paramPt));
  }

  static void closest_point_and_parametric_coords(const std::array<Vec3d,3> & triCoords, const Vec3d& queryPt, Vec3d& closestPt, Vec2d & paramPt)
  {
    closest_point_and_parametric_coords(triCoords[0], triCoords[1], triCoords[2], queryPt, closestPt, paramPt);
  }

  static bool is_projection_of_point_inside_triangle(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec3d& queryPt)
  {
    bool result = false;
    closest_point_projection<ClosestPointIsInside>(p0, p1, p2, queryPt, result);
    return result;
  }
};

template<typename REAL>
inline REAL CalcTriangle3<REAL>::distance_squared(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2, const Vec3d& queryPt)
{
  return (queryPt-closest_point(p0,p1,p2, queryPt)).length_squared();
}

template<typename REAL>
template<class ProjectionType>
void CalcTriangle3<REAL>::closest_point_projection(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d& p, typename ProjectionType::ResultType & result)
{
  const Vec3d ab(b - a);
  const Vec3d ac(c - a);
  const Vec3d ap(p - a);

  const REAL d1(stk::math::Dot(ab, ap));
  const REAL d2(stk::math::Dot(ac, ap));

  if(d1 <= 0.0 && d2 <= 0.0)
  {
    ProjectionType::vertex_projection(a, 0., 0., result);
    return;
  }

  const Vec3d bp(p - b);
  const REAL d3(stk::math::Dot(ab, bp));
  const REAL d4(stk::math::Dot(ac, bp));

  if(d3 >= 0.0 && d4 <= d3)
  {
    ProjectionType::vertex_projection(b, 1., 0., result);
    return;
  }

  const REAL vc(d1 * d4 - d3 * d2);

  if(vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
  {
    ProjectionType::edge_AB_projection(a,ab,d1,d3, result);
    return;
  }

  const Vec3d cp(p - c);
  const REAL d5(Dot(ab, cp));
  const REAL d6(Dot(ac, cp));

  if(d6 >= 0.0 && d5 <= d6)
  {
    ProjectionType::vertex_projection(c, 0., 1., result);
    return;
  }

  const REAL vb(d5 * d2 - d1 * d6);

  if(vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
  {
    ProjectionType::edge_AC_projection(a,ac,d2,d6, result);
    return;
  }

  const REAL va(d3 * d6 - d5 * d4);

  if(va <= 0.0 && d4 >= d3 && d5 >= d6)
  {
     ProjectionType::edge_BC_projection(b,c,d3,d4,d5,d6, result);
     return;
  }

  ProjectionType::face_projection(a,ab,ac,va,vb,vc, result);
}

}
#endif // Akri_Triangle_h
