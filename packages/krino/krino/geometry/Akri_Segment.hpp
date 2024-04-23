// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Segment_h
#define Akri_Segment_h

#include <stk_math/StkVector.hpp>
#include <cmath>
#include <iostream>
#include <array>

namespace krino {

template<class REAL>
class CalcSegment3 {
 public:
  using Vec3d = stk::math::Vec<REAL,3>;

  static REAL length(const Vec3d & p0, const Vec3d & p1) { return (p1-p0).length(); }
  static REAL length_squared(const Vec3d & p0, const Vec3d & p1) { return (p1-p0).length_squared(); }

  //
  //  Find the closest point projection between point and the face.
  //
  static bool is_projection_of_point_inside_segment(const Vec3d & p0, const Vec3d & p1, const Vec3d &queryPt)
  {
    Vec3d edge_dir(p1-p0);
    REAL dotValA = Dot(edge_dir,(queryPt-p0));
    if(dotValA < 0.0) return false;
    REAL dotValB = Dot(edge_dir,(queryPt-p1));
    if(dotValB > 0.0) return false;
    return true;
  }
  static REAL closest_parametric_location(const Vec3d & p0, const Vec3d & p1, const Vec3d &queryPt)
  {
    Vec3d edge_dir(p1-p0);
    REAL dotValA = Dot(edge_dir,(queryPt-p0));
    if(dotValA <= 0.0) { return 0.0; }
    REAL dotValB = Dot(edge_dir,(queryPt-p1));
    if(dotValB >= 0.0) { return 1.0; }
    REAL lenSquared = edge_dir.length_squared();
    if(lenSquared == 0.0) { return 0.0; }
    return dotValA / lenSquared;
  }

  static void closest_point_and_parametric_coord(const Vec3d & p0, const Vec3d & p1, const Vec3d &queryPt, Vec3d &closestPt, REAL &paramLoc)
  {
    paramLoc = closest_parametric_location(p0, p1, queryPt);
    closestPt = (1.0-paramLoc)*p0 + paramLoc*p1;
  }

  static void closest_point(const Vec3d & p0, const Vec3d & p1, const Vec3d &queryPt, Vec3d &closestPt)
  {
    const REAL paramLoc = closest_parametric_location(p0, p1, queryPt);
    closestPt = (1.0-paramLoc)*p0 + paramLoc*p1;
  }

  static Vec3d closest_point(const Vec3d & p0, const Vec3d & p1, const Vec3d& queryPt)
  {
    Vec3d closestPt(stk::math::MemberInit::NONE);
    closest_point(p0, p1, queryPt, closestPt);
    return closestPt;
  }

  static void closest_point_and_parametric_coord(const std::array<Vec3d,2> & coords, const Vec3d &queryPt, Vec3d &closestPt, REAL &paramLoc) { closest_point_and_parametric_coord(coords[0], coords[1], queryPt, closestPt, paramLoc); }
  static void closest_point(const std::array<Vec3d,2> & coords, const Vec3d &queryPt, Vec3d &closestPt) { closest_point(coords[0], coords[1], queryPt, closestPt); }
  static Vec3d closest_point(const std::array<Vec3d,2> & coords, const Vec3d &queryPt) { closest_point(coords[0], coords[1], queryPt); }

  static REAL distance_squared(const Vec3d & p0, const Vec3d & p1, const Vec3d& queryPt)
  {
    return (queryPt - closest_point(p0, p1, queryPt)).length_squared();
  }

  static REAL length(const std::array<Vec3d,2> & coords) { return (coords[1]-coords[0]).length(); }
  static REAL length_squared(const std::array<Vec3d,2> & coords) { return (coords[1]-coords[0]).length_squared(); }
  static REAL distance_squared(const std::array<Vec3d,2> & coords, const Vec3d& queryPt) { return distance_squared(coords[0], coords[1], queryPt); }
};

} // namespace krino

#endif // Akri_Segment_h
