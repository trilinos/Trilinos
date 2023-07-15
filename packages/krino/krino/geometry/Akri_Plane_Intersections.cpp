// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Plane_Intersections.hpp>
#include <stk_math/StkPlane.hpp>
#include <stk_math/StkVector.hpp>

namespace krino
{

bool
find_intersection_of_three_planes(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2, stk::math::Vector3d & point)
{
  const stk::math::Vector3d & normal0 = plane0.normal();
  const stk::math::Vector3d & normal1 = plane1.normal();
  const stk::math::Vector3d & normal2 = plane2.normal();

  const double a00 = normal0[0];
  const double a01 = normal0[1];
  const double a02 = normal0[2];
  const double a10 = normal1[0];
  const double a11 = normal1[1];
  const double a12 = normal1[2];
  const double a20 = normal2[0];
  const double a21 = normal2[1];
  const double a22 = normal2[2];
  const double b0 = plane0.constant();
  const double b1 = plane1.constant();
  const double b2 = plane2.constant();
  const double det = a00*(a22*a11-a21*a12)-a10*(a22*a01-a21*a02)+a20*(a12*a01-a11*a02);

  if (det == 0)
    return false;

  const double x =( b0*(a22*a11-a21*a12)-b1*(a22*a01-a21*a02)+b2*(a12*a01-a11*a02))/det;
  const double y =(-b0*(a22*a10-a20*a12)+b1*(a22*a00-a20*a02)-b2*(a12*a00-a10*a02))/det;
  const double z =( b0*(a21*a10-a20*a11)-b1*(a21*a00-a20*a01)+b2*(a11*a00-a10*a01))/det;

  point = stk::math::Vector3d{x,y,z};
  return true;
}

stk::math::Vector3d
triangle_parametric_coordinates_of_projected_point(const std::array<stk::math::Vector3d,3> & triCoords, const stk::math::Vector3d & pt)
{
  const stk::math::Vector3d v1 = triCoords[1] - triCoords[0];
  const stk::math::Vector3d v2 = triCoords[2] - triCoords[0];
  const stk::math::Vector3d normal = Cross(v1,v2);
  const double invNormalMag2 = 1./normal.length_squared();
  const stk::math::Vector3d diff = pt - triCoords[0];
  stk::math::Vector3d paramPt;
  paramPt[0] = Dot(Cross(diff,v2),normal) * invNormalMag2;
  paramPt[1] = Dot(Cross(v1,diff),normal) * invNormalMag2;
  paramPt[2] = 0.;
  return paramPt;
}

static bool within_tet_bounds(const stk::math::Vector3d & pt)
{
  const double tol = 1.e-13;
  return (pt[0] > -tol) &&
         (pt[1] > -tol) &&
         (pt[2] > -tol) &&
         ((1.-pt[0]-pt[1]-pt[2]) > -tol);
}

bool find_intersection_of_three_planes_within_tet(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2, stk::math::Vector3d & intersectionPoint)
{
  return find_intersection_of_three_planes(plane0, plane1, plane2, intersectionPoint) && within_tet_bounds(intersectionPoint);
}

static bool
find_intersection_of_two_planes_and_coordinate_plane(const int coordinatePlane, const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, stk::math::Vector3d & point)
{
  const stk::math::Vector3d & normal0 = plane0.normal();
  const stk::math::Vector3d & normal1 = plane1.normal();

  const std::array<std::array<int,2>,3> activeCoordsByPlane = {{ {{1,2}}, {{0,2}}, {{0,1}} }};
  const int i0 = activeCoordsByPlane[coordinatePlane][0];
  const int i1 = activeCoordsByPlane[coordinatePlane][1];

  const double a00 = normal0[i0];
  const double a01 = normal0[i1];
  const double a10 = normal1[i0];
  const double a11 = normal1[i1];
  const double b0 = plane0.constant();
  const double b1 = plane1.constant();
  const double det = a00*a11-a10*a01;

  if (det == 0)
    return false;

  point[i0] = (b0*a11-b1*a01)/det;
  point[i1] = (b1*a00-b0*a10)/det;
  point[coordinatePlane] = 0.;

  return true;
}

bool find_intersection_of_two_planes_and_side_of_tet(const int side, const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, stk::math::Vector3d & intersectionPoint)
{
  bool intersectsPlaneThatCoincidesWithSideOfTet = false;
  if (side == 1)
  {
    const stk::math::Plane3d side1Plane{stk::math::Vector3d{1.,0.,0.},stk::math::Vector3d{0.,0.,1.},stk::math::Vector3d{0.,1.,0.}};
    intersectsPlaneThatCoincidesWithSideOfTet = find_intersection_of_three_planes(plane0, plane1, side1Plane, intersectionPoint);
  }
  else
  {
    const int coordinatePlane = (side == 0) ? 1 : ((side == 2) ? 0 : 2);
    intersectsPlaneThatCoincidesWithSideOfTet = find_intersection_of_two_planes_and_coordinate_plane(coordinatePlane, plane0, plane1, intersectionPoint);
  }

  return intersectsPlaneThatCoincidesWithSideOfTet && within_tet_bounds(intersectionPoint);

}

bool find_intersection_of_two_2D_planes(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, stk::math::Vector3d & intersectionPoint)
{
  return find_intersection_of_two_planes_and_coordinate_plane(2, plane0, plane1, intersectionPoint);
}

bool within_tri_bounds(const stk::math::Vector3d & triangleParamCoords)
{
  const double tol = 1.e-13;
  return (triangleParamCoords[0] > -tol) &&
         (triangleParamCoords[1] > -tol) &&
         ((1.-triangleParamCoords[0]-triangleParamCoords[1]) > -tol);
}

bool find_intersection_of_two_2D_planes_within_tri(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, stk::math::Vector3d & intersectionPoint)
{
  return find_intersection_of_two_planes_and_coordinate_plane(2, plane0, plane1, intersectionPoint) && within_tri_bounds(intersectionPoint);
}

} // namespace krino



