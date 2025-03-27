/*
 * Akri_IntersectionUtils.cpp
 *
 *  Created on: Dec 10, 2024
 *      Author: drnoble
 */
#include <stk_math/StkVector.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_Sign.hpp>
#include <Akri_Triangle.hpp>
#include <iostream>

namespace krino {

template <typename VecType, size_t N>
bool are_all_components_lo(const VecType & bboxMin, const std::array<stk::math::Vector3d,N> & points, const unsigned comp)
{
  for (size_t i=0; i<N; ++i)
    if (points[i][comp] >= bboxMin[comp])
      return false;
  return true;
}

template <typename VecType, size_t N>
bool are_all_components_hi(const VecType & bboxMax, const std::array<stk::math::Vector3d,N> & points, const unsigned comp)
{
  for (size_t i=0; i<N; ++i)
    if (points[i][comp] <= bboxMax[comp])
      return false;
  return true;
}

template <size_t N>
bool does_bounding_box_intersect_simplex(const BoundingBox & bbox, const std::array<stk::math::Vector3d,N> & points, const unsigned ndim)
{
  for (unsigned comp=0; comp<ndim; ++comp)
    if (are_all_components_lo(bbox.get_min(), points, comp) ||
        are_all_components_hi(bbox.get_max(), points, comp))
      return false;
  return true;
}

bool does_bounding_box_intersect_tetrahedron(const BoundingBox & bbox, const std::array<stk::math::Vector3d,4> & tetNodes)
{
  return does_bounding_box_intersect_simplex(bbox, tetNodes, 3);
}

bool does_bounding_box_intersect_triangle_3d(const BoundingBox & bbox, const std::array<stk::math::Vector3d,3> & triNodes)
{
  return does_bounding_box_intersect_simplex(bbox, triNodes, 3);
}

bool does_bounding_box_intersect_edge_2d(const BoundingBox & bbox, const std::array<stk::math::Vector3d,2> & edgeNodes)
{
  return does_bounding_box_intersect_simplex(bbox, edgeNodes, 2);
}

bool does_tetrahedon_contain_point(const std::array<stk::math::Vector3d,4> & tetNodes, const stk::math::Vector3d & pt)
{
  const stk::math::Vector3d paramPt = get_parametric_coordinates_of_point(tetNodes, pt);
  for (int i=0; i<3; ++i)
    if (paramPt[i] < 0.)
      return false;
  return (paramPt[0]+paramPt[1]+paramPt[2]) <= 1.;
}

bool are_parametric_coordinates_within_simplex(const stk::math::Vector3d & paramCoords)
{
  for (int i=0; i<3; ++i)
    if (paramCoords[i] < 0.)
      return false;
  return (paramCoords[0]+paramCoords[1]+paramCoords[2]) <= 1.;
}

bool does_enlarged_tetrahedon_contain_point(const std::array<stk::math::Vector3d,4> & tetNodes, const stk::math::Vector3d & pt)
{
  constexpr double expand {1e-10};
  const stk::math::Vector3d centroid = 1./4.*(tetNodes[0]+tetNodes[1]+tetNodes[2]+tetNodes[3]);
  const std::array<stk::math::Vector3d,4> enlargedTetNodes{
    tetNodes[0] + expand*(tetNodes[0]-centroid),
    tetNodes[1] + expand*(tetNodes[1]-centroid),
    tetNodes[2] + expand*(tetNodes[2]-centroid),
    tetNodes[3] + expand*(tetNodes[3]-centroid)
  };
  return does_tetrahedon_contain_point(enlargedTetNodes, pt);
}

void append_tetrahedron_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,4> & tetNodes, const stk::math::Vector3d & candidate, std::vector<stk::math::Vector3d> & intParamCoords)
{
  const stk::math::Vector3d paramPt = get_parametric_coordinates_of_point(tetNodes, candidate);
  if (are_parametric_coordinates_within_simplex(paramPt))
  {
    intParamCoords.push_back(paramPt);
  }
}

std::pair<bool,stk::math::Vector3d> is_point_in_triangle_3d_and_parametric_coordinates(const std::array<stk::math::Vector3d,3> & triCoords, const double triArea, const stk::math::Vector3d& triNormal, const stk::math::Vector3d& pt)
{
  const double signedAreaABP = CalcTriangle3<double>::signed_area(triNormal, triCoords[0], triCoords[1], pt);
  const double signedAreaAPC = CalcTriangle3<double>::signed_area(triNormal, triCoords[0], pt, triCoords[2]);

  const stk::math::Vector3d paramCoords(signedAreaAPC / triArea, signedAreaABP / triArea, 0.);

  return std::make_pair(are_parametric_coordinates_within_simplex(paramCoords), paramCoords);
}

std::pair<bool,stk::math::Vector3d> is_point_in_triangle_3d_and_parametric_coordinates(const std::array<stk::math::Vector3d,3> & triCoords, const stk::math::Vector3d& pt)
{
  const auto [triArea, triNormal] = CalcTriangle3<double>::area_and_normal(triCoords);
  return is_point_in_triangle_3d_and_parametric_coordinates(triCoords, triArea, triNormal, pt);
}

std::pair<bool,stk::math::Vector3d> does_edge_intersect_triangle_and_parametric_coordinates_of_intersection(const std::array<stk::math::Vector3d,3> & triNodes,
    const double triArea, // passed in, rather than recomputed, for performance
    const stk::math::Vector3d & triNormal, // passed in, rather than recomputed, for performance
    const stk::math::Vector3d & edgePt0,
    const stk::math::Vector3d & edgePt1)
{
  const double scaledDist0 = Dot(triNormal, edgePt0-triNodes[0]);
  const double scaledDist1 = Dot(triNormal, edgePt1-triNodes[0]);
  if (!sign_change(scaledDist0, scaledDist1))
    return std::make_pair(false, stk::math::Vector3d::ZERO);

  const double loc = scaledDist0 / (scaledDist0-scaledDist1);
  const stk::math::Vector3d planeIntLoc = (1.-loc)*edgePt0 + loc*edgePt1;

  return is_point_in_triangle_3d_and_parametric_coordinates(triNodes, triArea, triNormal, planeIntLoc);
}

void append_triangle_edge_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,3> & triNodes,
    const double triArea, // passed in, rather than recomputed, for performance
    const stk::math::Vector3d & triNormal, // passed in, rather than recomputed, for performance
    const stk::math::Vector3d & edgePt0,
    const stk::math::Vector3d & edgePt1,
    std::vector<stk::math::Vector3d> & intParamCoords)
{
  const auto & [doesIntersect, paramCoords] = does_edge_intersect_triangle_and_parametric_coordinates_of_intersection(triNodes, triArea, triNormal, edgePt0, edgePt1);
  if (doesIntersect)
    intParamCoords.push_back(paramCoords);
}

}


