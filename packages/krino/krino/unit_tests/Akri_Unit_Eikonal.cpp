// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <stk_math/StkVector.hpp>
#include <Akri_Eikonal_Calc.hpp>
#include <stk_math/StkPlane.hpp>

namespace krino {

static constexpr double farDistance = std::numeric_limits<double>::max();

void write_cubit_tet(const std::array<stk::math::Vector3d,4> & coords)
{
  for (auto && x : coords)
    std::cout << "create vertex " << x[0] << " " << x[1] << " " << x[2] << std::endl;
  int nodeId = 1;
  for (int i=0; i<4; ++i)
    std::cout << "create node vertex " << i+nodeId << std::endl;
  std::cout << "create tet node " << nodeId << " " << nodeId+1 << " " << nodeId+2 << " " << nodeId+3 << std::endl;
}

void expect_triangle_eikonal_solution(const std::array<stk::math::Vector3d,3> & x, const std::array<double,2> & d, const double speed, const double gold, const double far = farDistance)
{
  const int sign = 1;
  const double resultCase = eikonal_solve_triangle(x, d, sign, far, speed);
  EXPECT_NEAR(gold, resultCase, 1.e-8);
}

void expect_triangle_distance(const std::array<stk::math::Vector3d,3> & x, const std::array<double,2> & d, const double gold)
{
  const double speed = 1.0;
  expect_triangle_eikonal_solution(x,d,speed,gold);
}

void expect_tetrahedron_eikonal_solution(const std::array<stk::math::Vector3d,4> & x, const std::array<double,3> & d, const double speed, const double gold, const double far = farDistance)
{
  const int sign = 1;
  const double resultCase = eikonal_solve_tetrahedron(x, d, sign, far, speed);
  EXPECT_NEAR(gold, resultCase, 1.e-8);
}

void expect_tetrahedron_distance_and_extension_speed(const std::array<stk::math::Vector3d,4> & x, const std::array<double,3> & d, const std::array<double,3> & extSpeed, const double goldDist, const double goldExtSpeed)
{
  const int sign = 1;
  const auto [nodeDist, nodeSpeed] = eikonal_solve_tetrahedron_for_distance_and_extension_speed(x, d, extSpeed, sign, farDistance);
  EXPECT_NEAR(goldDist, nodeDist, 1.e-8);
  EXPECT_NEAR(goldExtSpeed, nodeSpeed, 1.e-8);
}

void expect_triangle_distance_and_extension_speed(const std::array<stk::math::Vector3d,3> & x, const std::array<double,2> & d, const std::array<double,2> & extSpeed, const double goldDist, const double goldExtSpeed)
{
  const int sign = 1;
  const auto [nodeDist, nodeSpeed] = eikonal_solve_triangle_for_distance_and_extension_speed(x, d, extSpeed, sign, farDistance);
  EXPECT_NEAR(goldDist, nodeDist, 1.e-8);
  EXPECT_NEAR(goldExtSpeed, nodeSpeed, 1.e-8);
}

void expect_tetrahedron_distance(const std::array<stk::math::Vector3d,4> & x, const std::array<double,3> & d, const double gold)
{
  const double speed = 1.0;
  expect_tetrahedron_eikonal_solution(x,d,speed,gold);
}

std::array<stk::math::Vector3d,3> get_regular_triangle_coordinates()
{
  return std::array<stk::math::Vector3d,3>{{ {-0.5,0.,0.}, {0.5,0.,0.}, {0.,0.5*std::sqrt(3.0),0.} }};
}

std::array<stk::math::Vector3d,3> get_right_triangle_coordinates()
{
  return std::array<stk::math::Vector3d,3>{{ {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.} }};;
}

std::array<stk::math::Vector3d,4> get_regular_tetrahedron_coordinates()
{
  return std::array<stk::math::Vector3d,4>{{ {1.,1.,1.}, {-1.,1.,-1.}, {1.,-1.,-1.}, {-1.,-1.,1.} }};
}

std::array<stk::math::Vector3d,4> get_right_tetrahedron_coordinates()
{
  return std::array<stk::math::Vector3d,4>{{ {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} }};;
}

template<size_t SIZE>
double get_edge_length(const std::array<stk::math::Vector3d,SIZE> & coords, const int node0, const int node1)
{
  return (coords[node0]-coords[node1]).length();
}

double get_regular_triangle_edge_length()
{
  return get_edge_length(get_regular_triangle_coordinates(), 0, 1);
}

double get_right_triangle_edge_length()
{
  return get_edge_length(get_right_triangle_coordinates(), 0, 1);
}

double get_regular_tetrahedron_edge_length()
{
  return get_edge_length(get_regular_tetrahedron_coordinates(), 0, 1);
}

double get_right_tetrahedron_edge_length()
{
  return get_edge_length(get_right_tetrahedron_coordinates(), 0, 1);
}

TEST(Eikonal_Calc, whenDecreasingInputNbrValues_resultingDistanceDecreases)
{
  // This test was used to drive a fix for an issue where a decrease in the distance function for that neighbors
  // could produce an increased distance at a node.  This was caused by a causal check that required the
  // delta from the neighbor values to be positive definite or else it would fall back to side (and then edge)
  // computations.
  const std::array<stk::math::Vector3d,4> coords{{
    { 0.0000000000,  0.0000000000,  0.0000000000},
    { 0.3143942356, -0.7254582644, -0.5797545612},
    { 0.1419055462, -0.5553394556,  0.1665592194},
    {-0.1894032955, -0.4612201452, -0.4812267423}
  }};

  write_cubit_tet(coords);

  const std::array<double,3> distanceCase1 {{0.56866, 0.00317126,  0.24237}};
  const std::array<double,3> distanceCase2 {{0.56866, 0.00000000,  0.24237}};

  const double resultCase1 = eikonal_solve_tetrahedron(coords, distanceCase1, 1, 1.);
  const double resultCase2 = eikonal_solve_tetrahedron(coords, distanceCase2, 1, 1.);

  EXPECT_LE(resultCase2, resultCase1);
}

TEST(Eikonal_Calc, regularTet_nbrsAllZero_resultIsTetHeight)
{
  const std::array<double,3> distance{{0.,0.,0.}};

  const double edgeLen = get_regular_tetrahedron_edge_length();
  const double gold = std::sqrt(2./3.)*edgeLen;
  expect_tetrahedron_distance(get_regular_tetrahedron_coordinates(), distance, gold);
}

TEST(Eikonal_Calc, regularTet_nbrsAllNearZero_resultIsTetHeight)
{
  const double eps = 1.e-8;
  const std::array<double,3> distance{{0.,0.,eps}};

  const double edgeLen = get_regular_tetrahedron_edge_length();
  const double gold = std::sqrt(2./3.)*edgeLen;
  expect_tetrahedron_distance(get_regular_tetrahedron_coordinates(), distance, gold);
}

TEST(Eikonal_Calc, regularTet_nbrsAllZero_distanceIsHeightAndExtensionFromFaceCentroid)
{
  const std::array<double,3> distance{{0.,0.,0.}};
  const std::array<double,3> extSpeed{{0.,1.,2.}};

  const double edgeLen = get_regular_tetrahedron_edge_length();
  const double goldDist = std::sqrt(2./3.)*edgeLen;
  const double goldExtSpeed = (extSpeed[0]+extSpeed[1]+extSpeed[2])/3.;
  expect_tetrahedron_distance_and_extension_speed(get_regular_tetrahedron_coordinates(), distance, extSpeed, goldDist, goldExtSpeed);
}

TEST(Eikonal_Calc, regularTet_gradienAlong0to3_resultIsEdgeLength)
{
  const std::array<stk::math::Vector3d,4> coords = get_regular_tetrahedron_coordinates();
  const double edgeLen = get_regular_tetrahedron_edge_length();
  const std::array<double,3> distance{{0.,Dot(coords[1]-coords[0],coords[3]-coords[0])/edgeLen,Dot(coords[2]-coords[0],coords[3]-coords[0])/edgeLen}};

  const double gold = edgeLen;
  expect_tetrahedron_distance(coords, distance, gold);
}

TEST(Eikonal_Calc, regularTet_gradienAlong0to3_distanceIsEdgeLengthAndExtensionFromNode0)
{
  const std::array<stk::math::Vector3d,4> coords = get_regular_tetrahedron_coordinates();
  const double edgeLen = get_regular_tetrahedron_edge_length();
  const std::array<double,3> distance{{0.,Dot(coords[1]-coords[0],coords[3]-coords[0])/edgeLen,Dot(coords[2]-coords[0],coords[3]-coords[0])/edgeLen}};
  const std::array<double,3> extSpeed{{1.,2.,3.}};

  const double goldDist = edgeLen;
  const double goldExtSpeed = extSpeed[0];
  expect_tetrahedron_distance_and_extension_speed(get_regular_tetrahedron_coordinates(), distance, extSpeed, goldDist, goldExtSpeed);
}

TEST(Eikonal_Calc, regularTri_nbrsAllZero_resultIsTriHeight)
{
  const std::array<double,2> distance{{0.,0.}};

  const double gold = 0.5*std::sqrt(3.0)*get_regular_triangle_edge_length();
  expect_triangle_distance(get_regular_triangle_coordinates(), distance, gold);
}

TEST(Eikonal_Calc, regularTri_nbrsAllZero_distanceIsHeightAndExtensionFromEdgeMidpt)
{
  const std::array<double,2> distance{{0.,0.}};
  const std::array<double,2> extSpeed{{0.,1.}};

  const double goldDist = 0.5*std::sqrt(3.0)*get_regular_triangle_edge_length();
  const double goldExtSpeed = (extSpeed[0]+extSpeed[1])/2.;
  expect_triangle_distance_and_extension_speed(get_regular_triangle_coordinates(), distance, extSpeed, goldDist, goldExtSpeed);
}

TEST(Eikonal_Calc, regularTri_gradienAlong0to2_resultIsEdgeLength)
{
  const std::array<stk::math::Vector3d,3> coords = get_regular_triangle_coordinates();
  const double edgeLen = get_regular_triangle_edge_length();
  const std::array<double,2> distance{{0.,Dot(coords[1]-coords[0],coords[2]-coords[0])/edgeLen}};

  const double gold = edgeLen;
  expect_triangle_distance(coords, distance, gold);
}

TEST(Eikonal_Calc, regularTri_gradienAlong0to2_distanceIsEdgeLengthAndExtensionFromNode0)
{
  const std::array<stk::math::Vector3d,3> coords = get_regular_triangle_coordinates();
  const double edgeLen = get_regular_triangle_edge_length();
  const std::array<double,2> distance{{0.,Dot(coords[1]-coords[0],coords[2]-coords[0])/edgeLen}};
  const std::array<double,2> extSpeed{{1.,2.}};

  const double goldDist = edgeLen;
  const double goldExtSpeed = extSpeed[0];
  expect_triangle_distance_and_extension_speed(coords, distance, extSpeed, goldDist, goldExtSpeed);
}

TEST(Eikonal_Calc, rightTriangleWithOneNeighborNodeAsClosestPointAndOtherNbrOnLongerEdge_resultIsDistanceFromClosestPoint)
{
  const std::array<stk::math::Vector3d,3> coords{{ {0.,0.,0.}, {2.,0.,0.}, {0.,1.,0.} }};
  const std::array<double,2> distance{{0.,2.}};

  const double gold = 1.0;
  expect_triangle_distance(coords, distance, gold);
}

TEST(Eikonal_Calc, regularTriangle_withFrontEmanatingFromOutsideBase_resultIsEdgeLength)
{
  const std::array<stk::math::Vector3d,3> coords = get_regular_triangle_coordinates();
  const double gold = 1.0;

  {
    const std::array<double,2> distance{{0., 0.9}};
    expect_triangle_distance(coords, distance, gold);
  }

  {
    const std::array<double,2> distance{{0.9, 0.}};
    expect_triangle_distance(coords, distance, gold);
  }

  // detJ = 0 case
  {
    const std::array<double,2> distance{{0., 1.}};
    expect_triangle_distance(coords, distance, gold);
  }
}

TEST(Eikonal_Calc, regularTriangle_withFrontEmanatingFromOutsideBase_distanceIsEdgeLengthAndExtensionFromNearestNode)
{
  const std::array<stk::math::Vector3d,3> coords = get_regular_triangle_coordinates();
  const double goldDist = 1.0;
  const std::array<double,2> extSpeed{{1.,2.}};

  {
    const std::array<double,2> distance{{0., 0.9}};
    expect_triangle_distance_and_extension_speed(coords, distance, extSpeed, goldDist, extSpeed[0]);
  }

  {
    const std::array<double,2> distance{{0.9, 0.}};
    expect_triangle_distance_and_extension_speed(coords, distance, extSpeed, goldDist, extSpeed[1]);
  }

  // detJ = 0 case
  {
    const std::array<double,2> distance{{0., 1.}};
    expect_triangle_distance_and_extension_speed(coords, distance, extSpeed, goldDist, extSpeed[0]);
  }
}

TEST(Eikonal_Calc, rightTetrahedronWithOneNeighborNodeAsClosestPointAndAnotherNbrOnLongerEdge_resultIsDistanceFromClosestPoint)
{
  const std::array<stk::math::Vector3d,4> coords{{ {0.,0.,0.}, {2.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} }};
  const std::array<double,3> distance{{0.,2.,1.}};

  const double gold = 1.0;
  expect_tetrahedron_distance(coords, distance, gold);
}

TEST(Eikonal_Calc, rightTriangle_nbrsAllZeroAndNonUnitSpeed_resultIsTimeOfArrival)
{
  const std::array<double,2> distance{{0.,0.}};
  const double speed = 2.0;

  const double gold = 1.0/speed;
  expect_triangle_eikonal_solution(get_right_triangle_coordinates(), distance, speed, gold);
}

TEST(Eikonal_Calc, rightTetrahedron_nbrsAllZeroAndUnitSpeed_resultIsDistance)
{
  const std::array<double,3> distance{{0.,0.,0.}};
  const double speed = 1.0;

  const double gold = get_right_tetrahedron_edge_length()/speed;
  expect_tetrahedron_eikonal_solution(get_right_tetrahedron_coordinates(), distance, speed, gold);
}

TEST(Eikonal_Calc, rightTetrahedron_nbrsAllZeroAndNonUnitSpeed_resultIsTimeOfArrival)
{
  const std::array<double,3> distance{{0.,0.,0.}};
  const double speed = 2.0;

  const double gold = 1.0/speed;
  expect_tetrahedron_eikonal_solution(get_right_tetrahedron_coordinates(), distance, speed, gold);
}

template<size_t SIZE>
std::array<double,SIZE> compute_distance_to_plane(const stk::math::Plane3d & plane, const std::array<stk::math::Vector3d,SIZE> & coords)
{
  std::array<double,SIZE> result;
  for (size_t i=0; i<SIZE; ++i)
    result[i] = plane.signed_distance(coords[i]);
  return result;
}

void expect_regular_triangle_distance_correct_for_closest_point_on_opposite_side(const double s)
{
  const std::array<stk::math::Vector3d,3> coords = get_regular_triangle_coordinates();

  ASSERT_TRUE(s >= 0.);
  const stk::math::Vector3d origin = (1.-s)*coords[0] + s*coords[1];
  const int node = s < 0.5 ? 0 : 1;

  const stk::math::Vector3d & gradientDir = coords[2] - origin;
  const stk::math::Plane3d plane(gradientDir, coords[node]);
  const auto distanceToPlane = compute_distance_to_plane(plane, coords);
  const double gold = distanceToPlane[2];
  const std::array<double,2> distance = {{distanceToPlane[0], distanceToPlane[1]}};
  expect_triangle_distance(coords, distance, gold);
}

TEST(Eikonal_Calc, regularTri_variousGradients)
{
  expect_regular_triangle_distance_correct_for_closest_point_on_opposite_side(0.);
  expect_regular_triangle_distance_correct_for_closest_point_on_opposite_side(1.);
  expect_regular_triangle_distance_correct_for_closest_point_on_opposite_side(0.1);
  expect_regular_triangle_distance_correct_for_closest_point_on_opposite_side(0.5);
  expect_regular_triangle_distance_correct_for_closest_point_on_opposite_side(0.8);
}

void expect_regular_tetrahedron_distance_correct_for_closest_point_on_opposite_side(const double r, const double s)
{
  const std::array<stk::math::Vector3d,4> coords = get_regular_tetrahedron_coordinates();

  const double t = 1.-r-s;
  ASSERT_TRUE(r >= 0. && s >= 0. && t >= 0.);
  const int node = (t>r && t>s) ? 0 : ((r>s && r>t) ? 1 : 2);

  const stk::math::Vector3d origin = t*coords[0] + r*coords[1] + s*coords[2];

  const stk::math::Vector3d & gradientDir = coords[3] - origin;
  const stk::math::Plane3d plane(gradientDir, coords[node]);
  const auto distanceToPlane = compute_distance_to_plane(plane, coords);
  const double gold = distanceToPlane[3];
  const std::array<double,3> distance = {{distanceToPlane[0], distanceToPlane[1], distanceToPlane[2]}};
  expect_tetrahedron_distance(coords, distance, gold);
}

void expect_regular_tetrahedron_distance_correct_for_closest_point_outside_of_opposite_side(const int edge, const double sEdge)
{
  const std::array<stk::math::Vector3d,4> coords = get_regular_tetrahedron_coordinates();

  ASSERT_TRUE(edge>=0 && edge<=2);
  ASSERT_TRUE(sEdge>=0. && sEdge<=1.);

  const double rFace = (edge==0) ? sEdge : (edge==1) ? (1.-sEdge) : 0.;
  const double sFace = (edge==0) ? 0 : (edge==1) ? sEdge : (1.-sEdge);
  const double tFace = 1.-rFace-sFace;
  const int node = (tFace>rFace && tFace>sFace) ? 0 : ((rFace>sFace && rFace>tFace) ? 1 : 2);

  const stk::math::Vector3d originBaseFace = tFace*coords[0] + rFace*coords[1] + sFace*coords[2];

  const stk::math::Vector3d & gradientDir = coords[3] - originBaseFace;
  const stk::math::Plane3d plane(gradientDir, coords[node]);
  const auto distanceToPlane = compute_distance_to_plane(plane, coords);
  std::array<double,3> distance = {{distanceToPlane[0], distanceToPlane[1], distanceToPlane[2]}};

  // Move front outside of base triangle
  const double edgeLen = get_regular_tetrahedron_edge_length();
  const std::array<int,3> edgeOppositeNode{{ 2, 0, 1 }};
  distance[edgeOppositeNode[edge]] += 0.1*edgeLen;

  const double gold = distanceToPlane[3];
  expect_tetrahedron_distance(coords, distance, gold);
}

TEST(Eikonal_Calc, regularTet_variousGradients)
{
  expect_regular_tetrahedron_distance_correct_for_closest_point_on_opposite_side(0.,0.);
  expect_regular_tetrahedron_distance_correct_for_closest_point_on_opposite_side(1.,0.);
  expect_regular_tetrahedron_distance_correct_for_closest_point_on_opposite_side(0.,1.);
  expect_regular_tetrahedron_distance_correct_for_closest_point_on_opposite_side(0.333,0.333);
  expect_regular_tetrahedron_distance_correct_for_closest_point_on_opposite_side(0.25,0.7);

  expect_regular_tetrahedron_distance_correct_for_closest_point_outside_of_opposite_side(0, 0.25);
  expect_regular_tetrahedron_distance_correct_for_closest_point_outside_of_opposite_side(0, 0.55);
  expect_regular_tetrahedron_distance_correct_for_closest_point_outside_of_opposite_side(1, 0.1);
  expect_regular_tetrahedron_distance_correct_for_closest_point_outside_of_opposite_side(1, 0.7);
  expect_regular_tetrahedron_distance_correct_for_closest_point_outside_of_opposite_side(2, 0.3);
  expect_regular_tetrahedron_distance_correct_for_closest_point_outside_of_opposite_side(2, 0.8);
}

TEST(Eikonal_Calc, regularTetrahedron_withFrontEmanatingFromOutsideBase_resultIsEdgeLength)
{
  const std::array<stk::math::Vector3d,4> coords = get_regular_tetrahedron_coordinates();
  const double edgeLen = get_regular_tetrahedron_edge_length();
  const double gold = edgeLen;

  {
    const std::array<double,3> distance{{0., 0.6*edgeLen, 0.6*edgeLen}};
    expect_tetrahedron_distance(coords, distance, gold);
  }

  {
    const std::array<double,3> distance{{0.6*edgeLen, 0., 0.6*edgeLen}};
    expect_tetrahedron_distance(coords, distance, gold);
  }

  {
    const std::array<double,3> distance{{0.6*edgeLen, 0.6*edgeLen, 0.}};
    expect_tetrahedron_distance(coords, distance, gold);
  }
}

TEST(Eikonal_Calc, regularTetrahedron_withFrontEmanatingFromOutsideBase_distanceIsEdgeLengthAndExtensionFromNearestNode)
{
  const std::array<stk::math::Vector3d,4> coords = get_regular_tetrahedron_coordinates();
  const double edgeLen = get_regular_tetrahedron_edge_length();
  const double goldDist = edgeLen;
  const std::array<double,3> extSpeed{{1.,2.,3.}};

  {
    const std::array<double,3> distance{{0., 0.6*edgeLen, 0.6*edgeLen}};
    expect_tetrahedron_distance_and_extension_speed(coords, distance, extSpeed, goldDist, extSpeed[0]);
  }

  {
    const std::array<double,3> distance{{0.6*edgeLen, 0., 0.6*edgeLen}};
    expect_tetrahedron_distance_and_extension_speed(coords, distance, extSpeed, goldDist, extSpeed[1]);
  }

  {
    const std::array<double,3> distance{{0.6*edgeLen, 0.6*edgeLen, 0.}};
    expect_tetrahedron_distance_and_extension_speed(coords, distance, extSpeed, goldDist, extSpeed[2]);
  }
}

TEST(Eikonal_Calc, rightTriangle_usingMultipleSpeeds_resultIsLengthOverSpeed)
{
  const std::array<double,2> distance{{0.,0.}};

  {
    const double speed = 0.5;
    const double gold = get_right_triangle_edge_length()/speed;
    expect_triangle_eikonal_solution(get_right_triangle_coordinates(),distance,speed,gold);
  }

  {
    const double speed = 2.0;
    const double gold = get_right_triangle_edge_length()/speed;
    expect_triangle_eikonal_solution(get_right_triangle_coordinates(),distance,speed,gold);
  }
}

TEST(Eikonal_Calc, rightTetrahedron_usingMultipleSpeeds_resultIsLengthOverSpeed)
{
  const std::array<double,3> distance{{0.,0.,0.}};

  {
    const double speed = 0.5;
    const double gold = get_right_tetrahedron_edge_length()/speed;
    expect_tetrahedron_eikonal_solution(get_right_tetrahedron_coordinates(),distance,speed,gold);
  }

  {
    const double speed = 2.0;
    const double gold = get_right_tetrahedron_edge_length()/speed;
    expect_tetrahedron_eikonal_solution(get_right_tetrahedron_coordinates(),distance,speed,gold);
  }
}

TEST(Eikonal_Calc, rightTriangle_negativeSign_resultIsNegativeLength)
{
  const std::array<double,2> distance{{0.,0.}};

  const double speed = 1.0;
  const int sign = -1;
  const double gold = get_right_triangle_edge_length()/speed*sign;
  const double resultCase = eikonal_solve_triangle(get_right_triangle_coordinates(), distance, sign, farDistance, speed);
  EXPECT_NEAR(gold, resultCase, 1.e-8);
}

TEST(Eikonal_Calc, rightTetrahedron_negativeSign_resultIsNegativeLength)
{
  const std::array<double,3> distance{{0.,0.,0.}};

  const double speed = 1.0;
  const int sign = -1;
  const double gold = get_right_tetrahedron_edge_length()/speed*sign;
  const double resultCase = eikonal_solve_tetrahedron(get_right_tetrahedron_coordinates(), distance, sign, farDistance, speed);
  EXPECT_NEAR(gold, resultCase, 1.e-8);
}

TEST(Eikonal_Calc, rightTriangle_upstreamNodeIsFar_resultIsFar)
{
  const double far = 100.0;
  const double speed = 1.0;

  const std::array<double,2> distance{{0.,far}};
  const double gold = far;
  expect_triangle_eikonal_solution(get_right_triangle_coordinates(),distance,speed,gold,far);
}

TEST(Eikonal_Calc, rightTetrahedron_upstreamNodeIsFar_resultIsFar)
{
  const double far = 100.0;
  const double speed = 1.0;

  const std::array<double,3> distance{{0.,far, 0.}};
  const double gold = far;
  expect_tetrahedron_eikonal_solution(get_right_tetrahedron_coordinates(),distance,speed,gold,far);
}

} //krino
