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

#include <gtest/gtest.h>

#include <stk_math/StkVector.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/Plane.hpp>

#include <string>
#include <iostream>
#include <limits>

namespace {


TEST( stk_search_bounding_box, Point)
{
  using namespace stk::search;

  Point<double> p;

  {
    std::ostringstream out;
    out << p;
    EXPECT_EQ( out.str(), std::string("(0,0,0)"));
  }

  p = Point<double>(1,2,3);

  EXPECT_EQ(p, Point<double>(1,2,3));
  EXPECT_NE(p, Point<double>());
}


TEST( stk_search_bounding_box, Sphere)
{
  typedef stk::search::Sphere<double> Sphere;
  typedef stk::search::Point<double> Point;

  using stk::search::is_valid;

  Sphere s;
  EXPECT_FALSE( is_valid(s));

  {
    std::ostringstream out;
    out << s;
    EXPECT_EQ( out.str(), std::string("{(0,0,0):-1}"));
  }

  s = Sphere(Point(1,2,3),4);

  EXPECT_EQ( s.center(), Point(1,2,3) );
  EXPECT_EQ( s.radius(), 4 );

  EXPECT_TRUE( is_valid(s));
  EXPECT_EQ(s, Sphere(Point(1,2,3),4));
  EXPECT_NE(s, Sphere());
}

typedef stk::search::Sphere<double> Sphere;
typedef stk::search::Box<double> Box;
typedef stk::search::Point<double> Point;
typedef stk::math::Vector3d Vec;

bool is_coplanar(const Point& point, const std::vector<Point>& tri)
{
    Vec p = Vec(point[0], point[1], point[2]);
    Vec a = Vec(tri[0][0], tri[0][1], tri[0][2]);
    Vec b = Vec(tri[1][0], tri[1][1], tri[1][2]);
    Vec c = Vec(tri[2][0], tri[2][1], tri[2][2]);

    Vec normal = Cross(b - a, c - a);
    return (Dot(a - p, normal) <= std::numeric_limits<double>::epsilon());
}

double det(const Vec& a, const Vec& b)
{
    Vec c = Cross(a,b);
    return c.length();
}

bool is_interior_point(const Point& point, const std::vector<Point>& tri)
{
    Vec p = Vec(point[0], point[1], point[2]);
    Vec a = Vec(tri[0][0], tri[0][1], tri[0][2]);
    Vec b = Vec(tri[1][0], tri[1][1], tri[1][2]);
    Vec c = Vec(tri[2][0], tri[2][1], tri[2][2]);

    Vec v0 = a;
    Vec v1 = c - a;
    Vec v2 = b - a;

    double denom = det(v1, v2);
    ThrowRequireMsg(denom > 0, "denominator cannot be zero");
    double alpha =  (det(p, v2) - det(v0, v2)) / denom;
    double beta  = -(det(p, v1) - det(v0, v1)) / denom;
    return (alpha >= 0 && beta >= 0 && alpha + beta <= 1);
}

size_t compute_num_equidistant_points(const Point& point, const std::vector<Point>& tri)
{
    Vec p = Vec(point[0], point[1], point[2]);
    Vec a = Vec(tri[0][0], tri[0][1], tri[0][2]);
    Vec b = Vec(tri[1][0], tri[1][1], tri[1][2]);
    Vec c = Vec(tri[2][0], tri[2][1], tri[2][2]);

    Vec PA = a - p;
    Vec PB = b - p;
    Vec PC = c - p;

    double lengthPA = PA.length();
    double lengthPB = PB.length();
    double lengthPC = PC.length();

    size_t numEquidistantPoints = 0;
    if (lengthPA == lengthPB && lengthPA == lengthPC)
        numEquidistantPoints = 3;
    else if(lengthPA == lengthPB || lengthPA == lengthPC || lengthPB == lengthPC)
        numEquidistantPoints = 2;
    return numEquidistantPoints ;
}

TEST( stk_search_bounding_box, minimumSphereEnclosingEquilateralTriangle)
{
  using stk::search::is_valid;

  std::vector<Point> tri = {Point(-1,0,0),
                            Point(1,0,0),
                            Point(0,std::sqrt(3),0)};

  Sphere s = minimumBoundingSphere(tri[0], tri[1], tri[2]);
  EXPECT_TRUE(is_valid(s));
  EXPECT_TRUE(is_coplanar(s.center(), tri));
  EXPECT_TRUE(is_interior_point(s.center(), tri));

  size_t goldEquidistantPoints = 3;
  EXPECT_EQ(goldEquidistantPoints, compute_num_equidistant_points(s.center(), tri));
}

TEST( stk_search_bounding_box, minimumSphereEnclosingRightIsoscelesTriangle)
{
  using stk::search::is_valid;

  std::vector<Point> tri = {Point(1,0,0),
                            Point(0,1,0),
                            Point(0,0,0)};

  Sphere s = minimumBoundingSphere(tri[0], tri[1], tri[2]);
  EXPECT_TRUE(is_valid(s));
  EXPECT_TRUE(is_coplanar(s.center(), tri));

  size_t goldEquidistantPoints = 3;
  EXPECT_EQ(goldEquidistantPoints, compute_num_equidistant_points(s.center(), tri));
}

TEST( stk_search_bounding_box, minimumSphereEnclosingObtuseIsoscelesTriangle)
{
  using stk::search::is_valid;

  std::vector<Point> tri = {Point(-1,0,0),
                            Point(1,0,0),
                            Point(0,0.1,0)};

  Sphere s = minimumBoundingSphere(tri[0], tri[1], tri[2]);
  EXPECT_TRUE(is_valid(s));
  EXPECT_TRUE(is_coplanar(s.center(), tri));
  EXPECT_TRUE(is_interior_point(s.center(), tri));

  size_t goldEquidistantPoints = 2;
  EXPECT_EQ(goldEquidistantPoints, compute_num_equidistant_points(s.center(), tri));
}

TEST( stk_search_bounding_box, minimumSphereEnclosingLine)
{
  using stk::search::is_valid;

  std::vector<Point> tri = {Point(0,0,0),
                            Point(1,1,1),
                            Point(3,3,3)};

  Sphere s = minimumBoundingSphere(tri[0], tri[1], tri[2]);
  EXPECT_TRUE(is_valid(s));
  EXPECT_TRUE(is_coplanar(s.center(), tri));

  size_t goldEquidistantPoints = 2;
  EXPECT_EQ(goldEquidistantPoints, compute_num_equidistant_points(s.center(), tri));
}

TEST( stk_search_bounding_box, Box)
{
  typedef stk::search::Box<double> Box;
  typedef stk::search::Point<double> Point;

  using stk::search::is_valid;

  Box b;
  EXPECT_FALSE( is_valid(b));

  {
    std::ostringstream out;
    out << b;
    EXPECT_EQ( out.str(), std::string("{(1.79769e+308,1.79769e+308,1.79769e+308)->(-1.79769e+308,-1.79769e+308,-1.79769e+308)}"));
  }

  b = Box(Point(0,0,0),Point(1,1,1));

  EXPECT_EQ( b.min_corner(), Point(0,0,0) );
  EXPECT_EQ( b.max_corner(), Point(1,1,1) );

  EXPECT_TRUE( is_valid(b));
  EXPECT_EQ(b, Box(Point(0,0,0),Point(1,1,1)));
  EXPECT_NE(b, Box());
}

TEST( stk_search_bounding_box, min_max_center)
{
  typedef stk::search::Sphere<double> Sphere;
  typedef stk::search::Box<double> Box;
  typedef stk::search::Point<double> Point;

  using stk::search::min_corner;
  using stk::search::max_corner;
  using stk::search::center;

  Point p(0,0,0);
  EXPECT_EQ( Point(0,0,0), min_corner(p) );
  EXPECT_EQ( Point(0,0,0), max_corner(p) );
  EXPECT_EQ( Point(0,0,0), center(p) );

  Box b(Point(0,0,0),Point(1,1,1));
  EXPECT_EQ( Point(0,0,0), min_corner(b) );
  EXPECT_EQ( Point(1,1,1), max_corner(b) );
  EXPECT_EQ( Point(0.5,0.5,0.5), center(b) );

  Sphere s(Point(1,1,1),1.5);
  EXPECT_EQ( Point(-0.5, -0.5, -0.5), min_corner(s) );
  EXPECT_EQ( Point(2.5, 2.5, 2.5), max_corner(s) );
  EXPECT_EQ( Point(1, 1, 1), center(s) );

}

template <typename T1, typename T2>
void CheckIntersections(T1 obj1, T2 obj2, bool expectedResult) {
    using stk::search::intersects;
    EXPECT_EQ(intersects(obj1, obj2), expectedResult);
    EXPECT_EQ(intersects(obj2, obj1), expectedResult);
}


TEST( stk_search_bounding_box, intersects_point_point)
{
  typedef stk::search::Point<double> Point;


  Point a(0,0,0);
  Point b(0,0,0);
  Point c(1,0,0);
  CheckIntersections(a, b, true);
  CheckIntersections(a, c, false);
}


TEST( stk_search_bounding_box, intersects_point_sphere)
{
  typedef stk::search::Sphere<double> Sphere;
  typedef stk::search::Point<double> Point;

  Point a(0,0,0);
  Sphere b(Point(0,0,0),1);
  Point c(2, 0, 0);

  CheckIntersections(a, b, true);
  CheckIntersections(a, c, false);
}

TEST( stk_search_bounding_box, intersects_point_box)
{
  typedef stk::search::Box<double> Box;
  typedef stk::search::Point<double> Point;


  Point a(0.5,0.5,0.5);
  Box b(Point(0,0,0),Point(1,1,1));
  Point c(2, 0, 0);

  CheckIntersections(a, b, true);
  CheckIntersections(a, c, false);
}


TEST( stk_search_bounding_box, intersects_sphere_sphere)
{
  typedef stk::search::Sphere<double> Sphere;

  using stk::search::intersects;

  Sphere a(Point(0,0,0),2);
  Sphere b(Point(1,1,1),2);
  Sphere c(Point(-3, -3, -3), 2);

  CheckIntersections(a, b, true);
  CheckIntersections(a, c, false);

}


TEST( stk_search_bounding_box, intersects_sphere_box)
{
  typedef stk::search::Sphere<double> Sphere;
  typedef stk::search::Box<double> Box;
  using stk::search::intersects;

  Box unitBox(Point(1, 1, 1), Point(2, 2, 2));
  Box nullSetBox;

  Sphere wayLeft           (Point(-5, 1.5, 1.5), 0.5);
  Sphere wayRight          (Point( 5, 1.5, 1.5), 0.5);
  Sphere diagonallyOutside1(Point(0.5, 0.5, 0.5), 0.8);
  Sphere diagonallyInside1 (Point(0.5, 0.5, 0.5), 0.9);
  Sphere diagonallyOutside2(Point(2.5, 2.5, 2.5), 0.8);
  Sphere diagonallyInside2 (Point(2.5, 2.5, 2.5), 0.9);
  Sphere justOutsideLeft   (Point(0.5, 1.5, 1.5), 0.4999999);
  Sphere justTouchingLeft  (Point(0.5, 1.5, 1.5), 0.5);
  Sphere justTouchingRight (Point(2.5, 1.5, 1.5), 0.5);
  Sphere justOutsideRight  (Point(2.5, 1.5, 1.5), 0.4999999);
  Sphere encompassing      (Point(4, 3, 1), 20);
  Sphere inside            (Point(1.1, 1.2, 1.3), 0.01); 

  CheckIntersections(unitBox, wayLeft,            false);
  CheckIntersections(unitBox, wayRight,           false);
  CheckIntersections(unitBox, diagonallyOutside1, false);
  CheckIntersections(unitBox, diagonallyInside1,  true);
  CheckIntersections(unitBox, diagonallyOutside2, false);
  CheckIntersections(unitBox, diagonallyInside2,  true);
  CheckIntersections(unitBox, justOutsideLeft,    false);
  CheckIntersections(unitBox, justTouchingLeft,   true);
  CheckIntersections(unitBox, justTouchingRight,  true);
  CheckIntersections(unitBox, justOutsideRight,   false);
  CheckIntersections(unitBox, encompassing,       true);
  CheckIntersections(unitBox, inside,             true);

  CheckIntersections(nullSetBox, wayLeft,            false);
  CheckIntersections(nullSetBox, wayRight,           false);
  CheckIntersections(nullSetBox, diagonallyOutside1, false);
  CheckIntersections(nullSetBox, diagonallyInside1,  false);
  CheckIntersections(nullSetBox, diagonallyOutside2, false);
  CheckIntersections(nullSetBox, diagonallyInside2,  false);
  CheckIntersections(nullSetBox, justOutsideLeft,    false);
  CheckIntersections(nullSetBox, justTouchingLeft,   false);
  CheckIntersections(nullSetBox, justTouchingRight,  false);
  CheckIntersections(nullSetBox, justOutsideRight,   false);
  CheckIntersections(nullSetBox, encompassing,       false);
  CheckIntersections(nullSetBox, inside,             false);

}

TEST( stk_search_bounding_box, intersects_box_box)
{
  typedef stk::search::Box<double> Box;
  using stk::search::intersects;

  Box unitBox(Point(1, 1, 1), Point(2, 2, 2));
  Box nullSetBox;

  Box wayLeft          (Point(-5, 1, 1), Point(-4, 2, 2));
  Box wayRight         (Point( 5, 1, 1), Point( 6, 2, 2));
  Box diagonallyOutside(Point(0.5, 0.5, 0.5), Point(0.9, 0.9, 0.9));
  Box diagonallyInside (Point(0.5, 0.5, 0.5), Point(1.1, 1.1, 1.1));
  Box justOutsideLeft  (Point(0, 1, 1), Point(0.999999, 2, 2));
  Box justTouchingLeft (Point(0, 1, 1), Point(1, 2, 2));
  Box justTouchingRight(Point(2, 1, 1), Point(3, 2, 2));
  Box justOutsideRight (Point(2.000001, 1, 1), Point(3, 2, 2));
  Box encompassing     (Point(-7, -8, -9), Point(10, 11, 12));
  Box inside           (Point(1.1, 1.2, 1.3), Point(1.4, 1.5, 1.6)); 
  Box alsoNull;

  CheckIntersections(unitBox, wayLeft,           false);
  CheckIntersections(unitBox, wayRight,          false);
  CheckIntersections(unitBox, diagonallyOutside, false);
  CheckIntersections(unitBox, diagonallyInside,  true);
  CheckIntersections(unitBox, justOutsideLeft,   false);
  CheckIntersections(unitBox, justTouchingLeft,  true);
  CheckIntersections(unitBox, justTouchingRight, true);
  CheckIntersections(unitBox, justOutsideRight,  false);
  CheckIntersections(unitBox, encompassing,      true);
  CheckIntersections(unitBox, inside,            true);
  CheckIntersections(unitBox, alsoNull,          false);

  CheckIntersections(nullSetBox, wayLeft,           false);
  CheckIntersections(nullSetBox, wayRight,          false);
  CheckIntersections(nullSetBox, diagonallyOutside, false);
  CheckIntersections(nullSetBox, diagonallyInside,  false);
  CheckIntersections(nullSetBox, justOutsideLeft,   false);
  CheckIntersections(nullSetBox, justTouchingLeft,  false);
  CheckIntersections(nullSetBox, justTouchingRight, false);
  CheckIntersections(nullSetBox, justOutsideRight,  false);
  CheckIntersections(nullSetBox, encompassing,      false);
  CheckIntersections(nullSetBox, inside,            false);
  CheckIntersections(nullSetBox, alsoNull,          false);
}

TEST( stk_search_bounding_box, read_box_as_written) {
  std::stringstream ss;
  stk::search::Box<float> box({0,0,0}, {1,1,1});
  stk::search::Box<float> readBox, readBox2;
  ss << box << "\n" << box << "\n";
  ss >> readBox >> readBox2;

  EXPECT_EQ(box, readBox );
  EXPECT_EQ(box, readBox2 );
}

TEST( stk_search_bounding_box, intersects_box_plane)
{
  typedef stk::search::Box<double> Box;
  typedef stk::search::Plane<double> Plane;
  using stk::search::intersects;

  Box unitBox(Point(-0.5, -0.5, -0.5), Point(0.5, 0.5, 0.5));

  Plane nullPlane;
  Plane wayLeft           (Point(-2, 0.0, 0.0), Point(1.0, 0.0, 0.0));
  Plane wayRight          (Point( 5, 0.0, 0.0), Point( 1.0, 0.0, 0.0));
  Plane inside            (Point(0.15, 100.0, 0.0 ), Point(1.0, 0.0, 0.0));
  Plane inside2           (Point(0.15, 0.1, 0.0 ), Point(1.0, 0.0, 0.0));
  Plane justOutsideLeft   (Point(-0.500001, 1, 1), Point(-1.0, 0.0, 0.0));
  Plane justOutsideRight  (Point(0.500001, 1, 1), Point(-1.0, 0.0, 0.0));
  Plane justTouchingLeft  (Point(-0.5, 100.0, -99), Point(-1.0, 0.0, 0.0));
  Plane justTouchingRight (Point(0.5, 1, 0.36), Point(1.0, 0.0, 0.0));
  Plane diagonallyOutside (Point(-0.530330, 0.530330, 0.000000), Point(-0.707107, 0.707107, 0.000000));
  Plane diagonallyInside  (Point(-0.459619, 0.459619, 0.000000), Point(-0.707107, 0.707107, 0.000000));
  Plane cutCorner1        (Point( -0.325000, 0.459619, -0.325000), Point(-0.500000, 0.707107, -0.500000));
  Plane cutCorner2        (Point( -0.325000, 0.459619, 0.325000), Point(-0.500000, 0.707107, 0.500000));
  Plane cutCorner3        (Point(  0.325000, 0.459619, 0.325000), Point( 0.500000, 0.707107, 0.500000));
  Plane cutCorner4        (Point(  0.325000, 0.459619, -0.325000), Point( 0.500000, 0.707107, -0.500000));
  Plane cutCorner5        (Point( -0.459619, -0.325000, -0.325000), Point(-0.707107, -0.500000, -0.500000));
  Plane cutCorner6        (Point( -0.325000, -0.325000, 0.459619), Point(-0.500000, -0.500000, 0.707107));
  Plane cutCorner7        (Point( 0.459619, -0.325000, 0.325000), Point(0.707107, -0.500000, 0.500000));
  Plane cutCorner8        (Point( 0.325000, -0.325000, -0.459619), Point(0.500000, -0.500000, -0.707107));


  CheckIntersections(nullPlane,        unitBox, false);
  CheckIntersections(wayLeft,          unitBox, false);
  CheckIntersections(wayRight,         unitBox, false);
  CheckIntersections(inside,           unitBox, true);
  CheckIntersections(inside2,          unitBox, true);
  CheckIntersections(justOutsideLeft,  unitBox, false);
  CheckIntersections(justOutsideRight, unitBox, false);
  CheckIntersections(justTouchingLeft, unitBox, true);
  CheckIntersections(justTouchingRight,unitBox, true);
  CheckIntersections(diagonallyOutside,unitBox, false);
  CheckIntersections(diagonallyInside, unitBox, true);
  CheckIntersections(cutCorner1,       unitBox, true);
  CheckIntersections(cutCorner2,       unitBox, true);
  CheckIntersections(cutCorner3,       unitBox, true);
  CheckIntersections(cutCorner4,       unitBox, true);
  CheckIntersections(cutCorner5,       unitBox, true);
  CheckIntersections(cutCorner6,       unitBox, true);
  CheckIntersections(cutCorner7,       unitBox, true);
  CheckIntersections(cutCorner8,       unitBox, true);

}



} // unnamed namespace
