// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include <gtest/gtest.h>

#include <stk_search/BoundingBox.hpp>
#include <stk_math/stk_math/StkVector.hpp>

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

TEST( stk_search_bounding_box, intersects)
{
  typedef stk::search::Sphere<double> Sphere;
  typedef stk::search::Box<double> Box;
  typedef stk::search::Point<double> Point;

  using stk::search::intersects;

  //Point,Point
  {
    Point a(0,0,0);
    Point b(0,0,0);
    EXPECT_TRUE( intersects(a,b) );
    a[0] = 1;
    EXPECT_FALSE( intersects(a,b) );
  }

  //Point,Sphere
  {
    Point a(0,0,0);
    Sphere b(Point(0,0,0),1);
    EXPECT_TRUE( intersects(a,b) );
    EXPECT_TRUE( intersects(b,a) );
    a[0] = 2;
    EXPECT_FALSE( intersects(a,b) );
    EXPECT_FALSE( intersects(b,a) );
  }

  //Point,Box
  {
    Point a(0.5,0.5,0.5);
    Box b(Point(0,0,0),Point(1,1,1));
    EXPECT_TRUE( intersects(a,b) );
    EXPECT_TRUE( intersects(b,a) );
    a[0] = 2;
    EXPECT_FALSE( intersects(a,b) );
    EXPECT_FALSE( intersects(b,a) );
  }

  //Sphere,Sphere
  {
    Sphere a(Point(0,0,0),2);
    Sphere b(Point(1,1,1),2);
    EXPECT_TRUE( intersects(a,b) );
    EXPECT_TRUE( intersects(b,a) );
    a.set_center(Point(-3,-3,-3));
    EXPECT_FALSE( intersects(a,b) );
    EXPECT_FALSE( intersects(b,a) );
  }

  //Sphere,Box
  {
    Sphere a(Point(0,0,0),1);
    Box b(Point(.3,.3,.3),Point(2,2,2));
    EXPECT_TRUE( intersects(a,b) );
    EXPECT_TRUE( intersects(b,a) );
    b.set_min_corner(Point(.9,.9,.9));
    EXPECT_FALSE( intersects(a,b) );
    EXPECT_FALSE( intersects(b,a) );
  }

  //Box,Box
  {
    Box a(Point(0,0,0),Point(1.5,1.5,1.5));
    Box b(Point(1,1,1),Point(3,3,3));
    EXPECT_TRUE( intersects(a,b) );
    EXPECT_TRUE( intersects(b,a) );
    b.set_min_corner(Point(2,1,1));
    EXPECT_FALSE( intersects(a,b) );
    EXPECT_FALSE( intersects(b,a) );
  }
}


} // unnamed namespace
