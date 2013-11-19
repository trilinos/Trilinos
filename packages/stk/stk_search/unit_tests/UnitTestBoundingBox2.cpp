#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_search/BoundingBox2.hpp>

#include <string>
#include <iostream>

namespace {


STKUNIT_UNIT_TEST( stk_search_bounding_box, Point3)
{
  using namespace stk::search;

  Point3<double> p;

  {
    std::ostringstream out;
    out << p;
    STKUNIT_EXPECT_EQUAL( out.str(), std::string("(0,0,0)"));
  }

  p = Point3<double>(1,2,3);

  STKUNIT_EXPECT_EQ(p, Point3<double>(1,2,3));
  STKUNIT_EXPECT_NE(p, Point3<double>());
}


STKUNIT_UNIT_TEST( stk_search_bounding_box, Sphere)
{
  typedef stk::search::Sphere<double> Sphere;
  typedef stk::search::Point3<double> Point;

  Sphere s;
  STKUNIT_EXPECT_FALSE( s.is_valid());

  {
    std::ostringstream out;
    out << s;
    STKUNIT_EXPECT_EQUAL( out.str(), std::string("{(0,0,0):-1}"));
  }

  s = Sphere(Point(1,2,3),4);

  STKUNIT_EXPECT_EQ( s.center(), Point(1,2,3) );
  STKUNIT_EXPECT_EQ( s.radius(), 4 );

  STKUNIT_EXPECT_TRUE( s.is_valid());
  STKUNIT_EXPECT_EQ(s, Sphere(Point(1,2,3),4));
  STKUNIT_EXPECT_NE(s, Sphere());
}


STKUNIT_UNIT_TEST( stk_search_bounding_box, AABox)
{
  typedef stk::search::AABox<double> Box;
  typedef stk::search::Point3<double> Point;

  Box b;
  STKUNIT_EXPECT_FALSE( b.is_valid());

  {
    std::ostringstream out;
    out << b;
    STKUNIT_EXPECT_EQUAL( out.str(), std::string("{(1,1,1)->(0,0,0)}"));
  }

  b = Box(Point(0,0,0),Point(1,1,1));

  STKUNIT_EXPECT_EQ( b.min_corner(), Point(0,0,0) );
  STKUNIT_EXPECT_EQ( b.max_corner(), Point(1,1,1) );

  STKUNIT_EXPECT_TRUE( b.is_valid());
  STKUNIT_EXPECT_EQ(b, Box(Point(0,0,0),Point(1,1,1)));
  STKUNIT_EXPECT_NE(b, Box());
}

STKUNIT_UNIT_TEST( stk_search_bounding_box, min_max_center)
{
  typedef stk::search::Sphere<double> Sphere;
  typedef stk::search::AABox<double> Box;
  typedef stk::search::Point3<double> Point;

  using stk::search::min_corner;
  using stk::search::max_corner;
  using stk::search::center;

  Point p(0,0,0);
  STKUNIT_EXPECT_EQ( Point(0,0,0), min_corner(p) );
  STKUNIT_EXPECT_EQ( Point(0,0,0), max_corner(p) );
  STKUNIT_EXPECT_EQ( Point(0,0,0), center(p) );

  Box b(Point(0,0,0),Point(1,1,1));
  STKUNIT_EXPECT_EQ( Point(0,0,0), min_corner(b) );
  STKUNIT_EXPECT_EQ( Point(1,1,1), max_corner(b) );
  STKUNIT_EXPECT_EQ( Point(0.5,0.5,0.5), center(b) );

  Sphere s(Point(1,1,1),1.5);
  STKUNIT_EXPECT_EQ( Point(-0.5, -0.5, -0.5), min_corner(s) );
  STKUNIT_EXPECT_EQ( Point(2.5, 2.5, 2.5), max_corner(s) );
  STKUNIT_EXPECT_EQ( Point(1, 1, 1), center(s) );

}

STKUNIT_UNIT_TEST( stk_search_bounding_box, intersects)
{
  typedef stk::search::Sphere<double> Sphere;
  typedef stk::search::AABox<double> Box;
  typedef stk::search::Point3<double> Point;

  using stk::search::intersects;

  {
    Point a(0,0,0);
    Point b(0,0,0);
    STKUNIT_EXPECT_TRUE( intersects(a,b) );
    a[0] = 1;
    STKUNIT_EXPECT_FALSE( intersects(a,b) );
  }

  {
    Point a(0,0,0);
    Sphere b(Point(0,0,0),1);
    STKUNIT_EXPECT_TRUE( intersects(a,b) );
    STKUNIT_EXPECT_TRUE( intersects(b,a) );
    a[0] = 2;
    STKUNIT_EXPECT_FALSE( intersects(a,b) );
    STKUNIT_EXPECT_FALSE( intersects(b,a) );
  }

  {
    Point a(0.5,0.5,0.5);
    Box b(Point(0,0,0),Point(1,1,1));
    STKUNIT_EXPECT_TRUE( intersects(a,b) );
    STKUNIT_EXPECT_TRUE( intersects(b,a) );
    a[0] = 2;
    STKUNIT_EXPECT_FALSE( intersects(a,b) );
    STKUNIT_EXPECT_FALSE( intersects(b,a) );
  }


}


} // unnamed namespace
