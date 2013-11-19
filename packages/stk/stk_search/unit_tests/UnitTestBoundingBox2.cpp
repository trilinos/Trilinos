#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_search/BoundingBox2.hpp>

#include <string>
#include <iostream>

namespace {


STKUNIT_UNIT_TEST( stk_search, Point3)
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


STKUNIT_UNIT_TEST( stk_search, Sphere)
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

  STKUNIT_EXPECT_TRUE( s.is_valid());
  STKUNIT_EXPECT_EQ(s, Sphere(Point(1,2,3),4));
  STKUNIT_EXPECT_NE(s, Sphere());
}

} // unnamed namespace
