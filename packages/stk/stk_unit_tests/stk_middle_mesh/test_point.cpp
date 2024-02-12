#include <complex>
#include "gtest/gtest.h"

#include "stk_middle_mesh/point.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

using Complex = std::complex<double>;
using Point = utils::Point;
using PointC = utils::PointT<Complex>;

void test_float_eq(const Point& p1, const Point& p2)
{
  EXPECT_FLOAT_EQ(p1.x, p2.x);
  EXPECT_FLOAT_EQ(p1.y, p2.y);
  EXPECT_FLOAT_EQ(p1.z, p2.z);
}

void test_float_eq(const PointC& p1, const PointC& p2)
{
  for (int i=0; i < 3; ++i)
  {
    EXPECT_FLOAT_EQ(p1[i].real(), p2[i].real());
    EXPECT_FLOAT_EQ(p1[i].imag(), p2[i].imag());
  }
}

void expect_near(const Complex& a, const Complex& b, double tol)
{
  EXPECT_NEAR(a.real(), b.real(), tol);
  EXPECT_NEAR(a.imag(), b.imag(), tol);
}


} // namespace

TEST(Point, Operations)
{
  Point p1(1, 2, 3);
  Point p2(4, 5, 6);
  Point p3(6, 5, 4);

  test_float_eq(p1 + p2, Point(5, 7, 9));
  test_float_eq(p1 - p3, Point(-5, -3, -1));

  test_float_eq(p1 * 2, Point(2, 4, 6));
  test_float_eq(2 * p1, Point(2, 4, 6));

  test_float_eq(p1 / 2, Point(0.5, 1, 1.5));

  EXPECT_FLOAT_EQ(dot(p1, p2), 1 * 4 + 2 * 5 + 3 * 6);
  test_float_eq(cross(p1, p2), Point(-3, 6, -3));

  auto p4 = project(p1, p2);
  EXPECT_FLOAT_EQ(dot(p4, p2), std::sqrt(dot(p4, p4)) * std::sqrt(dot(p2, p2)));

  test_float_eq(-p1, Point(-p1.x, -p1.y, -p1.z));

  auto ptmp = p1;
  ptmp += p2;
  test_float_eq(ptmp, p1 + p2);

  ptmp = p1;
  ptmp -= p2;
  test_float_eq(ptmp, p1 - p2);
}

TEST(Point, Addition)
{
  Point p_real(1, 2, 3);
  PointC p_complex({4, 5}, {6, 7}, {8, 9});
  PointC p_complex2({10, 11}, {12, 13}, {14, 15});


  static_assert(std::is_same_v<decltype(p_real + p_complex), PointC>);
  test_float_eq(p_real + p_complex, PointC({5, 5}, {8, 7}, {11, 9}));

  static_assert(std::is_same_v<decltype(p_complex + p_complex2), PointC>);
  test_float_eq(p_complex + p_complex2, PointC({14, 16}, {18, 20}, {22, 24}));  
}

TEST(Point, Subtraction)
{
  Point p_real(1, 2, 3);
  PointC p_complex({4, 5}, {6, 7}, {8, 9});
  PointC p_complex2({10, 11}, {12, 13}, {14, 15});


  static_assert(std::is_same_v<decltype(p_real - p_complex), PointC>);
  test_float_eq(p_real - p_complex, PointC({-3, -5}, {-4, -7}, {-5, -9}));

  static_assert(std::is_same_v<decltype(p_complex - p_complex2), PointC>);
  test_float_eq(p_complex - p_complex2, PointC({-6, -6}, {-6, -6}, {-6, -6}));    
}


TEST(Point, MultiplyByScalar)
{
  Point p_real(1, 2, 3);
  PointC p_complex({4, 5}, {6, 7}, {8, 9});


  static_assert(std::is_same_v<decltype(p_real * 2.0), Point>);
  test_float_eq(p_real * 2.0, Point(2, 4, 6));

  static_assert(std::is_same_v<decltype(p_real * 2), Point>);
  test_float_eq(p_real * 2, Point(2, 4, 6));  


  static_assert(std::is_same_v<decltype(p_complex * 2.0), PointC>);
  test_float_eq(p_complex * 2.0, PointC({8, 10}, {12, 14}, {16, 18}));  

  static_assert(std::is_same_v<decltype(p_complex * 2), PointC>);
  test_float_eq(p_complex * 2, PointC({8, 10}, {12, 14}, {16, 18}));    
}

TEST(Point, Division)
{
  Point p_real(1, 2, 3);
  PointC p_complex({4, 5}, {6, 7}, {8, 9});

  static_assert(std::is_same_v<decltype(p_real / 2.0), Point>);
  test_float_eq(p_real / 2.0, Point(0.5, 1, 1.5));

  static_assert(std::is_same_v<decltype(p_real * 2), Point>);
  test_float_eq(p_real / 2, Point(0.5, 1, 1.5));  


  static_assert(std::is_same_v<decltype(p_complex / 2.0), PointC>);
  test_float_eq(p_complex / 2.0, PointC({2, 2.5}, {3, 3.5}, {4, 4.5}));  

  static_assert(std::is_same_v<decltype(p_complex / 2), PointC>);
  test_float_eq(p_complex / 2, PointC({2, 2.5}, {3, 3.5}, {4, 4.5}));    
}

TEST(Point, DotProduct)
{
  Point p_real(1, 2, 3);
  PointC p_complex({4, 5}, {6, 7}, {8, 9});

  static_assert(std::is_same_v<decltype(dot(p_real, p_complex)), Complex>);
  expect_near(dot(p_real, p_complex), Complex(40, 46), 1e-13);

  static_assert(std::is_same_v<decltype(dot(p_complex, p_real)), Complex>);
  expect_near(dot(p_complex, p_real), Complex(40, 46), 1e-13);  
}  

TEST(Point, Dot_rev)
{
  Point p1(1, 2, 3), p2(4, 6, 9), p1Bar, p2Bar;

  dot_rev(p1, p2, p1Bar, p2Bar, 1);

  // compute finite difference
  auto v     = dot(p1, p2);
  double eps = 1e-7;

  p1.x += eps;
  double dvDp1x = (dot(p1, p2) - v) / eps;
  p1.x -= eps;

  p1.y += eps;
  double dvDp1y = (dot(p1, p2) - v) / eps;
  p1.y -= eps;

  p1.z += eps;
  double dvDp1z = (dot(p1, p2) - v) / eps;
  p1.z -= eps;

  p2.x += eps;
  double dvDp2x = (dot(p1, p2) - v) / eps;
  p2.x -= eps;

  p2.y += eps;
  double dvDp2y = (dot(p1, p2) - v) / eps;
  p2.y -= eps;

  p2.z += eps;
  double dvDp2z = (dot(p1, p2) - v) / eps;
  p2.z -= eps;

  EXPECT_NEAR(dvDp1x, p1Bar.x, 1e-6);
  EXPECT_NEAR(dvDp1y, p1Bar.y, 1e-6);
  EXPECT_NEAR(dvDp1z, p1Bar.z, 1e-6);
  EXPECT_NEAR(dvDp2x, p2Bar.x, 1e-6);
  EXPECT_NEAR(dvDp2y, p2Bar.y, 1e-6);
  EXPECT_NEAR(dvDp2z, p2Bar.z, 1e-6);
}

TEST(Point, CrossProduct)
{
  Point p_real(1, 2, 3);
  PointC p_complex({4, 5}, {6, 7}, {8, 9});

  static_assert(std::is_same_v<decltype(cross(p_real, p_complex)), PointC>);
  test_float_eq(cross(p_real, p_complex), PointC({-2, -3}, {4, 6}, {-2, -3}));

  static_assert(std::is_same_v<decltype(cross(p_complex, p_real)), PointC>);
  test_float_eq(cross(p_complex, p_real), PointC({2, 3}, {-4, -6}, {2, 3}));
}  

}
}
}
