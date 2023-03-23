#include "gtest/gtest.h"

#include "projection.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

void test_float_eq(const utils::Point& p1, const utils::Point& p2)
{
  EXPECT_FLOAT_EQ(p1.x, p2.x);
  EXPECT_FLOAT_EQ(p1.y, p2.y);
  EXPECT_FLOAT_EQ(p1.z, p2.z);
}

} // namespace

TEST(Point, Operations)
{
  utils::Point p1(1, 2, 3);
  utils::Point p2(4, 5, 6);
  utils::Point p3(6, 5, 4);

  test_float_eq(p1 + p2, utils::Point(5, 7, 9));
  test_float_eq(p1 - p3, utils::Point(-5, -3, -1));

  test_float_eq(p1 * 2, utils::Point(2, 4, 6));
  test_float_eq(2 * p1, utils::Point(2, 4, 6));

  test_float_eq(p1 / 2, utils::Point(0.5, 1, 1.5));

  EXPECT_FLOAT_EQ(dot(p1, p2), 1 * 4 + 2 * 5 + 3 * 6);
  test_float_eq(cross(p1, p2), utils::Point(-3, 6, -3));

  auto p4 = project(p1, p2);
  EXPECT_FLOAT_EQ(dot(p4, p2), std::sqrt(dot(p4, p4)) * std::sqrt(dot(p2, p2)));

  test_float_eq(-p1, utils::Point(-p1.x, -p1.y, -p1.z));

  auto ptmp = p1;
  ptmp += p2;
  test_float_eq(ptmp, p1 + p2);

  ptmp = p1;
  ptmp -= p2;
  test_float_eq(ptmp, p1 - p2);
}

TEST(Projection, ProjectToPlane)
{
  // test projection onto xy planes
  utils::Point p1(0, 0, 0);
  utils::Point p2(2, 0, 0);
  utils::Point p3(0, 3, 0);

  utils::Point p4(1, 2, 3);
  utils::impl::Projection proj(p1, p2, p3);

  test_float_eq(proj.project_to_plane(p4), utils::Point(1, 2, 0));

  // test projection onto yz plane
  p2   = utils::Point(0, 2, 0);
  p3   = utils::Point(0, 0, 3);
  proj = utils::impl::Projection(p1, p2, p3);
  test_float_eq(proj.project_to_plane(p4), utils::Point(0, 2, 3));

  // test projection onto xz plane
  p2   = utils::Point(2, 0, 0);
  p3   = utils::Point(0, 0, 3);
  proj = utils::impl::Projection(p1, p2, p3);
  test_float_eq(proj.project_to_plane(p4), utils::Point(1, 0, 3));
}

TEST(Projection, ComputePlaneCoords)
{
  // test xy plane
  utils::Point p1(0, 0, 0);
  utils::Point p2(2, 0, 0);
  utils::Point p3(0, 3, 0);

  utils::Point p4(1, 2, 3);
  utils::impl::Projection proj(p1, p2, p3);

  test_float_eq(proj.project_plane_coords(p4), utils::Point(0.5, 2.0 / 3.0, 0));

  // test non-orthogonal input basis
  p2 = utils::Point(1, 1, 0);
  p3 = utils::Point(0, 2, 0);

  // p4 = utils::Point(-0.25, 1.25, 0)/std::sqrt(2);
  p4   = utils::Point(-1, 5, 0);
  proj = utils::impl::Projection(p1, p2, p3);

  // test_float_eq(proj.project_plane_coords(p4), utils::Point(0.5, 0.75, 0));
  test_float_eq(proj.project_plane_coords(p4), utils::Point(2, 3, 0));
}

TEST(Projection, ProjectToPlane_rev)
{
  // test projection onto xy planes
  utils::Point p1(1, 2, 3);
  utils::Point p2(2, 3, 4);
  utils::Point p3(5, 5, 6);

  utils::Point p4(7, 8, 9);
  utils::impl::Projection proj(p1, p2, p3);
  auto p5 = proj.project_to_plane(p4);

  utils::Point pBar(1, 2, 3);
  auto p4Bar = proj.project_to_plane_rev(p4, pBar);

  // compute finite difference
  double eps = 1e-7;
  p4.x += eps;
  auto dp5Dp4x = (proj.project_to_plane(p4) - p5) / eps;
  p4.x -= eps;

  p4.y += eps;
  auto dp5Dp4y = (proj.project_to_plane(p4) - p5) / eps;
  p4.y -= eps;

  p4.z += eps;
  auto dp5Dp4z = (proj.project_to_plane(p4) - p5) / eps;
  p4.z -= eps;

  utils::Point p4BarFd((dp5Dp4x.x * pBar.x + dp5Dp4x.y * pBar.y + dp5Dp4x.z * pBar.z),
                       (dp5Dp4y.x * pBar.x + dp5Dp4y.y * pBar.y + dp5Dp4y.z * pBar.z),
                       (dp5Dp4z.x * pBar.x + dp5Dp4z.y * pBar.y + dp5Dp4z.z * pBar.z));

  EXPECT_NEAR(p4Bar.x, p4BarFd.x, 1e-6);
  EXPECT_NEAR(p4Bar.y, p4BarFd.y, 1e-6);
  EXPECT_NEAR(p4Bar.z, p4BarFd.z, 1e-6);
}

TEST(Projection, ComputePlaneCoords_rev)
{
  // test projection onto xy planes
  utils::Point p1(1, 2, 3);
  utils::Point p2(2, 3, 4);
  utils::Point p3(5, 5, 6);

  utils::Point p4(7, 8, 9);
  utils::impl::Projection proj(p1, p2, p3);
  p4 = proj.project_to_plane(p4);

  utils::Point pBar(1, 2, 3);
  auto p4Bar = proj.compute_plane_coords_rev(p4, pBar);

  // compute finite difference
  double eps = 1e-7;
  auto p5    = proj.compute_plane_coords(p4);
  p4.x += eps;
  auto dp5Dp4x = (proj.compute_plane_coords(p4) - p5) / eps;
  p4.x -= eps;

  p4.y += eps;
  auto dp5Dp4y = (proj.compute_plane_coords(p4) - p5) / eps;
  p4.y -= eps;

  p4.z += eps;
  auto dp5Dp4z = (proj.compute_plane_coords(p4) - p5) / eps;
  p4.z -= eps;

  utils::Point p4BarFd((dp5Dp4x.x * pBar.x + dp5Dp4x.y * pBar.y + dp5Dp4x.z * pBar.z),
                       (dp5Dp4y.x * pBar.x + dp5Dp4y.y * pBar.y + dp5Dp4y.z * pBar.z),
                       (dp5Dp4z.x * pBar.x + dp5Dp4z.y * pBar.y + dp5Dp4z.z * pBar.z));

  EXPECT_NEAR(p4Bar.x, p4BarFd.x, 1e-6);
  EXPECT_NEAR(p4Bar.y, p4BarFd.y, 1e-6);
  EXPECT_NEAR(p4Bar.z, p4BarFd.z, 1e-6);
}

TEST(Projection, project_plane_coords_rev)
{
  utils::Point p1(1, 2, 3);
  utils::Point p2(2, 3, 4);
  utils::Point p3(5, 5, 6);

  utils::Point p4(7, 8, 9);
  utils::impl::Projection proj(p1, p2, p3);
  // proj.project_to_plane(p4);

  utils::Point pBar(1, 2, 3);
  auto p4Bar = proj.project_plane_coords_rev(p4, pBar);

  // compute finite difference
  double eps = 1e-7;
  auto p5    = proj.project_plane_coords(p4);
  p4.x += eps;
  auto dp5Dp4x = (proj.project_plane_coords(p4) - p5) / eps;
  p4.x -= eps;

  p4.y += eps;
  auto dp5Dp4y = (proj.project_plane_coords(p4) - p5) / eps;
  p4.y -= eps;

  p4.z += eps;
  auto dp5Dp4z = (proj.project_plane_coords(p4) - p5) / eps;
  p4.z -= eps;

  utils::Point p4BarFd((dp5Dp4x.x * pBar.x + dp5Dp4x.y * pBar.y + dp5Dp4x.z * pBar.z),
                       (dp5Dp4y.x * pBar.x + dp5Dp4y.y * pBar.y + dp5Dp4y.z * pBar.z),
                       (dp5Dp4z.x * pBar.x + dp5Dp4z.y * pBar.y + dp5Dp4z.z * pBar.z));

  EXPECT_NEAR(p4Bar.x, p4BarFd.x, 1e-6);
  EXPECT_NEAR(p4Bar.y, p4BarFd.y, 1e-6);
  EXPECT_NEAR(p4Bar.z, p4BarFd.z, 1e-6);
}

TEST(Projection, Dot_rev)
{
  utils::Point p1(1, 2, 3), p2(4, 6, 9), p1Bar, p2Bar;

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

} // namespace impl
} // namespace middle_mesh
} // namespace stk
