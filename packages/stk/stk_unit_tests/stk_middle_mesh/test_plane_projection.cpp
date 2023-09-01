#include "gtest/gtest.h"

#include "stk_middle_mesh/plane_projection.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

void expect_float_eq2(const utils::Point& pt, const utils::Point& pt2)
{
  EXPECT_FLOAT_EQ(pt.x, pt2.x);
  EXPECT_FLOAT_EQ(pt.y, pt2.y);
}

void expect_float_eq3(const utils::Point& pt, const utils::Point& pt2)
{
  EXPECT_FLOAT_EQ(pt.x, pt2.x);
  EXPECT_FLOAT_EQ(pt.y, pt2.y);
  EXPECT_FLOAT_EQ(pt.z, pt2.z);
}

void test_projection(utils::impl::PlaneProjection val, const utils::Point& pt)
{
  auto pts                      = utils::impl::get_plane_projection_points(val);
  utils::impl::Projection* proj = utils::impl::get_plane_projection(pts);

  auto pt2 = proj->compute_plane_coords(pt);
  auto pt3 = utils::impl::apply_plane_projection(val, pt);

  expect_float_eq2(pt2, pt3);

  // test inverse function
  auto pt4 = utils::impl::apply_inverse_plane_projection(val, pt3);
  expect_float_eq3(pt, pt4);

  delete proj;
}

} // namespace

TEST(PlaneProjection, apply_plane_projection)
{
  utils::Point pt(1, 2, 3);

  // test equivalence between using Projection and applyutils::impl::PlaneProjection
  test_projection(utils::impl::PlaneProjection::YzPos, pt);
  test_projection(utils::impl::PlaneProjection::YzNeg, pt);
  test_projection(utils::impl::PlaneProjection::XzPos, pt);
  test_projection(utils::impl::PlaneProjection::XzNeg, pt);
  test_projection(utils::impl::PlaneProjection::XyPos, pt);
  test_projection(utils::impl::PlaneProjection::XyNeg, pt);
}

TEST(PlaneProjection, getPlaneProjectionEnum)
{
  utils::Point pt1(-0.7513110251673053, -0.8136577495531562, 0);
  utils::Point pt2(-0.2847545178640647, -1.081679676310688, 0);
  utils::Point pt3(-0.2688256750539578, -1.003439971294373, 0);

  EXPECT_EQ(utils::impl::get_plane_projection_enum(pt1, pt2, pt3), utils::impl::PlaneProjection::XyPos);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
