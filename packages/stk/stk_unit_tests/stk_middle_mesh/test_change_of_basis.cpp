#include "gtest/gtest.h"

#include "stk_middle_mesh/change_of_basis.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace utils::impl;

namespace {
void expect_float_eq(const utils::Point& pt1, const utils::Point& pt2)
{
  EXPECT_FLOAT_EQ(pt1.x, pt2.x);
  EXPECT_FLOAT_EQ(pt1.y, pt2.y);
  EXPECT_FLOAT_EQ(pt1.z, pt2.z);
}

// tests that the basis vectors are orthonormal
void test_orthonormal(std::array<utils::Point, 3> pts)
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (i == j)
        EXPECT_FLOAT_EQ(dot(pts[i], pts[j]), 1);
      else
        EXPECT_FLOAT_EQ(dot(pts[i], pts[j]), 0);

  // also test right handed
  for (int i = 0; i < 3; ++i)
  {
    int j = (i + 1) % 3;
    int k = (i + 2) % 3;
    expect_float_eq(cross(pts[i], pts[j]), pts[k]);
  }
}

void test_parallel(const utils::Point& pt1, const utils::Point& pt2)
{
  EXPECT_FLOAT_EQ(dot(pt1, pt2) / (std::sqrt(dot(pt1, pt1) * std::sqrt(dot(pt2, pt2)))), 1);
}
} // namespace

TEST(ChangeOfBasis, singleAxis)
{
  // test cases where the normal vector is along one of the coordinate axes
  utils::Point norm, p0;
  std::array<utils::Point, 3> basis;

  norm  = utils::Point(1, 0, 0) * 2;
  p0    = utils::Point(1, 1, 1);
  basis = compute_basis(norm);
  expect_float_eq(basis[0], utils::Point(0, 1, 0));
  expect_float_eq(basis[1], utils::Point(0, 0, 1));
  expect_float_eq(basis[2], utils::Point(1, 0, 0));

  norm  = utils::Point(-1, 0, 0) * 2;
  p0    = utils::Point(1, 1, 1);
  basis = compute_basis(norm);
  expect_float_eq(basis[0], utils::Point(0, 0, 1));
  expect_float_eq(basis[1], utils::Point(0, 1, 0));
  expect_float_eq(basis[2], utils::Point(-1, 0, 0));

  norm  = utils::Point(0, 1, 0) * 2;
  p0    = utils::Point(1, 1, 1);
  basis = compute_basis(norm);
  expect_float_eq(basis[0], utils::Point(0, 0, 1));
  expect_float_eq(basis[1], utils::Point(1, 0, 0));
  expect_float_eq(basis[2], utils::Point(0, 1, 0));

  norm  = utils::Point(0, -1, 0) * 2;
  p0    = utils::Point(1, 1, 1);
  basis = compute_basis(norm);
  expect_float_eq(basis[0], utils::Point(1, 0, 0));
  expect_float_eq(basis[1], utils::Point(0, 0, 1));
  expect_float_eq(basis[2], utils::Point(0, -1, 0));

  norm  = utils::Point(0, 0, 1) * 2;
  p0    = utils::Point(1, 1, 1);
  basis = compute_basis(norm);
  expect_float_eq(basis[0], utils::Point(1, 0, 0));
  expect_float_eq(basis[1], utils::Point(0, 1, 0));
  expect_float_eq(basis[2], utils::Point(0, 0, 1));

  norm  = utils::Point(0, 0, -1) * 2;
  p0    = utils::Point(1, 1, 1);
  basis = compute_basis(norm);
  expect_float_eq(basis[0], utils::Point(0, 1, 0));
  expect_float_eq(basis[1], utils::Point(1, 0, 0));
  expect_float_eq(basis[2], utils::Point(0, 0, -1));
}

TEST(ChangeOfBasis, twoAxis)
{
  // test cases where the normal vector is along two coordinate axes
  utils::Point norm;
  std::array<utils::Point, 3> basis;

  // xy
  norm  = utils::Point(1, 1, 0) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  norm  = utils::Point(1, -1, 0) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  norm  = utils::Point(-1, 1, 0) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  norm  = utils::Point(-1, -1, 0) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  // yz
  norm  = utils::Point(0, 1, 1) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  norm  = utils::Point(0, 1, -1) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  norm  = utils::Point(0, -1, 1) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  norm  = utils::Point(0, -1, -1) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  // xz
  norm  = utils::Point(1, 0, 1) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  norm  = utils::Point(1, 0, -1) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  norm  = utils::Point(-1, 0, 1) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);

  norm  = utils::Point(-1, 0, -1) * 2;
  basis = compute_basis(norm);
  test_orthonormal(basis);
  test_parallel(norm, basis[2]);
}

TEST(ChangeOfBasis, threeAxis)
{
  // general case
  std::vector<utils::Point> pts{utils::Point(1, 1, 1) * 2,   utils::Point(-1, 1, 1) * 2,  utils::Point(1, -1, 1) * 2,
                                utils::Point(-1, -1, 1) * 2, utils::Point(1, 1, -1) * 2,  utils::Point(1, 1, -1) * 2,
                                utils::Point(1, -1, 1) * 2,  utils::Point(-1, -1, -1) * 2};

  for (auto& norm : pts)
  {
    auto basis = compute_basis(norm);
    test_orthonormal(basis);
    test_parallel(norm, basis[2]);
  }
}

TEST(ChangeOfBasis, projection)
{
  // general case
  std::vector<utils::Point> normals{
      utils::Point(1, 1, 1) * 2,  utils::Point(-1, 1, 1) * 2, utils::Point(1, -1, 1) * 2, utils::Point(-1, -1, 1) * 2,
      utils::Point(1, 1, -1) * 2, utils::Point(1, 1, -1) * 2, utils::Point(1, -1, 1) * 2, utils::Point(-1, -1, -1) * 2};

  utils::Point pt(2, 3, 4);
  for (auto& norm : normals)
  {
    auto basis = compute_basis(norm);
    utils::impl::ChangeOfBasis proj(basis);
    auto ptPrime = proj.project_forward(pt);
    auto pt2     = proj.project_back(ptPrime);

    expect_float_eq(pt, pt2);
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
