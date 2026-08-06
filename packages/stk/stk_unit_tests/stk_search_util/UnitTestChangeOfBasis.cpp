#include "gtest/gtest.h"
#include "stk_search_util/ChangeOfBasis.hpp"

#include <cmath>

namespace {

double dot(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
  double val = 0;
  for (unsigned i=0; i < 3; ++i)
  {
    val += a[i]*b[i];
  }

  return val;
}

void check_orthonormal(const stk::search::impl::ChangeOfBasis& basis)
{
  std::array<std::array<double, 3>, 3> vecs;
  vecs[0] = basis.change_basis_forward({1, 0, 0});
  vecs[1] = basis.change_basis_forward({0, 1, 0});
  vecs[2] = basis.change_basis_forward({0, 0, 1});

  for (unsigned i=0; i < 3; ++i)
  {
    for (unsigned j=0; j < 3; ++j)
    {
      if (i == j)
      {
        double mag = std::sqrt(dot(vecs[i], vecs[i]));
        EXPECT_NEAR(mag, 1, 1e-13);
      } else
      {
        EXPECT_NEAR(dot(vecs[i], vecs[j]), 0, 1e-13);
      }
    }
  }
}

}

TEST(MatMat, Exactness)
{
  std::array<double, 9> a = {1, 2, 3,
                             4, 5, 6,
                             7, 8, 9};
  std::array<double, 9> b = {2, 3, 4,
                             5, 6, 7,
                             8, 9, 10};

  std::array<double, 9> c = stk::search::impl::matmat(a, b);
  std::array<double, 9> c_expected = {36, 42, 48, 81, 96, 111, 126, 150, 174};
  for (unsigned i=0; i < 9; ++i)
  {
    EXPECT_DOUBLE_EQ(c[i], c_expected[i]);
  }
}

TEST(MatVec, Exactness)
{
  std::array<double, 9> a = {1, 2, 3,
                             4, 5, 6,
                             7, 8, 9};
  std::array<double, 3> b = {2, 5, 8};

  std::array<double, 3> c = stk::search::impl::matvec(a, b);
  std::array<double, 3> c_expected = {36, 81, 126};
  for (unsigned i=0; i < 3; ++i)
  {
    EXPECT_DOUBLE_EQ(c[i], c_expected[i]);
  }
}

TEST(MatVecTransposed, Exactness)
{
  std::array<double, 9> a = {1, 2, 3,
                             4, 5, 6,
                             7, 8, 9};
  std::array<double, 3> b = {2, 5, 8};

  std::array<double, 3> c = stk::search::impl::matvec_transposed(a, b);
  std::array<double, 3> c_expected = {78, 93, 108};
  for (unsigned i=0; i < 3; ++i)
  {
    EXPECT_DOUBLE_EQ(c[i], c_expected[i]);
  }
}

TEST(ComputeRotationMatrix, ZToZ)
{
  std::array<double, 9> rot = stk::search::impl::compute_rotation_matrix({0, 0, 1});
  for (unsigned i=0; i < 3; ++i)
  {
    for (unsigned j=0; j < 3; ++j)
    {
      if (i == j)
      {
        EXPECT_DOUBLE_EQ(rot[3*i + j], 1);
      } else
      {
        EXPECT_DOUBLE_EQ(rot[3*i + j], 0);
      }
    }
  }
}

TEST(ComputeRotationMatrix, ZToX)
{
  std::array<double, 9> rot = stk::search::impl::compute_rotation_matrix({1, 0, 0});
  std::array<double, 3> pt_prime = stk::search::impl::matvec(rot, {0, 0, 1});
  EXPECT_NEAR(pt_prime[0], 1, 1e-13);
  EXPECT_NEAR(pt_prime[1], 0, 1e-13);
  EXPECT_NEAR(pt_prime[2], 0, 1e-13);
}

TEST(ComputeRotationMatrix, ZToY)
{
  std::array<double, 9> rot = stk::search::impl::compute_rotation_matrix({0, 1, 0});
  std::array<double, 3> pt_prime = stk::search::impl::matvec(rot, {0, 0, 1});
  EXPECT_NEAR(pt_prime[0], 0, 1e-13);
  EXPECT_NEAR(pt_prime[1], 1, 1e-13);
  EXPECT_NEAR(pt_prime[2], 0, 1e-13);
}

TEST(ChangeOfBasis, CombinedRotation)
{
  std::array<double, 3> new_z = {1, 1, 0};
  stk::search::impl::ChangeOfBasis basis(stk::search::impl::compute_rotation_matrix(new_z));
  check_orthonormal(basis);

  std::array<double, 3> new_z_projected = basis.change_basis_forward(new_z);
  EXPECT_NEAR(new_z_projected[0], 0, 1e-13);
  EXPECT_NEAR(new_z_projected[1], 0, 1e-13);
  EXPECT_NEAR(new_z_projected[2], std::sqrt(2), 1e-13);

  std::array<double, 3> new_z2 = basis.change_basis_reverse(new_z_projected);
  EXPECT_NEAR(new_z2[0], new_z[0], 1e-13);
  EXPECT_NEAR(new_z2[1], new_z[1], 1e-13);
  EXPECT_NEAR(new_z2[2], new_z[2], 1e-13);
}

