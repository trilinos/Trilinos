// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests for the geometry and mechanics modules: element measures
// (lengths, areas, volumes), volumetric-deviatoric decomposition and
// mechanics transforms.

#include "gtest/gtest.h"
#include "MiniTensor.h"

int
main(int ac, char* av[])
{
  Kokkos::initialize();

  ::testing::GTEST_FLAG(print_time) = (ac > 1) ? true : false;

  ::testing::InitGoogleTest(&ac, av);

  auto const
  retval = RUN_ALL_TESTS();

  Kokkos::finalize();

  return retval;
}

namespace minitensor {

TEST(MiniTensor, VolumetricDeviatoric)
{
  Tensor<Real> A = 3.0 * eye<Real>(3);

  ASSERT_LE(norm(A - vol(A)), machine_epsilon<Real>());

  Tensor<Real> const B = dev(A);

  A(0, 0) = 0.0;
  A(1, 1) = 0.0;
  A(2, 2) = 0.0;

  ASSERT_LE(norm(A - B), machine_epsilon<Real>());
}

TEST(MiniTensor, MechanicsTransforms)
{
  Tensor<Real> const F(0.0, -6.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0 / 3.0);

  Tensor<Real> sigma(0.0, 0.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0);

  Tensor<Real> const P = piola(F, sigma);

  Real error = std::abs(P(1, 0) - 100.0) / 100.0;

  ASSERT_LE(error, machine_epsilon<Real>());

  sigma = piola_inverse(F, P);

  error = std::abs(sigma(1, 1) - 50.0) / 50.0;

  ASSERT_LE(error, machine_epsilon<Real>());

  Tensor<Real> const E = 0.5 * (t_dot(F, F) - eye<Real>(3));

  Tensor<Real> const e = 0.5 * (eye<Real>(3) - inverse(dot_t(F, F)));

  Tensor<Real> const g = push_forward_covariant(F, E);

  error = norm(g - e) / norm(e);

  ASSERT_LE(error, machine_epsilon<Real>());

  Tensor<Real> const G = pull_back_covariant(F, e);

  error = norm(G - E) / norm(E);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, SegmentLength)
{
  Vector<Real, 3> const u(0.0, 0.0, 0.0);
  Vector<Real, 3> const v(1.0, 2.0, 3.0);

  Real const l = length(u, v);
  Real const n = norm(v);
  Real const error = std::abs(l - n);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, TriangleArea)
{
  Vector<Real, 3> const a(0.0, 0.0, 0.0);
  Vector<Real, 3> const b(1.0, 0.0, 0.0);
  Vector<Real, 3> const c(0.0, 1.0, 0.0);

  Real const A = area(a, b, c);
  Real const error = std::abs(A - 0.5);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, QuadArea)
{
  Vector<Real, 3> const a(1.0, 0.0, 0.0);
  Vector<Real, 3> const b(2.0, 1.0, 0.0);
  Vector<Real, 3> const c(1.0, 2.0, 0.0);
  Vector<Real, 3> const d(0.0, 1.0, 0.0);

  Real const A = area(a, b, c, d);
  Real const error = std::abs(A - 2.0);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, TetVolume)
{
  Vector<Real, 3> const a(0.0, 0.0, 0.0);
  Vector<Real, 3> const b(1.0, 0.0, 0.0);
  Vector<Real, 3> const c(0.0, 1.0, 0.0);
  Vector<Real, 3> const d(0.0, 0.0, 1.0);

  Real const V = volume(a, b, c, d);
  Real const error = std::abs(V - 1.0 / 6.0);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, HexVolume)
{
  Vector<Real, 3> const a(0.0, 0.0, 0.0);
  Vector<Real, 3> const b(1.0, 0.0, 0.0);
  Vector<Real, 3> const c(1.0, 1.0, 0.0);
  Vector<Real, 3> const d(0.0, 1.0, 0.0);
  Vector<Real, 3> const e(0.0, 0.0, 1.0);
  Vector<Real, 3> const f(1.0, 0.0, 1.0);
  Vector<Real, 3> const g(1.0, 1.0, 1.0);
  Vector<Real, 3> const h(0.0, 1.0, 1.0);

  Real const V = volume(a, b, c, d, e, f, g, h);
  Real const error = std::abs(V - 1.0);

  ASSERT_LE(error, machine_epsilon<Real>());
}

} // namespace minitensor
