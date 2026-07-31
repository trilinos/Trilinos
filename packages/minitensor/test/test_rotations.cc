// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests for the rotations module: SO(3) logarithmic and exponential
// maps, quaternions, and conversions among rotation representations.

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

TEST(MiniTensor, LogRotation)
{
  // Identity rotation
  Tensor<Real>
  I = identity<Real>(3);

  Tensor<Real>
  i = log_rotation(I);

  Real const
  error_I = norm(i) / norm(I);

  ASSERT_LE(error_I, machine_epsilon<Real>());

  // Pi / 4 rotation about Z.
  Real const
  c = std::sqrt(2.0) / 2.0;

  Tensor<Real> const R(c, -c, 0.0, c, c, 0.0, 0.0, 0.0, 1.0);

  Tensor<Real> const r = log_rotation(R);

  Real const
  Pi = std::acos(-1.0);

  Real const
  error_R = std::abs(r(0,1) + Pi / 4.0);

  ASSERT_LE(error_R, machine_epsilon<Real>());

  Real const
  error_r = std::abs(r(0,1) + r(1,0));

  ASSERT_LE(error_r, machine_epsilon<Real>());

  // Pi rotation about Z
  I(0, 0) = -1.0;
  I(1, 1) = -1.0;

  i = log_rotation(I);

  Real const
  error_pi = 0.5 * (std::abs(i(0, 1) + Pi) + std::abs(i(1, 0) - Pi));

  ASSERT_LE(error_pi, machine_epsilon<Real>());
}

TEST(MiniTensor, QuaternionAlgebra)
{
  // Default construction is the identity quaternion.
  Quaternion<Real> const
  id;

  ASSERT_EQ(id.scalar(), 1.0);
  ASSERT_EQ(norm(id.vector()), 0.0);

  Quaternion<Real> const
  q(0.5, Vector<Real, 3>(0.5, -0.5, 0.5));

  // Identity is the neutral element of the Hamilton product.
  ASSERT_TRUE(q * id == q);
  ASSERT_TRUE(id * q == q);
  ASSERT_TRUE(q != id);

  // Indexed access is scalar-first and consistent with to_vector
  // and the 4-vector constructor.
  Vector<Real, 4> const
  c = q.to_vector();

  for (Index i = 0; i < 4; ++i) {
    ASSERT_EQ(c(i), q(i));
  }

  Quaternion<Real> const
  p(c);

  ASSERT_TRUE(p == q);

  // q conjugate(q) has zero vector part and scalar part |q|^2.
  Quaternion<Real> const
  qqbar = q * conjugate(q);

  Real const
  norm_q = norm(q);

  ASSERT_LE(
      std::abs(qqbar.scalar() - norm_q * norm_q),
      machine_epsilon<Real>());

  ASSERT_LE(norm(qqbar.vector()), machine_epsilon<Real>());

  // q inverse(q) is the identity.
  Quaternion<Real> const
  qqinv = q * inverse(q);

  ASSERT_LE(std::abs(qqinv.scalar() - 1.0), 2.0 * machine_epsilon<Real>());
  ASSERT_LE(norm(qqinv.vector()), 2.0 * machine_epsilon<Real>());

  // unit(q) has unit norm.
  ASSERT_LE(
      std::abs(norm(unit(q)) - 1.0),
      2.0 * machine_epsilon<Real>());

  // The Hamilton product corresponds to rotation composition.
  Vector<Real, 3> const
  rv1(0.4, -0.3, 0.8);

  Vector<Real, 3> const
  rv2(-0.2, 0.7, 0.1);

  Quaternion<Real> const
  q1 = q_of_rv(rv1);

  Quaternion<Real> const
  q2 = q_of_rv(rv2);

  Tensor<Real, 3> const
  R12 = rt_of_q(q1 * q2);

  Tensor<Real, 3> const
  R1R2 = rt_of_q(q1) * rt_of_q(q2);

  Real const
  error = norm(R12 - R1R2) / norm(R1R2);

  ASSERT_LE(error, 16 * machine_epsilon<Real>());
}

TEST(MiniTensor, QuaternionRoundTrips)
{
  Real const
  Pi = std::acos(-1.0);

  Vector<Real, 3> const
  axis = unit(Vector<Real, 3>(1.0, 1.0, 1.0));

  // Rotation angles exercising: zero, both asymptotic tiers of
  // sin(x)/x, generic angles, near pi and exactly pi.
  Real const
  angles[] = {0.0, 1.0e-8, 1.0e-4, 0.01, 0.5, 2.0, 3.0,
      Pi - 1.0e-6, Pi};

  for (Real const angle : angles) {

    Vector<Real, 3> const
    w = angle * axis;

    Tensor<Real, 3> const
    R = rt_of_rv(w);

    // Rotation matrix invariants.
    Real const
    error_orthogonal = norm(dot_t(R, R) - identity<Real, 3>(3));

    ASSERT_LE(error_orthogonal, 16 * machine_epsilon<Real>());

    ASSERT_LE(std::abs(det(R) - 1.0), 16 * machine_epsilon<Real>());

    // Unit quaternion with non-negative scalar part after
    // sign normalization.
    Quaternion<Real> const
    q = q_of_rt(R);

    ASSERT_LE(std::abs(norm(q) - 1.0), 16 * machine_epsilon<Real>());

    // Matrix round trip, valid for all angles including pi.
    Vector<Real, 3> const
    v = rv_of_rt(R);

    Tensor<Real, 3> const
    S = rt_of_rv(v);

    ASSERT_LE(norm(S - R), 64 * machine_epsilon<Real>());

    // Principal property.
    ASSERT_LE(norm(v), Pi + 16 * machine_epsilon<Real>());

    // Vector round trip, unique for angles below pi.
    if (angle < Pi) {
      ASSERT_LE(norm(v - w), 64 * machine_epsilon<Real>() * (1.0 + angle));
    } else {
      ASSERT_LE(std::abs(norm(v) - Pi), 64 * machine_epsilon<Real>());
    }
  }

  // Exactly pi about each coordinate axis: the rotation matrix
  // round trip must be exact even though the rotation vector is
  // determined only up to sign.
  for (Index i = 0; i < 3; ++i) {

    Vector<Real, 3>
    w(0.0, 0.0, 0.0);

    w(i) = Pi;

    Tensor<Real, 3> const
    R = rt_of_rv(w);

    Tensor<Real, 3> const
    S = rt_of_rv(rv_of_rt(R));

    ASSERT_LE(norm(S - R), 64 * machine_epsilon<Real>());
  }
}

TEST(MiniTensor, QuaternionSpurrierBranches)
{
  Real const
  Pi = std::acos(-1.0);

  // Small rotation: the trace is dominant.
  // Rotations by pi about each coordinate axis: the corresponding
  // diagonal entry is dominant.
  Vector<Real, 3> const
  rvs[] = {
      Vector<Real, 3>(0.1, 0.2, -0.1),
      Vector<Real, 3>(Pi, 0.0, 0.0),
      Vector<Real, 3>(0.0, Pi, 0.0),
      Vector<Real, 3>(0.0, 0.0, Pi)};

  for (Vector<Real, 3> const & rv : rvs) {

    Tensor<Real, 3> const
    R = rt_of_rv(rv);

    Quaternion<Real> const
    q = q_of_rt(R);

    ASSERT_LE(std::abs(norm(q) - 1.0), 16 * machine_epsilon<Real>());

    Real const
    error = norm(rt_of_q(q) - R);

    ASSERT_LE(error, 16 * machine_epsilon<Real>());
  }
}

TEST(MiniTensor, QuaternionRotationConsistency)
{
  Vector<Real, 3> const
  w(0.4, -0.3, 0.8);

  // The quaternion path agrees with the exponential map ...
  Tensor<Real, 3> const
  R = rt_of_rv(w);

  Tensor<Real, 3> const
  E = exp_skew_symmetric(skew(w));

  ASSERT_LE(norm(R - E) / norm(E), 16 * machine_epsilon<Real>());

  // ... and with the logarithmic map away from pi.
  Tensor<Real, 3> const
  r = log_rotation(R);

  ASSERT_LE(norm(skew(rv_of_rt(R)) - r), 16 * machine_epsilon<Real>());

  // vee inverts skew from vector.
  ASSERT_TRUE(vee(skew(w)) == w);
}

TEST(MiniTensor, RotationVectorContinuation)
{
  Real const
  Pi = std::acos(-1.0);

  Vector<Real, 3> const
  zero(0.0, 0.0, 0.0);

  Vector<Real, 3> const
  old(0.1, 0.2, 0.3);

  Vector<Real, 3> const
  direction = unit(old);

  // Close previous vector: no continuation.
  ASSERT_TRUE(rv_continue(old, old) == old);

  // Zero previous vector: no continuation.
  ASSERT_TRUE(rv_continue(old, zero) == old);

  // Zero rotation with a previous full turn: the full turn
  // is recovered.
  Vector<Real, 3> const
  turn = 2.0 * Pi * Vector<Real, 3>(0.0, 1.0, 0.0);

  Real const
  error_zero = norm(rv_continue(zero, turn) - turn);

  ASSERT_LE(error_zero, 16 * machine_epsilon<Real>());

  // Previous vector one and two full turns ahead: the equivalent
  // rotation vector closest to it is recovered exactly.
  for (Real const turns : {1.0, 2.0}) {

    Vector<Real, 3> const
    prev = old + turns * 2.0 * Pi * direction;

    Vector<Real, 3> const
    continued = rv_continue(old, prev);

    Real const
    error = norm(continued - prev) / norm(prev);

    ASSERT_LE(error, 16 * machine_epsilon<Real>());

    // The continued vector represents the same rotation.
    Real const
    error_rotation = norm(rt_of_rv(continued) - rt_of_rv(old));

    ASSERT_LE(error_rotation, 256 * machine_epsilon<Real>());
  }

  // Previous vector a full turn behind, anti-parallel.
  Vector<Real, 3> const
  prev = old - 2.0 * Pi * direction;

  Vector<Real, 3> const
  continued = rv_continue(old, prev);

  Real const
  error = norm(continued - prev) / norm(prev);

  ASSERT_LE(error, 16 * machine_epsilon<Real>());

  Real const
  error_rotation = norm(rt_of_rv(continued) - rt_of_rv(old));

  ASSERT_LE(error_rotation, 256 * machine_epsilon<Real>());

  // Continuity across the pi and 2 pi crossings: sweep the rotation
  // angle, continuing each principal rotation vector from the
  // previous continued one, and recover the swept vector throughout.
  Vector<Real, 3> const
  sweep_axis = unit(Vector<Real, 3>(1.0, -2.0, 2.0));

  Vector<Real, 3>
  swept(0.0, 0.0, 0.0);

  for (Real angle = 0.05; angle < 7.0; angle += 0.1) {

    Vector<Real, 3> const
    principal = rv_of_rt(rt_of_rv(angle * sweep_axis));

    swept = rv_continue(principal, swept);

    Real const
    error_sweep = norm(swept - angle * sweep_axis);

    ASSERT_LE(error_sweep, 1.0e-12);
  }
}

TEST(MiniTensor, QuaternionFad)
{
  using FAD = Sacado::Fad::DFad<Real>;

  Real const
  w0[3] = {0.1, -0.2, 0.3};

  Vector<FAD, 3>
  w(3);

  for (Index i = 0; i < 3; ++i) {
    w(i) = FAD(3, i, w0[i]);
  }

  // The whole conversion cycle instantiates with a FAD scalar.
  Tensor<FAD, 3> const
  R = rt_of_q(q_of_rv(w));

  Vector<FAD, 3> const
  v = rv_of_q(q_of_rt(R));

  FAD const
  error = norm(v - w) / norm(w);

  ASSERT_LE(error.val(), 64 * machine_epsilon<Real>());

  // Derivative spot check of rt_of_rv against central differences.
  Real const
  h = 1.0e-6;

  for (Index k = 0; k < 3; ++k) {

    Vector<Real, 3>
    wp(w0[0], w0[1], w0[2]);

    Vector<Real, 3>
    wm(w0[0], w0[1], w0[2]);

    wp(k) += h;
    wm(k) -= h;

    Tensor<Real, 3> const
    Rp = rt_of_rv(wp);

    Tensor<Real, 3> const
    Rm = rt_of_rv(wm);

    for (Index i = 0; i < 3; ++i) {
      for (Index j = 0; j < 3; ++j) {

        Real const
        finite_difference = (Rp(i, j) - Rm(i, j)) / (2.0 * h);

        ASSERT_NEAR(R(i, j).dx(k), finite_difference, 1.0e-9);
      }
    }
  }
}

} // namespace minitensor
