// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests for the matrix-functions module: tensor exponential,
// logarithm and their variants, and the Baker-Campbell-Hausdorff
// series.

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

TEST(MiniTensor, Exponential)
{
  Tensor<Real> const A = eye<Real>(3) + Tensor<Real>(3, Filler::ONES);

  Tensor<Real> const B = exp_pade(A);

  Tensor<Real> const C = exp_taylor(A);

  Tensor<Real> const D = B - C;

  Real const error = norm(D) / norm(B);

  ASSERT_LE(error, 2.0 * machine_epsilon<Real>());
}

TEST(MiniTensor, Log)
{
  Tensor<Real>
  I = identity<Real>(3);

  Tensor<Real>
  r = sqrt(I);

  Real const
  error_sqrt = norm(r - I);

  ASSERT_LE(error_sqrt, machine_epsilon<Real>());

  Tensor<Real>
  i = log(I);

  Real const
  error_I = norm(i) / norm(I);

  ASSERT_LE(error_I, machine_epsilon<Real>());

  Tensor<Real>
  F(-0.16777540263807703,   1.2889030921484332,    0.09298444646599896,
    -0.718646161000825955,  0.02120989960140519,  -0.039217352050714333,
    -0.073802850037046119,  0.036685806156092855,  1.0456450021778172);

  Tensor<Real>
  f = exp(log(F));

  Real const
  error_F = norm(f - F);

  ASSERT_LE(error_F, 8 * machine_epsilon<Real>());
}

TEST(MiniTensor, ExpSym)
{
  Tensor<Real> const
  O(3, Filler::ZEROS);

  // exp of the zero tensor is the identity
  Tensor<Real> const
  I = exp_sym(O);

  Real const
  error_I = norm(I - identity<Real>(3));

  ASSERT_LE(error_I, machine_epsilon<Real>());

  // symmetric tensor with eigenvalues of both signs
  Tensor<Real> const
  A(0.7, -0.3,  0.2,
   -0.3,  0.1, -0.5,
    0.2, -0.5, -0.9);

  // agree with the general scaling-and-squaring exponential
  Tensor<Real> const
  E = exp_sym(A);

  Real const
  error_exp = norm(E - exp(A)) / norm(E);

  ASSERT_LE(error_exp, 16 * machine_epsilon<Real>());

  // round trip through the symmetric logarithm recovers the argument
  Tensor<Real> const
  R = log_sym(E);

  Real const
  error_rt = norm(R - A) / norm(A);

  ASSERT_LE(error_rt, 16 * machine_epsilon<Real>());
}

TEST(MiniTensor, LogIllConditioned)
{
  // SPD tensor with eigenvalues (1e6, 1, 1e-6) and unit determinant. The
  // expanded 3x3 determinant of this tensor cancels catastrophically to
  // exactly zero in double precision, which used to make the determinant
  // scaling in sqrt_dbp non-finite and abort (or corrupt memory in release
  // builds) inside log().
  Tensor<Real>
  A(3.05423526260858198e+04, 1.44662229379329743e+05, 9.31790062083673984e+04,
    1.44662229379329743e+05, 6.85186842144699185e+05, 4.41337014755186858e+05,
    9.31790062083673984e+04, 4.41337014755186858e+05, 2.84271805230215017e+05);

  Tensor<Real> const
  L = log(A);

  for (Index i = 0; i < 3; ++i) {
    for (Index j = 0; j < 3; ++j) {
      ASSERT_TRUE(std::isfinite(L(i, j)));
    }
  }

  // The eigendecomposition path is independent of the scaling-and-squaring
  // path; at this conditioning both are accurate to ~1e-5 relative.
  Tensor<Real> const
  L_sym = log_sym(A);

  Real const
  error_log = norm(L - L_sym) / norm(L_sym);

  ASSERT_LE(error_log, 1.0e-4);
}

TEST(MiniTensor, LogGeneralComplexPair)
{
  // Nonsymmetric logarithm with a complex-conjugate eigenvalue pair.
  // The eigenvalue imaginary parts of L are bounded by its norm, which is
  // well inside (-pi, pi), so log(exp(L)) must recover L exactly.
  Tensor<Real> const
  L(0.20, -1.10,  0.30,
    1.10,  0.20, -0.40,
    0.10,  0.25, -0.30);

  Tensor<Real> const
  A = exp(L);

  // the general exponential must agree with the Taylor series
  Real const
  error_exp = norm(A - exp_taylor(L)) / norm(A);

  ASSERT_LE(error_exp, 16 * machine_epsilon<Real>());

  Tensor<Real> const
  R = log(A);

  Real const
  error_rt = norm(R - L) / norm(L);

  ASSERT_LE(error_rt, 256 * machine_epsilon<Real>());
}

TEST(MiniTensor, LogScaledRotation2x2)
{
  // 2x2 scaled rotation: log(s R(t)) = [log s, -t; t, log s] for |t| < pi.
  Real const s = 2.0;
  Real const t = 2.5;

  Tensor<Real> const
  A(s * std::cos(t), -s * std::sin(t),
    s * std::sin(t),  s * std::cos(t));

  Tensor<Real> const
  L = log(A);

  Tensor<Real> const
  L_exact(std::log(s), -t, t, std::log(s));

  Real const
  error = norm(L - L_exact) / norm(L_exact);

  ASSERT_LE(error, 16 * machine_epsilon<Real>());
}

TEST(MiniTensor, BakerCampbellHausdorff)
{
  Real const
  Pi = std::acos(-1.0);

  Real const
  gamma = 0.1;

  Tensor<Real> const u(0, gamma, 0, gamma, 0, 0, 0, 0, 0);

  Tensor<Real> const r(0, -Pi/4, 0, Pi/4, 0, 0, 0, 0, 0);

  Tensor<Real> const f_bch = bch(r, u);

  Tensor<Real> const U = exp(u);

  Tensor<Real> const R = exp(r);

  Tensor<Real> const F = R * U;

  Tensor<Real> const f = log(F);

  Real const
  error = norm(f_bch - f) / norm(F);

  // Our implementation of the Baker-Campbell-Hausdorff
  // formula uses only 4 terms, so we expect some error here.
  Real const
  tolerance = 1.0e-3;

  ASSERT_LE(error, tolerance);
}

TEST(MiniTensor, BinaryPowering)
{
  Tensor<Real> const A(1.1, 0.2, 0.3, 0.4, 1.2, 0.5, 0.6, 0.7, 1.3);

  // Qualified on purpose: binary_powering is part of the public
  // interface and external packages call it as minitensor::binary_powering.
  ASSERT_LE(norm(minitensor::binary_powering(A, 0) - eye<Real>(3)),
      machine_epsilon<Real>());

  Tensor<Real> B = eye<Real>(3);

  for (Index exponent = 1; exponent <= 9; ++exponent) {

    B = B * A;

    Real const
    error = norm(minitensor::binary_powering(A, exponent) - B) / norm(B);

    ASSERT_LE(error, 128.0 * machine_epsilon<Real>());

  }
}

} // namespace minitensor
