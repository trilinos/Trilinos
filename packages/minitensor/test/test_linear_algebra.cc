// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests for the linear-algebra modules (Norms, Inverse,
// Factorizations): invariants, inversion, eigendecompositions, SVD,
// polar decompositions and Cholesky, including the property-based
// (randomized) decomposition tests.

#include <random>

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

namespace {

// Number of random samples per dimension.
Index const NUM_SAMPLES = 256;

// Fixed seed: reproducible, order-independent failures.
unsigned const SEED = 20260624u;

// Reconstruction / invariant tolerances. A decomposition is a few O(1)
// operations or a tightly converged iteration, so the residuals sit at a small
// multiple of machine epsilon. The bugs these tests target (sign flips,
// eigenvector/eigenvalue mispairing) produce O(1) residuals, so the exact
// constant is not delicate; these are loose enough to never flake yet far
// below any real defect.
Real const TOL_DIRECT    = 1.0e-11;  // eig, inverse: near machine precision
Real const TOL_ITERATIVE = 1.0e-9;   // svd, polar: Newton/Jacobi iterations
Real const TOL_SERIES    = 1.0e-8;   // exp / log: truncated series

//
// Random tensor builders.
//

Tensor<Real>
random_general(Index const dimension, std::mt19937_64 & rng)
{
  std::normal_distribution<Real> normal(0.0, 1.0);
  Tensor<Real>                   A(dimension);
  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      A(i, j) = normal(rng);
    }
  }
  return A;
}

Tensor<Real>
random_symmetric(Index const dimension, std::mt19937_64 & rng)
{
  return sym(random_general(dimension, rng));
}

// Symmetric positive definite: G G^T is positive semidefinite; the shift makes
// it strictly positive definite and keeps the conditioning reasonable.
Tensor<Real>
random_spd(Index const dimension, std::mt19937_64 & rng)
{
  Tensor<Real> const G = random_general(dimension, rng);
  return G * transpose(G) + Real(dimension) * eye<Real>(dimension);
}

// Deformation-gradient-like tensor: a moderate perturbation of the identity,
// guaranteed invertible with positive determinant and good conditioning. This
// keeps the polar/inverse tests independent of the matrix exponential.
Tensor<Real>
random_deformation(Index const dimension, std::mt19937_64 & rng)
{
  return eye<Real>(dimension) + 0.25 * random_general(dimension, rng);
}

//
// Invariant checks (return a nonnegative error; 0 is perfect).
//

Real
orthonormality_error(Tensor<Real> const & Q)
{
  Index const dimension = Q.get_dimension();
  return norm(transpose(Q) * Q - eye<Real>(dimension));
}

Real
symmetry_error(Tensor<Real> const & A)
{
  return norm(skew(A)) / std::max(norm(A), machine_epsilon<Real>());
}

// Largest off-diagonal magnitude relative to the largest diagonal magnitude.
Real
off_diagonal_error(Tensor<Real> const & D)
{
  Index const dimension = D.get_dimension();
  Real        off       = 0.0;
  Real        diag      = machine_epsilon<Real>();
  for (Index i = 0; i < dimension; ++i) {
    diag = std::max(diag, std::abs(D(i, i)));
    for (Index j = 0; j < dimension; ++j) {
      if (i != j) off = std::max(off, std::abs(D(i, j)));
    }
  }
  return off / diag;
}

bool
is_descending(Tensor<Real> const & D)
{
  Index const dimension = D.get_dimension();
  for (Index i = 0; i + 1 < dimension; ++i) {
    if (D(i, i) < D(i + 1, i + 1)) return false;
  }
  return true;
}

}  // anonymous namespace

TEST(MiniTensor, TensorManipulation)
{
  Tensor<Real> A = eye<Real>(3);
  Tensor<Real> B(3);
  Tensor<Real> C(3);
  Vector<Real> u(3);

  A = 2.0 * A;
  A(1, 0) = A(0, 1) = 1.0;
  A(2, 1) = A(1, 2) = 1.0;

  B = inverse(A);

  C = A * B;

  ASSERT_LE(norm(C - eye<Real>(3)), machine_epsilon<Real>());

  Real I1_A = I1(A);
  Real I2_A = I2(A);
  Real I3_A = I3(A);

  u(0) = I1_A - 6;
  u(1) = I2_A - 10;
  u(2) = I3_A - 4;

  Real const error = norm(u);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, Determinant)
{
  // The determinant dispatches on the run-time dimension, and the
  // branches that a statically sized tensor can never take are discarded
  // at compile time. Exercise every branch that remains, both statically
  // and dynamically sized, against determinants known in closed form.

  Tensor<Real, 1>
  A1(Filler::ZEROS);

  A1(0, 0) = 3.0;

  ASSERT_LE(std::abs(det(A1) - 3.0), machine_epsilon<Real>());

  Tensor<Real, 2> const
  A2(1.0, 2.0, 3.0, 4.0);

  ASSERT_LE(std::abs(det(A2) + 2.0), 8 * machine_epsilon<Real>());

  Tensor<Real, 3> const
  A3(2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0);

  ASSERT_LE(std::abs(det(A3) - 4.0), 8 * machine_epsilon<Real>());

  // Dimensions above three go through the Laplace expansion, which
  // recurses on subtensors that have the same static dimension but a
  // smaller run-time dimension. L and U are triangular, so the
  // determinant of their product is the product of all their diagonal
  // entries, here 24 * 1/24 = 1.
  Index const
  N = 4;

  Tensor<Real, N>
  L(Filler::ZEROS);

  Tensor<Real, N>
  U(Filler::ZEROS);

  for (Index i = 0; i < N; ++i) {

    L(i, i) = static_cast<Real>(i + 1);

    U(i, i) = 1.0 / static_cast<Real>(i + 1);

    for (Index j = 0; j < i; ++j) {
      L(i, j) = 0.5 * static_cast<Real>(i + j);
      U(j, i) = 0.25 * static_cast<Real>(i - j);
    }

  }

  Tensor<Real, N> const
  A4 = L * U;

  ASSERT_LE(std::abs(det(A4) - 1.0), 64 * machine_epsilon<Real>());

  // Dynamically sized tensors keep every branch, and must agree with
  // their statically sized counterparts.
  for (Index dimension = 1; dimension <= N; ++dimension) {

    Tensor<Real>
    A(dimension);

    Tensor<Real, N> const
    B = A4;

    for (Index i = 0; i < dimension; ++i) {
      for (Index j = 0; j < dimension; ++j) {
        A(i, j) = B(i, j);
      }
    }

    Tensor<Real, N>
    C(Filler::ZEROS);

    C.set_dimension(dimension);

    for (Index i = 0; i < dimension; ++i) {
      for (Index j = 0; j < dimension; ++j) {
        C(i, j) = B(i, j);
      }
    }

    ASSERT_LE(std::abs(det(A) - det(C)), 64 * machine_epsilon<Real>());

  }
}

TEST(MiniTensor, Inverse2x2)
{
  Index const
  N = 2;

  Tensor<Real, N> const
  A = 2.0 * eye<Real, N>() + Tensor<Real, N>(Filler::RANDOM_UNIFORM);

  Tensor<Real, N> const
  B = inverse(A);

  Tensor<Real, N> const
  C = A * B;

  Real const
  error = norm(C - eye<Real, N>()) / norm(A);

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (N - 1) * machine_epsilon<Real>();

  ASSERT_LE(error, tolerance);
}

TEST(MiniTensor, Inverse3x3)
{
  Index const
  N = 3;

  Tensor<Real, N> const
  A = 2.0 * eye<Real, N>() + Tensor<Real, N>(Filler::RANDOM_UNIFORM);

  Tensor<Real, N> const
  B = inverse(A);

  Tensor<Real, N> const
  C = A * B;

  Real const
  error = norm(C - eye<Real, N>()) / norm(A);

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (N - 1) * machine_epsilon<Real>();

  ASSERT_LE(error, tolerance);
}

TEST(MiniTensor, InverseNxN)
{
  Index const
  N = 11;

  Tensor<Real, N> const
  A = 2.0 * eye<Real, N>() + Tensor<Real, N>(Filler::RANDOM_UNIFORM);

  Tensor<Real, N> const
  B = inverse(A);

  Tensor<Real, N> const
  C = A * B;

  Real const
  error = norm(C - eye<Real, N>()) / norm(A);

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (N - 1) * machine_epsilon<Real>();

  ASSERT_LE(error, tolerance);
}

TEST(MiniTensor, Inverse_4th_NxN)
{
  Index const
  N = 4;

  Tensor4<Real, N> const
  A = 2.0 * identity_1<Real, N>() + Tensor4<Real, N>(Filler::RANDOM_UNIFORM);

  Tensor4<Real, N> const
  B = inverse(A);

  Tensor4<Real, N> const
  C = dotdot(A, B);

  Real const
  error = norm_f(C - identity_1<Real, N>()) / norm_f(A);

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (2 * N - 1) * machine_epsilon<Real>();

  ASSERT_LE(error, tolerance);
}

TEST(MiniTensor, SymmetricEigen)
{
  Tensor<Real> A = eye<Real>(3);
  A(0, 1) = 0.1;
  A(1, 0) = 0.1;

  Tensor<Real> V(3);
  Tensor<Real> D(3);

  std::tie(V, D) = eig_sym(A);

  ASSERT_LE(std::abs(D(0, 0) - 1.1), machine_epsilon<Real>());
  ASSERT_LE(std::abs(D(1, 1) - 1.0), machine_epsilon<Real>());
  ASSERT_LE(std::abs(D(2, 2) - 0.9), machine_epsilon<Real>());
}

TEST(MiniTensor, LeftPolarDecomposition)
{
  Tensor<Real> const X(1.1, 0.2, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0, 1.2);

  Real const
  c = std::sqrt(2.0) / 2.0;

  Tensor<Real> const Y(c, -c, 0.0, c, c, 0.0, 0.0, 0.0, 1.0);

  Tensor<Real> const F = X * Y;
  Tensor<Real> V(3);
  Tensor<Real> R(3);

  std::tie(V, R) = polar_left(F);

  Real const
  error_x = norm(V - X) / norm(X);

  Real const
  error_y = norm(R - Y) / norm(Y);

  ASSERT_LE(error_x, machine_epsilon<Real>());
  ASSERT_LE(error_y, machine_epsilon<Real>());
}

TEST(MiniTensor, PolarLeftLog)
{
  Real const
  gamma = 0.1;

  Tensor<Real> const x(0, gamma, 0, gamma, 0, 0, 0, 0, 0);

  Tensor<Real> const X = exp(x);

  Real const
  c = std::sqrt(2.0) / 2.0;

  Tensor<Real> const Y(c, -c, 0.0, c, c, 0.0, 0.0, 0.0, 1.0);

  Tensor<Real> const F = X * Y;

  Tensor<Real> V(3), R(3), v(3);

  std::tie(V, R, v) = polar_left_logV(F);

  Real const error = norm(v - x) / norm(x);

  Real const
  tolerance = 16.0 * machine_epsilon<Real>();

  ASSERT_LE(error, tolerance);
}

TEST(MiniTensor, SVD2x2)
{
  Real const phi = 1.0;

  Real const psi = 2.0;

  Real const s0 = std::sqrt(3.0);

  Real const s1 = std::sqrt(2.0);

  Real const cl = cos(phi);

  Real const sl = sin(phi);

  Real const cr = cos(psi);

  Real const sr = sin(psi);

  Tensor<Real> const X(cl, -sl, sl, cl);

  Tensor<Real> const Y(cr, -sr, sr, cr);

  Tensor<Real> const D(s0, 0.0, 0.0, s1);

  Tensor<Real> const A = X * D * transpose(Y);

  Tensor<Real> U(2), S(2), V(2);

  std::tie(U, S, V) = svd(A);

  Tensor<Real> B = U * S * transpose(V);

  Real const error = norm(A - B) / norm(A);

  ASSERT_LE(error, 2.0 * machine_epsilon<Real>());
}

TEST(MiniTensor, SVD3x3)
{
  Tensor<Real> const A(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

  Tensor<Real> U(3), S(3), V(3);

  std::tie(U, S, V) = svd(A);

  Tensor<Real> const B = U * S * transpose(V);

  Real const error = norm(A - B) / norm(A);

  ASSERT_LE(error, 2.0 * machine_epsilon<Real>());
}

TEST(MiniTensor, SVD3x3Fad)
{
  Tensor<Sacado::Fad::DFad<Real>> const
  A(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

  Tensor<Sacado::Fad::DFad<Real>> U(3), S(3), V(3);

  std::tie(U, S, V) = svd(A);

  Tensor<Sacado::Fad::DFad<Real>> const
  B = U * S * transpose(V);

  Sacado::Fad::DFad<Real> const
  error = norm(B - A) / norm(A);

  ASSERT_LE(error, 2.0 * machine_epsilon<Real>());
}

TEST(MiniTensor, SymmetricEigen2x2)
{
  // Reconstruction A = V * D * transpose(V) must hold for general symmetric
  // 2x2 tensors, including unequal diagonals with shear. The unequal-diagonal
  // cases below regress Trilinos issue #15389, where the 2x2 path returned an
  // inconsistent (V, D) that flipped the off-diagonal sign on reconstruction.
  std::vector<Tensor<Real>> const tensors = {
      Tensor<Real>(2.0, 1.0, 1.0, 2.0),        // equal diagonals
      Tensor<Real>(1.0025, 0.1, 0.1, 1.0025),  // equal diagonals, shear
      Tensor<Real>(1.168, 0.06, 0.06, 0.922),  // f > h, shear
      Tensor<Real>(0.922, 0.06, 0.06, 1.168),  // f < h, shear
      Tensor<Real>(2.0, 0.0, 0.0, 3.0),        // diagonal (g == 0)
      Tensor<Real>(-1.0, 0.5, 0.5, -2.0)       // negative eigenvalues
  };

  for (Tensor<Real> const & A : tensors) {
    Tensor<Real> V(2), D(2);

    std::tie(V, D) = eig_sym(A);

    Tensor<Real> const B = V * D * transpose(V);

    Real const error = norm(A - B) / norm(A);

    ASSERT_LE(error, 2.0 * machine_epsilon<Real>());

    // Eigenvalues are returned in descending order (matches eig_sym_NxN).
    ASSERT_GE(D(0, 0), D(1, 1));
  }
}

TEST(MiniTensor, SymmetricEigen3x3)
{
  Tensor<Real> const A(2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0);

  Tensor<Real> V(3), D(3);

  std::tie(V, D) = eig_sym(A);

  Tensor<Real> const B = V * D * transpose(V);

  Real const error = norm(A - B) / norm(A);

  ASSERT_LE(error, 4.0 * machine_epsilon<Real>());
}

TEST(MiniTensor, Polar3x3)
{
  Tensor<Real> const A(2.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0);

  Tensor<Real> R(3), U(3);

  std::tie(R, U) = polar_right(A);

  Tensor<Real> X(3), D(3), Y(3);

  std::tie(X, D, Y) = svd(A);

  Tensor<Real> const B = R - X * transpose(Y) + U - Y * D * transpose(Y);

  Real const error = norm(B) / norm(A);

  ASSERT_LE(error, 8.0 * machine_epsilon<Real>());
}

TEST(MiniTensor, Cholesky)
{
  Tensor<Real> const A(1.0, 1.0, 1.0, 1.0, 5.0, 3.0, 1.0, 3.0, 3.0);

  Tensor<Real> G(3);

  bool is_spd;

  std::tie(G, is_spd) = cholesky(A);

  Tensor<Real> const B(1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 1.0, 1.0);

  Real const error = norm(G - B) / norm(A);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, EigSymProperties)
{
  for (Index const dimension : {2, 3, 4}) {
    std::mt19937_64 rng(SEED);
    for (Index sample = 0; sample < NUM_SAMPLES; ++sample) {
      Tensor<Real> const A = random_symmetric(dimension, rng);

      Tensor<Real> V(dimension), D(dimension);
      std::tie(V, D) = eig_sym(A);

      Real const scale = std::max(norm(A), machine_epsilon<Real>());

      ASSERT_LE(norm(A - V * D * transpose(V)) / scale, TOL_DIRECT)
          << "eig_sym reconstruction, dim " << dimension << ", sample " << sample;
      ASSERT_LE(orthonormality_error(V), TOL_DIRECT)
          << "eig_sym eigenvectors not orthonormal, dim " << dimension;
      ASSERT_LE(off_diagonal_error(D), TOL_DIRECT)
          << "eig_sym D not diagonal, dim " << dimension;
      ASSERT_TRUE(is_descending(D))
          << "eig_sym eigenvalues not descending, dim " << dimension;
    }
  }
}

TEST(MiniTensor, EigSpdProperties)
{
  for (Index const dimension : {2, 3, 4}) {
    std::mt19937_64 rng(SEED);
    for (Index sample = 0; sample < NUM_SAMPLES; ++sample) {
      Tensor<Real> const A = random_spd(dimension, rng);

      Tensor<Real> V(dimension), D(dimension);
      std::tie(V, D) = eig_spd(A);

      Real const scale = std::max(norm(A), machine_epsilon<Real>());

      ASSERT_LE(norm(A - V * D * transpose(V)) / scale, TOL_DIRECT)
          << "eig_spd reconstruction, dim " << dimension << ", sample " << sample;
      ASSERT_LE(orthonormality_error(V), TOL_DIRECT)
          << "eig_spd eigenvectors not orthonormal, dim " << dimension;
      for (Index i = 0; i < dimension; ++i) {
        ASSERT_GT(D(i, i), 0.0)
            << "eig_spd eigenvalue not positive, dim " << dimension;
      }
    }
  }
}

TEST(MiniTensor, SvdProperties)
{
  for (Index const dimension : {2, 3, 4}) {
    std::mt19937_64 rng(SEED);
    for (Index sample = 0; sample < NUM_SAMPLES; ++sample) {
      Tensor<Real> const A = random_general(dimension, rng);

      Tensor<Real> U(dimension), S(dimension), V(dimension);
      std::tie(U, S, V) = svd(A);

      Real const scale = std::max(norm(A), machine_epsilon<Real>());

      ASSERT_LE(norm(A - U * S * transpose(V)) / scale, TOL_ITERATIVE)
          << "svd reconstruction, dim " << dimension << ", sample " << sample;
      ASSERT_LE(orthonormality_error(U), TOL_ITERATIVE)
          << "svd U not orthonormal, dim " << dimension;
      ASSERT_LE(orthonormality_error(V), TOL_ITERATIVE)
          << "svd V not orthonormal, dim " << dimension;
      ASSERT_LE(off_diagonal_error(S), TOL_ITERATIVE)
          << "svd S not diagonal, dim " << dimension;
      ASSERT_TRUE(is_descending(S))
          << "svd singular values not descending, dim " << dimension << ", sample " << sample;
      for (Index i = 0; i < dimension; ++i) {
        ASSERT_GE(S(i, i), -TOL_ITERATIVE)
            << "svd singular value negative, dim " << dimension << ", sample " << sample;
      }
    }
  }
}

TEST(MiniTensor, PolarLeftProperties)
{
  for (Index const dimension : {2, 3}) {
    std::mt19937_64 rng(SEED);
    for (Index sample = 0; sample < NUM_SAMPLES; ++sample) {
      Tensor<Real> const A = random_deformation(dimension, rng);

      Tensor<Real> V(dimension), R(dimension);
      std::tie(V, R) = polar_left(A);

      Real const scale = std::max(norm(A), machine_epsilon<Real>());

      ASSERT_LE(norm(A - V * R) / scale, TOL_ITERATIVE)
          << "polar_left reconstruction, dim " << dimension << ", sample " << sample;
      ASSERT_LE(orthonormality_error(R), TOL_ITERATIVE)
          << "polar_left R not orthogonal, dim " << dimension;
      ASSERT_NEAR(det(R), 1.0, TOL_ITERATIVE)
          << "polar_left R not a proper rotation, dim " << dimension;
      ASSERT_LE(symmetry_error(V), TOL_ITERATIVE)
          << "polar_left V not symmetric, dim " << dimension;
    }
  }
}

TEST(MiniTensor, PolarRightProperties)
{
  for (Index const dimension : {2, 3}) {
    std::mt19937_64 rng(SEED);
    for (Index sample = 0; sample < NUM_SAMPLES; ++sample) {
      Tensor<Real> const A = random_deformation(dimension, rng);

      Tensor<Real> R(dimension), U(dimension);
      std::tie(R, U) = polar_right(A);

      Real const scale = std::max(norm(A), machine_epsilon<Real>());

      ASSERT_LE(norm(A - R * U) / scale, TOL_ITERATIVE)
          << "polar_right reconstruction, dim " << dimension << ", sample " << sample;
      ASSERT_LE(orthonormality_error(R), TOL_ITERATIVE)
          << "polar_right R not orthogonal, dim " << dimension;
      ASSERT_NEAR(det(R), 1.0, TOL_ITERATIVE)
          << "polar_right R not a proper rotation, dim " << dimension;
      ASSERT_LE(symmetry_error(U), TOL_ITERATIVE)
          << "polar_right U not symmetric, dim " << dimension;
    }
  }
}

TEST(MiniTensor, ExpLogInverseProperties)
{
  for (Index const dimension : {2, 3}) {
    std::mt19937_64 rng(SEED);
    for (Index sample = 0; sample < NUM_SAMPLES; ++sample) {
      Tensor<Real> const A = random_spd(dimension, rng);

      Real const scale_a = std::max(norm(A), machine_epsilon<Real>());
      ASSERT_LE(norm(exp(log(A)) - A) / scale_a, TOL_SERIES)
          << "exp(log(A)) != A, dim " << dimension << ", sample " << sample;

      Tensor<Real> const B = 0.5 * random_symmetric(dimension, rng);

      Real const scale_b = std::max(norm(B), machine_epsilon<Real>());
      ASSERT_LE(norm(log(exp(B)) - B) / scale_b, TOL_SERIES)
          << "log(exp(B)) != B, dim " << dimension << ", sample " << sample;
    }
  }
}

TEST(MiniTensor, InverseProperties)
{
  for (Index const dimension : {2, 3, 4}) {
    std::mt19937_64 rng(SEED);
    for (Index sample = 0; sample < NUM_SAMPLES; ++sample) {
      Tensor<Real> const A     = random_deformation(dimension, rng);
      Tensor<Real> const A_inv = inverse(A);
      Tensor<Real> const I     = eye<Real>(dimension);

      ASSERT_LE(norm(A * A_inv - I), TOL_DIRECT)
          << "A * inverse(A) != I, dim " << dimension << ", sample " << sample;
      ASSERT_LE(norm(A_inv * A - I), TOL_DIRECT)
          << "inverse(A) * A != I, dim " << dimension << ", sample " << sample;
    }
  }
}

} // namespace minitensor
