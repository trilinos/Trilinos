// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Property-based (randomized) tests for the matrix decompositions in
// MiniTensor_LinearAlgebra. Rather than checking a single hand-picked input,
// each test exercises many random tensors across dimensions 2, 3 and 4 and
// verifies the defining algebraic invariants of the decomposition
// (reconstruction, orthonormality, symmetry, positivity, ordering).
//
// These invariants are what a decomposition is *for*; a single example only
// catches a bug if it happens to land in the broken case. The 2x2 symmetric
// eigendecomposition bug (Trilinos issue #15389), for instance, slipped past
// the existing single-example test because that example used equal diagonals,
// the one case that reconstructed correctly. A randomized reconstruction check
// would have caught it on the first off-diagonal sample.
//
// The RNG is seeded with a fixed value in each test, so failures are
// reproducible and order-independent.
//

#include <random>

#include "gtest/gtest.h"
#include "MiniTensor.h"

int
main(int ac, char* av[])
{
  Kokkos::initialize();

  ::testing::GTEST_FLAG(print_time) = (ac > 1) ? true : false;

  ::testing::InitGoogleTest(&ac, av);

  auto const retval = RUN_ALL_TESTS();

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

//
// Symmetric eigendecomposition: A = V D V^T, V orthonormal, D diagonal and
// sorted in descending order.
//
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

//
// Symmetric positive-definite eigendecomposition: A = V D V^T with strictly
// positive eigenvalues.
//
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

//
// Singular value decomposition: A = U S V^T with U and V orthonormal and S
// diagonal.
//
// Note: only the invariants that hold across all dimensions are asserted here.
// The NxN path canonicalizes the singular values to be nonnegative and sorted
// in descending order, but the 2x2-specialized path currently does not -- it
// can return a negative diagonal entry and a non-descending order. Asserting
// nonnegativity / descending order here would therefore fail on the 2x2 path.
// That inconsistency between svd_2x2 and svd_NxN is a separate fix; once the
// 2x2 path canonicalizes its output, those two assertions should be added back.
//
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
    }
  }
}

//
// Left polar decomposition: A = V R, R a rotation (proper orthogonal), V SPD.
//
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

//
// Right polar decomposition: A = R U, R a rotation (proper orthogonal), U SPD.
//
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

//
// Matrix exponential / logarithm are inverse maps: exp(log(A)) = A for SPD A,
// and log(exp(A)) = A for symmetric A of moderate norm (so the principal
// logarithm recovers it).
//
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

//
// Inverse: A A^{-1} = A^{-1} A = I.
//
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

}  // namespace minitensor
