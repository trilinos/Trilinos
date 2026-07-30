// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Factorizations_h)
#define MiniTensor_Factorizations_h

// Matrix factorizations: eigen, SVD, polar, Cholesky, Givens.
#include "MiniTensor_Inverse.h"
#include "MiniTensor_Norms.h"

namespace minitensor {

///
/// Apply Givens-Jacobi rotation on the left in place.
/// \param c and s for a rotation G in form [c, s; -s, c]
/// \param i and k indices for rows and columns where rotation is applied.
/// \param A tensor to rotate
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
givens_left(T const & c, T const & s, Index i, Index k, Tensor<T, N> & A);

///
/// Apply Givens-Jacobi rotation on the right in place.
/// \param c and s for a rotation G in form [c, s; -s, c]
/// \param i and k indices for rows and columns where rotation is applied.
/// \param A tensor to rotate
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
givens_right(T const & c, T const & s, Index i, Index k, Tensor<T, N> & A);

///
/// Sort and index. Useful for ordering singular values
/// and eigenvalues and corresponding vectors in the
/// respective decompositions.
/// \param u vector to sort
/// \return v P sorted vector, permutation matrix such that v = P^T u
///
template <typename T, Index N>
std::pair<Vector<T, N>, Tensor<T, N>> sort_permutation(Vector<T, N> const &u);

///
/// Singular value decomposition (SVD)
/// \return \f$ A = USV^T\f$
///
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>> svd(Tensor<T, N> const &A);

///
/// Project to O(N) (Orthogonal Group) using a Newton-type algorithm.
/// See Higham's Functions of Matrices p210 [2008]
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ R = \arg min_Q \|A - Q\|\f$
/// This algorithm projects a given tensor in GL(N) to O(N).
/// The rotation/reflection obtained through this projection is
/// the orthogonal component of the real polar decomposition
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
polar_rotation(Tensor<T, N> const & A);

///
/// Left polar decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_left(Tensor<T, N> const &A);

///
/// Right polar decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ RU = A \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_right(Tensor<T, N> const &A);

///
/// Left polar decomposition computed with eigenvalue decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_left_eig(Tensor<T, N> const &A);

///
/// R^3 right polar decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ RU = F \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_right_eig(Tensor<T, N> const &A);

///
/// Left polar decomposition with matrix logarithm for V
/// \param F tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD, and log V
///
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
polar_left_logV(Tensor<T, N> const &F);

template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
polar_left_logV_eig(Tensor<T, N> const &F);

///
/// Left polar decomposition with matrix logarithm for V using eig_spd_cos
/// \param F tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD, and log V
///
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
polar_left_logV_lame(Tensor<T, N> const &F);

///
/// Symmetric Schur algorithm for R^2.
/// \param \f$ A = [f, g; g, h] \in S(2) \f$
/// \return \f$ c, s \rightarrow [c, -s; s, c]\f$ diagonalizes A$
///
template <typename T>
std::pair<T, T> schur_sym(const T f, const T g, const T h);

///
/// Givens rotation. [c, -s; s, c] [a; b] = [r; 0]
/// \return c and s
///
template <typename T> std::pair<T, T> givens(T const &a, T const &b);

///
/// Eigenvalue decomposition for symmetric 2nd-order tensor
/// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_sym(Tensor<T, N> const &A);

///
/// Eigenvalue decomposition for SPD 2nd-order tensor
/// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_spd(Tensor<T, N> const &A);

///
/// Eigenvalue decomposition for SPD 2nd-order tensor
/// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
/// This algorithm comes from the journal article
/// Scherzinger and Dohrmann, CMAME 197 (2008) 4007-4015
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_spd_cos(Tensor<T, N> const &A);

///
/// Cholesky decomposition, rank-1 update algorithm
/// (Matrix Computations 3rd ed., Golub & Van Loan, p145)
/// \param A assumed symmetric tensor
/// \return G Cholesky factor A = GG^T and completed (bool)
/// algorithm ran to completion
///
template <typename T, Index N>
std::pair<Tensor<T, N>, bool> cholesky(Tensor<T, N> const &A);

///
/// Condition number: ratio of largest to smalest singular values.
///
template <typename T, Index N> T cond(Tensor<T, N> const &A);

///
/// Reciprocal condition number: ratio of smallest to largest singular values.
///
template <typename T, Index N> T inv_cond(Tensor<T, N> const &A);

//
// Condition number.
//
template <typename T, Index N> T cond(Tensor<T, N> const &A) {
  Index const
  dimension = A.get_dimension();

  Tensor<T, N> const S = std::get<1>(svd(A));

  T const
  k = S(0, 0) / S(dimension - 1, dimension - 1);

  return k;
}

//
// Reciprocal condition number.
//
template <typename T, Index N> T inv_cond(Tensor<T, N> const &A) {
  Index const
  dimension = A.get_dimension();

  Tensor<T, N> const S = std::get<1>(svd(A));

  T const
  k = S(dimension - 1, dimension - 1) / S(0, 0);

  return k;
}

//
// Sort and index in descending order. Useful for ordering singular values
// and eigenvalues and corresponding vectors in the respective decompositions.
//
template <typename T, Index N>
std::pair<Vector<T, N>, Tensor<T, N>> sort_permutation(Vector<T, N> const &u) {

  Index const
  dimension = u.get_dimension();

  std::vector <std::pair<T, Index>>
  s(dimension);

  for (Index i = 0; i < dimension; ++i) {
    s[i].first = u(i);
    s[i].second = i;
  }

  std::sort(s.begin(), s.end(), greater_than<std::pair<T, Index>>);

  Vector<T, N> v(dimension);

  Tensor<T, N>
  P = zero<T, N>(dimension);

  for (Index i = 0; i < dimension; ++i) {
    v(i) = s[i].first;
    P(s[i].second, i) = 1.0;
  }

  return std::make_pair(v, P);
}

//
// Apply Givens-Jacobi rotation on the left in place.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
givens_left(T const & c, T const & s, Index i, Index k, Tensor<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index j = 0; j < dimension; ++j) {
    T const t1 = A(i,j);
    T const t2 = A(k,j);
    A(i,j) = c * t1 - s * t2;
    A(k,j) = s * t1 + c * t2;
  }
  return;
}

//
// Apply Givens-Jacobi rotation on the right in place.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
givens_right(T const & c, T const & s, Index i, Index k, Tensor<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index j = 0; j < dimension; ++j) {
    T const t1 = A(j,i);
    T const t2 = A(j,k);
    A(j,i) = c * t1 - s * t2;
    A(j,k) = s * t1 + c * t2;
  }
  return;
}

namespace {

//
// Singular value decomposition (SVD) for 2x2
// bidiagonal matrix. Used for general 2x2 SVD.
// Adapted from LAPAPCK's DLASV2, Netlib's dlasv2.c
// and LBNL computational crystallography toolbox
// \param f, g, h where A = [f, g; 0, h]
// \return \f$ A = USV^T\f$
//
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>> svd_bidiagonal(T f, T g,
                                                                    T h) {
  T fa = std::abs(f);
  T ga = std::abs(g);
  T ha = std::abs(h);

  T s0 = 0.0;
  T s1 = 0.0;

  T cu = 1.0;
  T su = 0.0;
  T cv = 1.0;
  T sv = 0.0;

  bool swap_diag = (ha > fa);

  if (swap_diag == true) {
    std::swap(fa, ha);
    std::swap(f, h);
  }

  // diagonal matrix
  if (ga == 0.0) {
    s1 = ha;
    s0 = fa;
  } else if (ga > fa && fa / ga < machine_epsilon<T>()) {
    // case of very large ga
    s0 = ga;
    s1 = ha > 1.0 ?
        T(fa / (ga / ha)) :
        T((fa / ga) * ha);
    cu = 1.0;
    su = h / g;
    cv = f / g;
    sv = 1.0;
  } else {
    // normal case
    T d = fa - ha;
    T l = d / fa; // l \in [0,1]
    T m = g / f; // m \in (-1/macheps, 1/macheps)
    T t = 2.0 - l; // t \in [1,2]
    T mm = m * m;
    T tt = t * t;
    T s = std::sqrt(tt + mm); // s \in [1,1 + 1/macheps]
    T r = l != 0.0 ?
        T(std::sqrt(l * l + mm)) :
        T(std::abs(m)); // r \in [0,1 + 1/macheps]
    T a = 0.5 * (s + r); // a \in [1,1 + |m|]
    s1 = ha / a;
    s0 = fa * a;

    // Compute singular vectors
    T tau; // second assignment to T in DLASV2
    if (mm != 0.0) {
      tau = (m / (s + t) + m / (r + l)) * (1.0 + a);
    } else {
      // note that m is very tiny
      tau = l == 0.0 ?
          T(copysign(T(2.0), f) * copysign(T(1.0), g)) :
          T(g / copysign(d, f) + m / t);
    }
    T lv = std::sqrt(tau * tau + 4.0); // second assignment to L in DLASV2
    cv = 2.0 / lv;
    sv = tau / lv;
    cu = (cv + sv * m) / a;
    su = (h / f) * sv / a;
  }

  // Fix signs of singular values in accordance to sign of singular vectors
  s0 = copysign(s0, f);
  s1 = copysign(s1, h);

  if (swap_diag == true) {
    std::swap(cu, sv);
    std::swap(su, cv);
  }

  Tensor<T, N> U(cu, -su, su, cu);
  Tensor<T, N> S(s0, 0.0, 0.0, s1);
  Tensor<T, N> V(cv, -sv, sv, cv);

  return std::make_tuple(U, S, V);
}

//
// R^2 singular value decomposition (SVD)
// \param A tensor
// \return \f$ A = USV^T\f$
//
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
svd_2x2(Tensor<T, N> const &A) {
  assert(A.get_dimension() == 2);

  // First compute a givens rotation to eliminate 1,0 entry in tensor
  T c = 1.0;
  T s = 0.0;
  std::tie(c, s) = givens(A(0, 0), A(1, 0));

  Tensor<T, N>
  R(c, -s, s, c);

  Tensor<T, N>
  B = R * A;

  // B is bidiagonal. Use specialized algorithm to compute its SVD
  Tensor<T, N>
  X(2), S(2), V(2);

  std::tie(X, S, V) = svd_bidiagonal<T, N>(B(0, 0), B(0, 1), B(1, 1));

  // Complete general 2x2 SVD with givens rotation calculated above
  Tensor<T, N>
  U = transpose(R) * X;

  return std::make_tuple(U, S, V);
}

//
// R^N singular value decomposition (SVD)
// \param A tensor
// \return \f$ A = USV^T\f$
//
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
svd_NxN(Tensor<T, N> const &A) {
  // Scale first
  T const
  norm_a = norm(A);

  T const
  scale = norm_a > 0.0 ? norm_a : T(1.0);

  Tensor<T, N>
  S = A / scale;

  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  U = identity<T, N>(dimension);

  Tensor<T, N>
  V = identity<T, N>(dimension);

  T
  off = norm_off_diagonal(S);

  T const
  tol = machine_epsilon<T>();

  Index const
  max_iter = 2048;

  Index
  num_iter = 0;

  while (off > tol && num_iter < max_iter) {

    // Find largest off-diagonal entry
    Index
    p = 0;

    Index
    q = 0;

    std::tie(p, q) = arg_max_off_diagonal(S);

    if (p > q) {
      std::swap(p, q);
    }

    // Obtain left and right Givens rotations by using 2x2 SVD
    Tensor <T, 2>
    Spq(S(p,p), S(p,q), S(q,p), S(q,q));

    Tensor <T, 2>
    L(2), D(2), R(2);

    std::tie(L, D, R) = svd_2x2(Spq);

    T const &
    cl = L(0,0);

    T const &
    sl = L(0,1);

    T const &
    cr = R(0,0);

    T const &
    sr = (sgn(R(0,1)) == sgn(R(1,0))) ? T(-R(0,1)) : T(R(0,1));

    // Apply both Givens rotations to matrices
    // that are converging to singular values and singular vectors
    givens_left(cl, sl, p, q, S);
    givens_right(cr, sr, p, q, S);

    givens_right(cl, sl, p, q, U);
    givens_left(cr, sr, p, q, V);

    off = norm_off_diagonal(S);
    num_iter++;
  }

  if (num_iter == max_iter) {
    MT_WARNING("SVD iteration did not converge.");
  }

  // Fix signs for entries in the diagonal matrix S
  // that are negative
  for (Index i = 0; i < dimension; ++i) {
    if (S(i,i) < 0.0) {
      S(i,i) = -S(i,i);
      for (Index j = 0; j < dimension; ++j) {
        U(j,i) = -U(j,i);
      }
    }
  }

  Vector<T, N> s(dimension);
  Tensor<T, N> P(dimension);

  std::tie(s, P) = sort_permutation(diag(S));
  S = scale * diag(s);
  U = U * P;
  V = V * P;

  return std::make_tuple(U, diag(diag(S)), transpose(V));
}

} // anonymous namespace

//
// R^N singular value decomposition (SVD)
// \param A tensor
// \return \f$ A = USV^T\f$
//
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
svd(Tensor<T, N> const &A) {
  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  U(dimension), S(dimension), V(dimension);

  switch (dimension) {

    default:
      std::tie(U, S, V) = svd_NxN(A);
      break;

    case 2:
      std::tie(U, S, V) = svd_2x2(A);
      // svd_2x2 doubles as a building block inside svd_NxN's Jacobi sweep, so
      // it returns the raw 2x2 factorization without the sign/order
      // canonicalization that svd_NxN applies to its own result. Apply that
      // canonicalization here, for the top-level 2D SVD only, so svd() returns
      // the standard convention (nonnegative singular values in descending
      // order) regardless of dimension. Folding a negative sign into the
      // corresponding column of U preserves A = U S V^T.
      for (Index i = 0; i < dimension; ++i) {
        if (S(i, i) < 0.0) {
          S(i, i) = -S(i, i);
          for (Index j = 0; j < dimension; ++j) {
            U(j, i) = -U(j, i);
          }
        }
      }
      // Sort descending. For 2x2 this is a single compare-and-swap of the two
      // columns of U and V (and the two singular values), avoiding the heap
      // allocation and matrix multiply that the general sort_permutation
      // incurs.
      if (S(0, 0) < S(1, 1)) {
        std::swap(S(0, 0), S(1, 1));
        std::swap(U(0, 0), U(0, 1));
        std::swap(U(1, 0), U(1, 1));
        std::swap(V(0, 0), V(0, 1));
        std::swap(V(1, 0), V(1, 1));
      }
      break;

  }

  return std::make_tuple(U, S, V);
}

//
// Project to O(N) (Orthogonal Group) using a Newton-type algorithm.
// See Higham's Functions of Matrices p210 [2008]
// \param A tensor (often a deformation-gradient-like tensor)
// \return \f$ R = \argmin_Q \|A - Q\|\f$
// This algorithm projects a given tensor in GL(N) to O(N).
// The rotation/reflection obtained through this projection is
// the orthogonal component of the real polar decomposition
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
polar_rotation(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  bool
  scale = true;

  T const
  tol_scale = 0.01;

  T const tol_conv =
      static_cast<Index>(Kokkos::sqrt(static_cast<double>(dimension))) *
      machine_epsilon<T>();

  Tensor<T, N>
  X = A;

  T
  gamma = 2.0;

  Index const
  max_iter = 128;

  Index
  num_iter = 0;

  while (num_iter < max_iter) {

    Tensor<T, N>
    Y = inverse(X);

    T
    mu = 1.0;

    if (scale == true) {
      mu = (norm_1(Y) * norm_infinity(Y)) / (norm_1(X) * norm_infinity(X));
      mu = std::sqrt(std::sqrt(mu));
    }

    Tensor<T, N>
    Z = 0.5 * (mu * X + transpose(Y) / mu);

    Tensor<T, N>
    D = Z - X;

    T
    delta = norm(D) / norm(Z);

    if (scale == true && delta < tol_scale) {
      scale = false;
    }

    bool
    end_iter =
        norm(D) <= std::sqrt(tol_conv) ||
        (delta > 0.5 * gamma && scale == false);

    X = Z;
    gamma = delta;

    if (end_iter == true) {
      break;
    }

    num_iter++;

  }

  if (num_iter == max_iter) {
    MT_WARNING("Polar iteration did not converge.");
  }

  return X;
}

//
// R^N Left polar decomposition
// \param A tensor (often a deformation-gradient-like tensor)
// \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_left(Tensor<T, N> const &A) {
  Tensor<T, N>
  R = polar_rotation(A);

  Tensor<T, N>
  V = sym(A * transpose(R));

  return std::make_pair(V, R);
}

//
// R^N Right polar decomposition
// \param A tensor (often a deformation-gradient-like tensor)
// \return \f$ RU = A \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_right(Tensor<T, N> const &A) {
  Tensor<T, N>
  R = polar_rotation(A);

  Tensor<T, N>
  U = sym(transpose(R) * A);

  return std::make_pair(R, U);
}

//
// R^3 left polar decomposition with eigenvalue decomposition
// \param F tensor (often a deformation-gradient-like tensor)
// \return \f$ VR = F \f$ with \f$ R \in SO(3) \f$ and V SPD(3)
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_left_eig(Tensor<T, N> const &F) {
  assert(F.get_dimension() == 3);

  // set up return tensors
  Tensor<T, N>
  R(3);

  Tensor<T, N>
  V(3);

  // temporary tensor used to compute R
  Tensor<T, N>
  Vinv(3);

  // compute spd tensor
  Tensor<T, N>
  b = F * transpose(F);

  // get eigenvalues/eigenvectors
  Tensor<T, N>
  eVal(3);

  Tensor<T, N>
  eVec(3);
  std::tie(eVec, eVal) = eig_spd(b);

  // compute sqrt() and inv(sqrt()) of eigenvalues
  Tensor<T, N>
  x = zero<T, N>(3);

  for (Index i = 0; i < 3; ++i) {
    if (eVal(i, i) < T(0)) {
      MT_ERROR_EXIT("Non-SPD input: negative eigenvalue.");
    }
  }

  x(0,0) = std::sqrt(eVal(0,0));
  x(1,1) = std::sqrt(eVal(1,1));
  x(2,2) = std::sqrt(eVal(2,2));

  Tensor<T, N>
  xi = zero<T, N>(3);

  xi(0,0) = 1.0 / x(0,0);
  xi(1,1) = 1.0 / x(1,1);
  xi(2,2) = 1.0 / x(2,2);

  // compute V, Vinv, and R
  V    = eVec * x * transpose(eVec);
  Vinv = eVec * xi * transpose(eVec);
  R    = Vinv * F;
  return std::make_pair(V, R);
}

//
// R^3 right polar decomposition with eigenvalue decomposition
// \param F tensor (often a deformation-gradient-like tensor)
// \return \f$ RU = F \f$ with \f$ R \in SO(3) \f$ and U SPD(3)
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_right_eig(Tensor<T, N> const &F) {
  Index const
  dimension = F.get_dimension();

  assert(dimension == 3);

  Tensor<T, N>
  R(dimension);

  Tensor<T, N>
  U(dimension);

  // temporary tensor used to compute R
  Tensor<T, N>
  Uinv(dimension);

  // compute spd tensor
  Tensor<T, N>
  C = transpose(F) * F;

  // get eigenvalues/eigenvectors
  Tensor<T, N>
  eVal(dimension);

  Tensor<T, N>
  eVec(dimension);

  std::tie(eVec, eVal) = eig_spd(C);

  // compute sqrt() and inv(sqrt()) of eigenvalues
  Tensor<T, N>
  x = zero<T, N>(dimension);

  for (Index i = 0; i < dimension; ++i) {
    if (eVal(i, i) < T(0)) {
      MT_ERROR_EXIT("Non-SPD input: negative eigenvalue.");
    }
  }

  x(0,0) = std::sqrt(eVal(0,0));
  x(1,1) = std::sqrt(eVal(1,1));
  x(2,2) = std::sqrt(eVal(2,2));

  Tensor<T, N>
  xi = zero<T, N>(dimension);

  xi(0,0) = 1.0 / x(0,0);
  xi(1,1) = 1.0 / x(1,1);
  xi(2,2) = 1.0 / x(2,2);

  // compute U, Uinv, and R
  U    = eVec * x * transpose(eVec);
  Uinv = eVec * xi * transpose(eVec);
  R    = F * Uinv;

  return std::make_pair(R, U);
}

//
// R^N left polar decomposition with matrix logarithm for V
// \param F tensor (often a deformation-gradient-like tensor)
// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD(N), and log V
//
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
polar_left_logV(Tensor<T, N> const &F) {
  Index const
  dimension = F.get_dimension();

  Tensor<T, N>
  X(dimension), S(dimension), Y(dimension);

  std::tie(X, S, Y) = svd(F);

  Tensor<T, N>
  R = X * transpose(Y);

  Tensor<T, N>
  V = X * S * transpose(X);

  Tensor<T, N>
  s = S;

  for (Index i = 0; i < dimension; ++i) {
    s(i,i) = std::log(s(i,i));
  }

  Tensor<T, N>
  v = X * s * transpose(X);

  return std::make_tuple(V, R, v);
}

template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
polar_left_logV_eig(Tensor<T, N> const &F) {
  Index const
  dimension = F.get_dimension();

  Tensor<T, N> const
  b = dot_t(F, F);

  Tensor<T, N>
  V(dimension), D(dimension);

  std::tie(V, D) = eig_sym(b);

  Tensor<T, N>
  DQ(dimension, Filler::ZEROS), DI(dimension, Filler::ZEROS), DL(dimension, Filler::ZEROS);

  for (Index i = 0; i < dimension; ++i) {
    if (D(i,i) < T(0)) {
      MT_ERROR_EXIT("Non-SPD input: negative eigenvalue.");
    }
    DQ(i,i) = std::sqrt(D(i,i));
    DI(i,i) = 1.0 / DQ(i,i);
    DL(i,i) = std::log(DQ(i,i));
  }

  Tensor<T, N> const
  R = dot(V, DI) * t_dot(V, F);

  Tensor<T, N> const
  X = V * dot_t(DQ, V);

  Tensor<T, N> const
  x = V * dot_t(DL, V);

  return std::make_tuple(X, R, x);
}

//
// R^N left polar decomposition with matrix logarithm for V
// \param F tensor (often a deformation-gradient-like tensor)
// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD(N), and log V
//
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
polar_left_logV_lame(Tensor<T, N> const &F) {
  Index const
  dimension = F.get_dimension();

  // set up return tensors
  Tensor<T, N> R(dimension), V(dimension), v(dimension), Vinv(dimension);

  // compute spd tensor
  Tensor<T, N> b = F*transpose(F);

  // get eigenvalues/eigenvectors
  Tensor<T, N> eVal(dimension);
  Tensor<T, N> eVec(dimension);
  std::tie(eVec, eVal) = eig_spd_cos(b);

  // compute sqrt() and inv(sqrt()) of eigenvalues
  for (Index i = 0; i < 3; ++i) {
    if (eVal(i,i) < T(0)) {
      MT_ERROR_EXIT("Non-SPD input: negative eigenvalue.");
    }
  }
  Tensor<T, N> x = zero<T, N>(3);
  x(0,0) = std::sqrt(eVal(0,0));
  x(1,1) = std::sqrt(eVal(1,1));
  x(2,2) = std::sqrt(eVal(2,2));
  Tensor<T, N> xi = zero<T, N>(3);
  xi(0,0) = 1.0/x(0,0);
  xi(1,1) = 1.0/x(1,1);
  xi(2,2) = 1.0/x(2,2);
  Tensor<T, N> lnx = zero<T, N>(3);
  lnx(0,0) = std::log(x(0,0));
  lnx(1,1) = std::log(x(1,1));
  lnx(2,2) = std::log(x(2,2));
  // compute V, Vinv, log(V)=v, and R
  V    = eVec*x*transpose(eVec);
  Vinv = eVec*xi*transpose(eVec);
  v    = eVec*lnx*transpose(eVec);
  R    = Vinv*F;

  return std::make_tuple(V, R, v);
}

//
// Symmetric Schur algorithm for R^2.
// \param \f$ A = [f, g; g, h] \in S(2) \f$
// \return \f$ c, s \rightarrow [c, -s; s, c]\f diagonalizes A$
//
template <typename T>
std::pair<T, T> schur_sym(T const f, T const g, T const h) {
  T c = 1.0;
  T s = 0.0;

  if (g != 0.0) {
    T t = (h - f) / (2.0 * g);

    if (t >= 0.0) {
      t = 1.0 / (std::sqrt(1.0 + t * t) + t);
    } else {
      t = -1.0 / (std::sqrt(1.0 + t * t) - t);
    }
    c = 1.0 / std::sqrt(1.0 + t * t);
    s = t * c;
  }

  return std::make_pair(c, s);
}

//
// Givens rotation. [c, -s; s, c] [a; b] = [r; 0]
// \param a, b
// \return c, s
//
template <typename T> std::pair<T, T> givens(T const &a, T const &b) {
  T c = 1.0;
  T s = 0.0;

  if (b != 0.0) {
    if (std::abs(b) > std::abs(a)) {
      T const t = - a / b;
      s = 1.0 / std::sqrt(1.0 + t * t);
      c = t * s;
    } else {
      T const t = - b / a;
      c = 1.0 / std::sqrt(1.0 + t * t);
      s = t * c;
    }
  }

  return std::make_pair(c, s);
}

namespace {

//
// R^N eigenvalue decomposition for symmetric 2nd-order tensor
// \param A tensor
// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
// See algorithm 8.4.2 in Matrix Computations, Golub & Van Loan 1996
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_sym_NxN(Tensor<T, N> const &A) {
  Tensor<T, N>
  D = sym(A);

  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  V = identity<T, N>(dimension);

  T
  off = norm_off_diagonal(D);

  T
  tol = machine_epsilon<T>() * norm(A);

  // Estimate based on random generation and linear regression.
  // Golub & Van Loan p 429 expect ~ dimension * log(dimension)
  Index const
  max_iter = 5 * dimension * dimension / 2;

  Index
  num_iter = 0;

  while (off > tol && num_iter < max_iter) {

    // Find largest off-diagonal entry
    Index
    p = 0;

    Index
    q = 0;

    std::tie(p, q) = arg_max_off_diagonal(D);
    if (p > q) {
      std::swap(p,q);
    }

    // Obtain Givens rotations by using 2x2 symmetric Schur algorithm
    T const &
    f = D(p,p);

    T const &
    g = D(p,q);

    T const &
    h = D(q,q);

    T
    c, s;

    std::tie(c, s) = schur_sym(f, g, h);

    // Apply Givens rotation to matrices
    // that are converging to eigenvalues and eigenvectors
    givens_left(c, s, p, q, D);
    givens_right(c, s, p, q, D);

    givens_right(c, s, p, q, V);

    off = norm_off_diagonal(D);
    num_iter++;
  }

  Vector<T, N> d(dimension);
  Tensor<T, N> P(dimension);

  std::tie(d, P) = sort_permutation(diag(D));
  D = diag(d);
  V = V * P;

  return std::make_pair(V, D);
}

//
// R^2 eigenvalue decomposition for symmetric 2nd-order tensor
// \param A tensor
// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_sym_2x2(Tensor<T, N> const &A) {
  assert(A.get_dimension() == 2);

  T const f = A(0,0);
  T const g = 0.5 * (A(0,1) + A(1,0));
  T const h = A(1,1);

  //
  // Eigenvalues, based on LAPACK's dlae2
  //
  T const sum = f + h;
  T const dif = std::abs(f - h);
  T const g2 = std::abs(g + g);

  T fhmax = f;
  T fhmin = h;

  const bool swap_diag = std::abs(h) > std::abs(f);

  if (swap_diag == true) {
    std::swap(fhmax, fhmin);
  }

  T r = 0.0;
  if (dif > g2) {
    T const t = g2 / dif;
    r = dif * std::sqrt(1.0 + t * t);
  } else if (dif < g2) {
    T const t = dif / g2;
    r = g2 * std::sqrt(1.0 + t * t);
  } else {
    // dif == g2, including zero
        r = g2 * std::sqrt(2.0);
  }

  T s0 = 0.0;
  T s1 = 0.0;

  if (sum != 0.0) {
    s0 = 0.5 * (sum + copysign(r, sum));
    // Order of execution important.
    // To get fully accurate smaller eigenvalue,
    // next line needs to be executed in higher precision.
    s1 = (fhmax / s0) * fhmin - (g / s0) * g;
  } else {
    // s0 == s1, including zero
    s0 = 0.5 * r;
    s1 = -0.5 * r;
  }

  //
  // Eigenvectors. schur_sym returns the Jacobi rotation J = [c, s; -s, c]
  // whose columns are the eigenvectors and which diagonalizes A as
  // A = J * diag(l0, l1) * transpose(J), with the eigenvalue attached to
  // column 0 of J being l0 = f - tan(theta) * g and to column 1 being
  // l1 = h + tan(theta) * g (tan(theta) = s / c).
  //
  // The previous implementation used transpose(J) as V and paired the
  // magnitude-ordered eigenvalues (s0, s1) to the columns through swap_diag, so
  // V * D * transpose(V) did not reconstruct A for general off-diagonal cases
  // (Trilinos issue #15389). Build V = J directly, and attach the accurate
  // (s0, s1) eigenvalues to their matching columns.
  //
  T
  c, s;

  std::tie(c, s) = schur_sym(f, g, h);

  Tensor<T, N>
  V(c, s, -s, c);

  T const l0 = f - (s / c) * g;

  T d0, d1;
  if (std::abs(l0 - s0) <= std::abs(l0 - s1)) {
    d0 = s0;
    d1 = s1;
  } else {
    d0 = s1;
    d1 = s0;
  }

  //
  // Return eigenvalues in descending order, permuting the eigenvectors to
  // match, for consistency with the eig_sym_NxN path. For 2x2 this is a single
  // compare-and-swap; done inline to avoid the heap allocation and matrix
  // multiply of the general sort_permutation, as eig_sym is called per
  // integration point in 2D mechanics.
  //
  if (d0 < d1) {
    std::swap(d0, d1);
    std::swap(V(0, 0), V(0, 1));
    std::swap(V(1, 0), V(1, 1));
  }

  Tensor<T, N>
  D(d0, 0.0, 0.0, d1);

  return std::make_pair(V, D);
}

} // anonymous namespace

//
// R^N eigenvalue decomposition for symmetric 2nd-order tensor
// \param A tensor
// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_sym(Tensor<T, N> const &A) {
  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  V(dimension), D(dimension);

  switch (dimension) {

    default:
      std::tie(V, D) = eig_sym_NxN(A);
      break;

    case 2:
      std::tie(V, D) = eig_sym_2x2(A);
      break;

  }

  return std::make_pair(V, D);
}

//
// R^N eigenvalue decomposition for SPD 2nd-order tensor
// \param A tensor
// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_spd(Tensor<T, N> const &A) {
  return eig_sym(A);
}

//
// R^3 eigenvalue decomposition for SPD 2nd-order tensor
// \param A tensor
// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_spd_cos(Tensor<T, N> const &A) {
  Index const
  dimension = A.get_dimension();

  assert(dimension == 3);

  // This algorithm comes from the journal article
  // Scherzinger and Dohrmann, CMAME 197 (2008) 4007-4015

  // this algorithm will return the eigenvalues in D
  // and the eigenvectors in V
  Tensor<T, N>
  D = zero<T, N>(dimension);

  Tensor<T, N>
  V = zero<T, N>(dimension);

  // not sure if this is necessary...
  T
  pi = std::acos(-1);

  // convenience operators
  Tensor<T, N> const
  I = identity<T, N>(dimension);

  int
  ii[3][2] = { { 1, 2 }, { 2, 0 }, { 0, 1 } };

  Tensor<T, N>
  rm = zero<T, N>(dimension);

  // scale the matrix to reduce the characteristic equation
  T
  trA = (1.0/3.0) * I1(A);

  Tensor<T, N>
  Ap(A - trA*I);

  // compute other invariants
  T
  J2 = I2(Ap);

  T
  J3 = det(Ap);

  // deal with volumetric tensors
  if (-J2 <= 1.e-30)
  {
    D(0,0) = trA;
    D(1,1) = trA;
    D(2,2) = trA;

    V(0,0) = 1.0;
    V(1,0) = 0.0;
    V(2,0) = 0.0;

    V(0,1) = 0.0;
    V(1,1) = 1.0;
    V(2,1) = 0.0;

    V(0,2) = 0.0;
    V(1,2) = 0.0;
    V(2,2) = 1.0;
  }
  else
  {
    // first things first, find the most dominant e-value
    // Need to solve cos(3 theta)=rhs for theta
    T
    t1 = 3.0 / -J2;

    T
    rhs = (J3 / 2.0) * T(std::sqrt(t1 * t1 * t1));

    T
    theta = pi / 2.0 * (1.0 - (rhs < 0 ? -1.0 : 1.0));

    if (std::abs(rhs) <= 1.0) theta = std::acos(rhs);

    T
    thetad3 = theta / 3.0;

    if (thetad3 > pi / 6.0) thetad3 += 2.0 * pi / 3.0;

    // most dominant e-value
    D(2,2) = 2.0 * std::cos(thetad3) * std::sqrt(-J2 / 3.0);

    // now reduce the system
    Tensor<T, N>
    R = Ap - D(2,2) * I;

    // QR factorization with column pivoting
    Vector<T, N> a(dimension);
    a(0) = R(0,0)*R(0,0) + R(1,0)*R(1,0) + R(2,0)*R(2,0);
    a(1) = R(0,1)*R(0,1) + R(1,1)*R(1,1) + R(2,1)*R(2,1);
    a(2) = R(0,2)*R(0,2) + R(1,2)*R(1,2) + R(2,2)*R(2,2);

    // find the most dominant column
    int k = 0;
    T max = a(0);
    if (a(1) > max)
    {
      k = 1;
      max = a(1);
    }
    if (a(2) > max)
    {
      k = 2;
    }

    // normalize the most dominant column to get s1
    a(k) = std::sqrt(a(k));
    for (int i(0); i < dimension; ++i)
      R(i,k) /= a(k);

    // dot products of dominant column with other two columns
    T d0 = 0.0;
    T d1 = 0.0;
    for (int i(0); i < dimension; ++i)
    {
      d0 += R(i,k) * R(i,ii[k][0]);
      d1 += R(i,k) * R(i,ii[k][1]);
    }

    // projection
    for (int i(0); i < dimension; ++i)
    {
      R(i,ii[k][0]) -= d0 * R(i,k);
      R(i,ii[k][1]) -= d1 * R(i,k);
    }

    // now finding next most dominant column
    a.clear();
    for (int i(0); i < dimension; ++i)
    {
      a(0) += R(i,ii[k][0]) * R(i,ii[k][0]);
      a(1) += R(i,ii[k][1]) * R(i,ii[k][1]);
    }

    int p = 0;
    if (std::abs(a(1)) > std::abs(a(0))) p = 1;

    // normalize next most dominant column to get s2
    a(p) = std::sqrt(a(p));
    int k2 = ii[k][p];

    for (int i(0); i < dimension; ++i)
      R(i,k2) /= a(p);

    // set first eigenvector as cross product of s1 and s2
    V(0,2) = R(1,k) * R(2,k2) - R(2,k) * R(1,k2);
    V(1,2) = R(2,k) * R(0,k2) - R(0,k) * R(2,k2);
    V(2,2) = R(0,k) * R(1,k2) - R(1,k) * R(0,k2);

    // normalize
    T
    mag = std::sqrt(V(0,2) * V(0,2) + V(1,2) * V(1,2) + V(2,2) * V(2,2));

    V(0,2) /= mag;
    V(1,2) /= mag;
    V(2,2) /= mag;

    // now for the other two eigenvalues, extract vectors
    Vector<T, N>
    rk(R(0,k), R(1,k), R(2,k));

    Vector<T, N>
    rk2(R(0,k2), R(1,k2), R(2,k2));

    // compute projections
    Vector<T, N>
    ak = Ap * rk;

    Vector<T, N>
    ak2 = Ap * rk2;

    // set up reduced remainder matrix
    rm(0,0) = dot(rk,ak);
    rm(0,1) = dot(rk,ak2);
    rm(1,1) = dot(rk2,ak2);

    // compute eigenvalues 2 and 3
    T
    b = 0.5 * (rm(0,0) - rm(1,1));

    T
    fac = (b < 0 ? -1.0 : 1.0);

    T
    arg = b * b + rm(0,1) * rm(0,1);

    if (arg == 0)
      D(0,0) = rm(1,1) + b;
    else
      D(0,0) = rm(1,1) + b - fac * std::sqrt(b * b + rm(0,1) * rm(0,1));

    D(1,1) = rm(0,0) + rm(1,1) - D(0,0);

    // update reduced remainder matrix
    rm(0,0) -= D(0,0);
    rm(1,0) = rm(0,1);
    rm(1,1) -= D(0,0);

    // again, find most dominant column
    a.clear();
    a(0) = rm(0,0) * rm(0,0) + rm(0,1) * rm(0,1);
    a(1) = rm(0,1) * rm(0,1) + rm(1,1) * rm(1,1);

    int k3 = 0;
    if (a(1) > a(0)) k3 = 1;
    if (a(k3) == 0.0)
    {
      rm(0,k3) = 1.0;
      rm(1,k3) = 0.0;
    }

    // set 2nd eigenvector via cross product
    V(0,0) = rm(0,k3) * rk2(0) - rm(1,k3) * rk(0);
    V(1,0) = rm(0,k3) * rk2(1) - rm(1,k3) * rk(1);
    V(2,0) = rm(0,k3) * rk2(2) - rm(1,k3) * rk(2);

    // normalize
    mag = std::sqrt(V(0,0) * V(0,0) + V(1,0) * V(1,0) + V(2,0) * V(2,0));
    V(0,0) /= mag;
    V(1,0) /= mag;
    V(2,0) /= mag;

    // set last eigenvector as cross product of other two
    V(0,1) = V(1,0) * V(2,2) - V(2,0) * V(1,2);
    V(1,1) = V(2,0) * V(0,2) - V(0,0) * V(2,2);
    V(2,1) = V(0,0) * V(1,2) - V(1,0) * V(0,2);

    // normalize
    mag = std::sqrt(V(0,1) * V(0,1) + V(1,1) * V(1,1) + V(2,1) * V(2,1));
    V(0,1) /= mag;
    V(1,1) /= mag;
    V(2,1) /= mag;

    // add back in the offset
    for (int i(0); i < dimension; ++i)
      D(i,i) += trA;
  }

  return std::make_pair(V, D);
}

//
// Cholesky decomposition, rank-1 update algorithm
// (Matrix Computations 3rd ed., Golub & Van Loan, p145)
// \param A assumed symmetric tensor
// \return G Cholesky factor A = GG^T
// \return completed (bool) algorithm ran to completion
//
template <typename T, Index N>
std::pair<Tensor<T, N>, bool> cholesky(Tensor<T, N> const &A) {
  Tensor<T, N>
  G = sym(A);

  Index const
  dimension = A.get_dimension();

  for (Index k = 0; k < dimension; ++k) {

    // Zeros above the diagonal
    for (Index j = k + 1; j < dimension; ++j) {
      G(k,j) = 0.0;
    }

    T
    s = G(k,k);

    if (s <= 0.0) {
      return std::make_pair(G, false);
    }

    s = std::sqrt(s);

    for (Index j = k + 1; j < dimension; ++j) {
      G(j,k) /= s;
    }

    G(k,k) = s;

    for (Index j = k + 1; j < dimension; ++j) {
      for (Index i = j; i < dimension; ++i) {
        G(i,j) -= G(i,k) * G(j,k);
      }
    }

  }

  return std::make_pair(G, true);
}

} // namespace minitensor

#endif // MiniTensor_Factorizations_h
