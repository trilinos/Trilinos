// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_ArithTraits.hpp"
#if !defined(MiniTensor_LinearAlgebra_t_h)
#define MiniTensor_LinearAlgebra_t_h

namespace minitensor {

//
// Inverse defaults to fast inverse for 2 and 3 dimensions, otherwise
// use full piviting version
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
inverse(Tensor<T, N> const & A)
{
  return inverse_fast23(A);
}

//
// R^N 2nd-order tensor inverse
// Gauss-Jordan elimination. Warning: full pivoting for small tensors.
// Use Teuchos LAPACK interface for more efficient and robust techniques.
// \param A nonsingular tensor
// \return \f$ A^{-1} \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
inverse_full_pivot(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  B = identity<T, N>(dimension);

  return solve_full_pivot(A, B);
}

//
// R^N 2nd-order tensor inverse
// Fast analytic expressions for 2 and 3 dimensions
// \param A nonsingular tensor
// \return \f$ A^{-1} \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
inverse_fast23(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  switch (dimension) {

  case 3:
    {
      T const determinant = det(A);
      assert(determinant != 0.0);
      return Tensor<T, N>(
        -A(1,2)*A(2,1) + A(1,1)*A(2,2),
         A(0,2)*A(2,1) - A(0,1)*A(2,2),
        -A(0,2)*A(1,1) + A(0,1)*A(1,2),
         A(1,2)*A(2,0) - A(1,0)*A(2,2),
        -A(0,2)*A(2,0) + A(0,0)*A(2,2),
         A(0,2)*A(1,0) - A(0,0)*A(1,2),
        -A(1,1)*A(2,0) + A(1,0)*A(2,1),
         A(0,1)*A(2,0) - A(0,0)*A(2,1),
        -A(0,1)*A(1,0) + A(0,0)*A(1,1)
        ) / determinant;
    }

  case 2:
    {
      T const determinant = det(A);
      assert(determinant != 0.0);
      return Tensor<T, N>(A(1,1), -A(0,1), -A(1,0), A(0,0)) / determinant;
    }

  case 1:
    return Tensor<T, N>(1, Filler::ONES) / A(0,0);

  default:
    break;
  }

  return inverse_full_pivot(A);
}

//
//
//
template<typename T, Index N, typename RHS>
KOKKOS_INLINE_FUNCTION
RHS
solve_full_pivot(Tensor<T, N> const & A, RHS const & b)
{
  Index const
  dimension{A.get_dimension()};

  Index const
  maximum_dimension{INDEX_SIZE};

  if (dimension > maximum_dimension) {
    MT_ERROR_EXIT("Max dim (%d) exceeded: %d.", dimension, maximum_dimension);
  }

  RHS
  B{b};

  Index const
  num_rhs{B.get_num_cols()};

  switch (dimension) {

  case 1:
    for (Index i{0}; i < num_rhs; ++i) {
      B(0, i) = b(0, i) / A(0, 0);
    }
    return B;

  default:
    break;
  }

  Tensor<T, N>
  S{A};

  // Set 1 ... dimension bits to one.
  Index
  intact_rows{static_cast<Index>((1UL << dimension) - 1)};

  Index
  intact_cols{static_cast<Index>((1UL << dimension) - 1)};

  // Gauss-Jordan elimination with full pivoting
  for (Index k{0}; k < dimension; ++k) {

    // Determine full pivot
    T
    pivot{0.0};

    Index
    pivot_row{dimension};

    Index
    pivot_col{dimension};

    for (Index row{0}; row < dimension; ++row) {

      if (!(intact_rows & (1 << row))) continue;

      for (Index col{0}; col < dimension; ++col) {

        if (!(intact_cols & (1 << col))) continue;

        T
        s{std::abs(S(row, col))};
        if (s > pivot) {

          pivot_row = row;
          pivot_col = col;
          pivot = s;

        }

      }

    }

    assert(pivot_row < dimension);
    assert(pivot_col < dimension);

    // Gauss-Jordan elimination
    T const
    t{S(pivot_row, pivot_col)};

    assert(t != 0.0);

    for (Index j{0}; j < dimension; ++j) {
      S(pivot_row, j) /= t;
    }
    for (Index j{0}; j < num_rhs; ++j) {
      B(pivot_row, j) /= t;
    }

    for (Index i{0}; i < dimension; ++i) {
      if (i == pivot_row) continue;

      T const
      c{S(i, pivot_col)};

      for (Index j = 0; j < dimension; ++j) {
        S(i, j) -= c * S(pivot_row, j);
      }
      for (Index j = 0; j < num_rhs; ++j) {
        B(i, j) -= c * B(pivot_row, j);
      }
    }

    // Eliminate current row and col from intact rows and cols
    intact_rows &= ~(1 << pivot_row);
    intact_cols &= ~(1 << pivot_col);

  }

  RHS const
  X = t_dot(S, B);

  return X;
}

//
// R^N Subtensor
// \param A tensor
// \param i index
// \param j index
// \return Subtensor with i-row and j-col deleted.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
subtensor(Tensor<T, N> const & A, Index const i, Index const j)
{
  Index const
  dimension = A.get_dimension();

  assert(i < dimension);
  assert(j < dimension);

  Tensor<T, N>
  B(dimension - 1);

  Index p = 0;
  for (Index m = 0; m < dimension; ++m) {
    if (m == i) continue;
    Index q = 0;
    for (Index n = 0; n < dimension; ++n) {
      if (n == j) continue;
      B(p, q) = A(m, n);
      ++q;
    }
    ++p;
  }

  return B;
}

//
// Exponential map
//
template <typename T, Index N> Tensor<T, N> exp(Tensor<T, N> const &A) {
  return exp_pade(A);
}

//
// R^N exponential map by Taylor series, radius of convergence is infinity
// \param A tensor
// \return \f$ \exp A \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
exp_taylor(Tensor<T, N> const & A)
{
  Index const
  max_iter = 128;

  T const
  tol = machine_epsilon<T>();

  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  term = identity<T, N>(dimension);

  // Relative error taken wrt to the first term, which is I and norm = 1
  T
  relative_error = 1.0;

  Tensor<T, N>
  B = term;

  Index
  k = 0;

  while (relative_error > tol && k < max_iter) {
    term = static_cast<T>(1.0 / (k + 1.0)) * term * A;
    B = B + term;
    relative_error = norm_1(term);
    ++k;
  }

  return B;
}

namespace {

//
// Scaling parameter theta for scaling and squaring exponential.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
scaling_squaring_theta(Index const order)
{
  assert(order > 0 && order < 22);

  T const theta[] =
  {
      0.0e-0, 3.7e-8, 5.3e-4, 1.5e-2, 8.5e-2, 2.5e-1, 5.4e-1, 9.5e-1,
      1.5e-0, 2.1e-0, 2.8e-0, 3.6e-0, 4.5e-0, 5.4e-0, 6.3e-0, 7.3e-0,
      8.4e-0, 9,4e-0, 1.1e+1, 1.2e+1, 1.3e+1, 1.4e+1
  };

  return theta[order];
}

//
// Polynomial coefficients for Padé approximants.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
polynomial_coefficient(Index const order, Index const index)
{
  assert(index <= order);

  T
  c = 0.0;

  switch (order) {

    default:
      MT_ERROR_EXIT("Wrong order in Pade' polynomial coefficient: ");
      break;

    case 3:
    {
      T const
      b[] = {120.0, 60.0, 12.0, 1.0};

      c = b[index];
    }
    break;

    case 5:
    {
      T const
      b[] = {30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0};

      c = b[index];
    }
    break;

    case 7:
    {
      T const
      b[] = {17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0,
          56.0, 1.0};

      c = b[index];
    }
    break;

    case 9:
    {
      T const
      b[] = {17643225600.0, 8821612800.0, 2075673600.0, 302702400.0,
          30270240.0, 2162160.0, 110880.0, 3960.0, 90.0, 1.0};

      c = b[index];
    }
    break;

    case 13:
    {
      T const
      b[] = {64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
          1187353796428800.0, 129060195264000.0, 10559470521600.0,
          670442572800.0, 33522128640.0, 1323241920.0, 40840800.0,
          960960.0, 16380.0, 182.0, 1.0};

      c = b[index];
    }
    break;

  }

  return c;
}

//
// Padé approximant polynomial odd and even terms.
//
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>>
pade_polynomial_terms(Tensor<T, N> const &A, Index const order) {
  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  B = identity<T, N>(dimension);

  Tensor<T, N>
  U = polynomial_coefficient<Real>(order, 1) * B;

  Tensor<T, N>
  V = polynomial_coefficient<Real>(order, 0) * B;

  Tensor<T, N> const
  A2 = A * A;

  for (Index i = 3; i <= order; i += 2) {

    B = B * A2;

    Tensor<T, N> const
    O = polynomial_coefficient<Real>(order, i) * B;

    Tensor<T, N> const
    E = polynomial_coefficient<Real>(order, i - 1) * B;

    U += O;

    V += E;

  }

  U = A * U;

  return std::make_pair(U, V);
}

//
// Compute a non-negative integer power of a tensor by binary manipulation.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
binary_powering(Tensor<T, N> const & A, Index const exponent)
{
  if (exponent == 0) return eye<T, N>(A.get_dimension());

  Index const
  rightmost_bit = 1;

  Index const
  number_digits = INDEX_SIZE;

  Index const
  leftmost_bit = rightmost_bit << (number_digits - 1);

  Index
  t = 0;

  for (Index j = 0; j < number_digits; ++j) {

    if (((exponent << j) & leftmost_bit) != 0) {

      t = number_digits - j - 1;
      break;

    }

  }

  Tensor<T, N>
  P = A;

  Index
  i = 0;

  Index
  m = exponent;

  while ((m & rightmost_bit) == 0) {
    P = P * P;
    ++i;
    m = m >> 1;
  }

  Tensor<T, N>
  X = P;

  for (Index j = i + 1; j <= t; ++j) {
    P = P * P;

    if (((exponent >> j) & rightmost_bit) != 0) {
      X = X * P;
    }
  }

  return X;
}

} // anonymous namespace

//
// Exponential map by squaring and scaling and Padé approximants.
// See algorithm 10.20 in Functions of Matrices, N.J. Higham, SIAM, 2008.
// \param A tensor
// \return \f$ \exp A \f$
//
template <typename T, Index N> Tensor<T, N> exp_pade(Tensor<T, N> const &A) {
  Index const
  dimension = A.get_dimension();

  Index const
  orders[] = {3, 5, 7, 9, 13};

  Index const
  number_orders = 5;

  Index const
  highest_order = orders[number_orders - 1];

  Tensor<T, N>
  B;

  Real const
  norm = Sacado::ScalarValue<T>::eval((norm_1(A)));

  for (Index i = 0; i < number_orders; ++i) {

    Index const
    order = orders[i];

    Real const
    theta = scaling_squaring_theta<Real>(order);

    if (order < highest_order && norm < theta) {

      Tensor<T, N>
      U;

      Tensor<T, N>
      V;

      std::tie(U, V) = pade_polynomial_terms(A, order);

      B = inverse(V - U) * (U + V);

      break;

    } else if (order == highest_order) {

      Real const
      theta_highest = scaling_squaring_theta<Real>(order);

      int const
      signed_power = static_cast<int>(std::ceil(std::log2(norm / theta_highest)));
      Index const
      power_two = signed_power > 0 ? static_cast<Index>(signed_power) : 0;

      Real
      scale = 1.0;

      for (Index j = 0; j < power_two; ++j) {
        scale /= 2.0;
      }

      Tensor<T, N> const
      I = identity<T, N>(dimension);

      Tensor<T, N> const
      A1 = scale * A;

      Tensor<T, N> const
      A2 = A1 * A1;

      Tensor<T, N> const
      A4 = A2 * A2;

      Tensor<T, N> const
      A6 = A2 * A4;

      Real const b0  = polynomial_coefficient<Real>(order, 0);
      Real const b1  = polynomial_coefficient<Real>(order, 1);
      Real const b2  = polynomial_coefficient<Real>(order, 2);
      Real const b3  = polynomial_coefficient<Real>(order, 3);
      Real const b4  = polynomial_coefficient<Real>(order, 4);
      Real const b5  = polynomial_coefficient<Real>(order, 5);
      Real const b6  = polynomial_coefficient<Real>(order, 6);
      Real const b7  = polynomial_coefficient<Real>(order, 7);
      Real const b8  = polynomial_coefficient<Real>(order, 8);
      Real const b9  = polynomial_coefficient<Real>(order, 9);
      Real const b10 = polynomial_coefficient<Real>(order, 10);
      Real const b11 = polynomial_coefficient<Real>(order, 11);
      Real const b12 = polynomial_coefficient<Real>(order, 12);
      Real const b13 = polynomial_coefficient<Real>(order, 13);

      Tensor<T, N> const
      U = A1 * (
          (A6 * (b13 * A6 + b11 * A4 + b9 * A2) +
              b7 * A6 + b5 * A4 + b3 * A2 + b1 * I));

      Tensor<T, N> const
      V = A6 * (b12 * A6 + b10 * A4 + b8 * A2) +
      b6 * A6 + b4 * A4 + b2 * A2 + b0 * I;

      Tensor<T, N> const
      R = inverse(V - U) * (U + V);

      Index const
      exponent = (1U << power_two);

      B = binary_powering(R, exponent);

    }

  }

  return B;
}

//
// Logarithmic map by Taylor series.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_taylor(Tensor<T, N> const & A)
{
  Index const
  max_iter = 128;

  T const
  tol = machine_epsilon<T>();

  T const
  norm_tensor = norm_1(A);

  Index const
  dimension = A.get_dimension();

  Tensor<T, N> const
  A_minus_I = A - identity<T, N>(dimension);

  Tensor<T, N>
  term = A_minus_I;

  T
  norm_term = norm_1(term);

  T
  relative_error = norm_term / norm_tensor;

  Tensor<T, N>
  B = term;

  Index
  k = 1;

  while (relative_error > tol && k <= max_iter) {
    term = static_cast<T>(- (k / (k + 1.0))) * term * A_minus_I;
    B = B + term;
    norm_term = norm_1(term);
    relative_error = norm_term / norm_tensor;
    ++k;
  }

  return B;
}

//
// Logarithmic map.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log(Tensor<T, N> const & A)
{
  return log_gregory(A);
}

//
// Logarithmic map by Gregory series.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_gregory(Tensor<T, N> const & A)
{
  Index const
  max_iter = 8192;

  T const
  tol = machine_epsilon<T>();

  T const
  norm_tensor = norm_1(A);

  Index const
  dimension = A.get_dimension();

  Tensor<T, N> const
  I_minus_A = identity<T, N>(dimension) - A;

  Tensor<T, N> const
  I_plus_A = identity<T, N>(dimension) + A;

  Tensor<T, N>
  term = I_minus_A * inverse(I_plus_A);

  T
  norm_term = norm_1(term);

  T
  relative_error = norm_term / norm_tensor;

  Tensor<T, N> const
  C = term * term;

  Tensor<T, N>
  B = term;

  Index
  k = 1;

  while (relative_error > tol && k <= max_iter + 1) {
    term = static_cast<T>((2 * k - 1.0) / (2 * k + 1.0)) * term * C;
    B = B + term;
    norm_term = norm_1(term);
    relative_error = norm_term / norm_tensor;
    ++k;
  }

  B = - 2.0 * B;

  return B;
}

//
// Logarithmic map for symmetric tensor.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_sym(Tensor<T, N> const & A)
{
  return log_eig_sym(A);
}

//
// Logarithmic map for symmetric tensor using eigenvalue decomposition.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_eig_sym(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  V(dimension);

  Tensor<T, N>
  D(dimension);

  std::tie(V, D) = eig_sym(A);

  for (Index i = 0; i < dimension; ++i) {
    D(i, i) = std::log(D(i, i));
  }

  Tensor<T, N> const
  B = dot_t(dot(V, D), V);

  return B;
}

//
// R^N logarithmic map of a rotation.
// \param R with \f$ R \in SO(N) \f$
// \return \f$ r = \log R \f$ with \f$ r \in so(N) \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_rotation(Tensor<T, N> const & R)
{
  Index const
  dimension = R.get_dimension();

  //firewalls, make sure R \in SO(N)
  assert(norm(dot_t(R,R) - eye<T, N>(dimension)) <
         std::max(1.0e-12 * norm(R), 1.0e-12));
  assert(std::abs(det(R) - 1.0) <
         std::max(1.0e-12 * norm(R), 1.0e-12));
  // acos requires input between -1 and +1
  T
  cosine = 0.5 * (trace(R) - 1.0);

  if (cosine < -1.0) {
    cosine = -1.0;
  } else if(cosine > 1.0) {
    cosine = 1.0;
  }
  T
  theta = std::acos(cosine);

  Tensor<T, N>
  r(dimension);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Logarithm of SO(N) N != 2,3 not implemented.");
      break;

    case 3:
      if (theta == 0.0) {

        r = zero<T, N>(3);

      } else if (std::abs(cosine + 1.0) < 10.0 * machine_epsilon<T>())  {

        r = log_rotation_pi(R);

      } else {

        r = theta / std::sin(theta) * skew(R);

      }
      break;

    case 2:
      r(0,0) = 0.0;
      r(0,1) = -theta;
      r(1,0) = theta;
      r(1,1) = 0.0;
      break;

    case 1:
      r(0,0) = 0.0;
      break;

  }

  return r;
}

// R^N Logarithmic map of a 180-degree rotation.
// \param R with \f$ R \in SO(N) \f$
// \return \f$ r = \log R \f$ with \f$ r \in so(N) \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_rotation_pi(Tensor<T, N> const & R)
{
  Index const
  dimension = R.get_dimension();

  // set firewall to make sure the rotation is indeed 180 degrees
  assert(std::abs(trace(R) + 1.0) < 10.0 * machine_epsilon<T>());

  Tensor<T, N>
  r(dimension);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Logarithm of SO(N) N != 2,3 not implemented.");
      break;

    case 3:
    {
      Vector<T, N>
      normal(3);

      Tensor<T, N> const
      B = R - identity<T, N>(3);

      Vector<T, N> const
      u = row(B, 0);

      Vector<T, N> const
      v = row(B, 1);

      normal = cross(u, v);

      if (norm(normal) < machine_epsilon<T>()) {

        Vector<T, N> const
        w = row(B, 2);

        normal = cross(u, w);

        if (norm(normal) < machine_epsilon<T>()) {
          MT_ERROR_EXIT("Cannot determine rotation vector of rotation.");
        }

      }

      normal = unit(normal);

      r.fill(Filler::ZEROS);
      r(0,1) = -normal(2);
      r(0,2) =  normal(1);
      r(1,0) =  normal(2);
      r(1,2) = -normal(0);
      r(2,0) = -normal(1);
      r(2,1) =  normal(0);

      T const
      pi = std::acos(-1.0);

      r = pi * r;
    }
    break;

    case 2:
    {
      T theta = std::acos(-1.0);
      if (R(0,0) > 0.0) {
        theta = -theta;
      }

      r(0,0) = 0.0;
      r(0,1) = -theta;
      r(1,0) = theta;
      r(1,1) = 0.0;
    }
    break;

  }

  return r;
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

//
// Apply rank-one update on the left in place
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
rank_one_left(T const & beta, Vector<T, N> const & v, Tensor<T, N> & A)
{
  A -= beta * dyad(v, dot(v, A));
  return;
}

//
// Apply rank-one update on the right in place
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
rank_one_right(T const & beta, Vector<T, N> const & v, Tensor<T, N> & A)
{
  A -= beta * dyad(dot(A, v), v);
  return;
}

//
// R^N exponential map of a skew-symmetric tensor.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
exp_skew_symmetric(Tensor<T, N> const & r)
{
  // Check whether skew-symmetry holds
  assert(norm(sym(r)) < std::max(1.0e-12 * norm(r), 1.0e-12));

  Index const
  dimension = r.get_dimension();

  Tensor<T, N>
  R = identity<T, N>(dimension);

  T
  theta = 0.0;

  switch (dimension) {

    default:
      R = exp(r);
      break;

    case 3:
      theta = sqrt(r(2,1)*r(2,1)+r(0,2)*r(0,2)+r(1,0)*r(1,0));

      //Check whether norm == 0. If so, return identity.
      if (theta >= machine_epsilon<T>()) {
        R += sin(theta) / theta * r +
            (1.0 - cos(theta)) / (theta * theta) * r * r;
      }
      break;

    case 2:
      theta = r(1,0);

      {
        T const
        c = std::cos(theta);

        T const
        s = std::sin(theta);

        R(0,0) = c;
        R(0,1) = -s;
        R(1,0) = s;
        R(1,1) = c;
      }

      break;

    case 1:
      R(0,0) = 1.0;
      break;

 }

  return R;
}

//
// R^N off-diagonal norm. Useful for SVD and other algorithms
// that rely on Jacobi-type procedures.
// \param A
// \return \f$ \sqrt(\sum_i \sum_{j, j\neq i} a_{ij}^2) \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_off_diagonal(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

    default:
      for (Index i = 0; i < dimension; ++i) {
        for (Index j = 0; j < dimension; ++j) {
          if (i != j) s += A(i,j)*A(i,j);
        }
      }
      break;

    case 3:
      s = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2) +
      A(1,0)*A(1,0) + A(2,0)*A(2,0) + A(2,1)*A(2,1);
      break;

    case 2:
      s = A(0,1)*A(0,1) + A(1,0)*A(1,0);
      break;

    case 1:
      s = 0.0;
      break;

  }
  return std::sqrt(s);
}

//
// R^N arg max abs. Useful for inverse and other algorithms
// that rely on Jacobi-type procedures.
// \param A
// \return \f$ (p,q) = arg max_{i,j} |a_{ij}| \f$
//
template <typename T, Index N>
std::pair<Index, Index> arg_max_abs(Tensor<T, N> const &A) {

  Index p = 0;
  Index q = 0;

  T
  s = std::abs(A(p,q));

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      if (std::abs(A(i,j)) > s) {
        p = i;
        q = j;
        s = std::abs(A(i,j));
      }
    }
  }

  return std::make_pair(p,q);
}

//
// R^N arg max off-diagonal. Useful for SVD and other algorithms
// that rely on Jacobi-type procedures.
// \param A
// \return \f$ (p,q) = arg max_{i \neq j} |a_{ij}| \f$
//
template <typename T, Index N>
std::pair<Index, Index> arg_max_off_diagonal(Tensor<T, N> const &A) {
  Index p = 0;
  Index q = 1;

  T s = std::abs(A(p,q));

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      if (i != j && std::abs(A(i,j)) > s) {
        p = i;
        q = j;
        s = std::abs(A(i,j));
      }
    }
  }

  return std::make_pair(p,q);
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
      Kokkos::ArithTraits<Index>::sqrt(dimension) * machine_epsilon<T>();

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
// R^N logarithmic map using BCH expansion (4 terms)
// \param x tensor
// \param y tensor
// \return Baker-Campbell-Hausdorff series up to 4 terms
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
bch(Tensor<T, N> const & x, Tensor<T, N> const & y)
{
  return
      // first order term
      x + y
      +
      // second order term
      0.5*(x*y - y*x)
      +
      // third order term
      1.0/12.0 *
      (x*x*y - 2.0*x*y*x + x*y*y + y*x*x - 2.0*y*x*y + y*y*x)
      +
      // fourth order term
      1.0/24.0 *
      (x*x*y*y - 2.0*x*y*x*y + 2.0*y*x*y*x - y*y*x*x);
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

  Tensor<T, N>
  D(s0, 0.0, 0.0, s1);

  //
  // Eigenvectors
  //
  T
  c, s;

  std::tie(c, s) = schur_sym(f, g, h);

  Tensor<T, N>
  V(c, -s, s, c);

  if (swap_diag == true) {
    // swap eigenvectors if eigenvalues were swapped
    std::swap(V(0, 0), V(0, 1));
    std::swap(V(1, 0), V(1, 1));
  }

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

// Auxiliary functions for precondioners.
namespace {

//
//
//
template <typename T, Index N, typename RHS>
std::pair<Tensor<T, N>, RHS> identity_precon(Tensor<T, N> const &A,
                                             RHS const &B) {
  return std::make_pair(A, B);
}

//
//
//
template <typename T, Index N, typename RHS>
std::pair<Tensor<T, N>, RHS> diagonal_precon(Tensor<T, N> const &A,
                                             RHS const &B) {
  Vector<T, N> const
  d = diag(A);

  Vector<T, N> const
  v = 1.0 / d;

  Tensor<T, N> const
  P = diag(v);

  return std::make_pair(P * A, P * B);
}

//
//
//
template <typename T, Index N, typename RHS>
std::pair<Tensor<T, N>, RHS> maxabsrow_precon(Tensor<T, N> const &A, RHS &B) {
  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  P(dimension, Filler::ZEROS);

  for (Index i{0}; i < dimension; ++i) {
    P(i, i) = 1.0 / norm_infinity(row(A, i));
  }

  return std::make_pair(P * A, P * B);
}

} // anonymous namespace

//
//
//
template <typename T, Index N, typename RHS>
std::pair<Tensor<T, N>, RHS> precon(PreconditionerType const pt,
                                    Tensor<T, N> const &A, RHS const &B) {
  switch (pt) {
  default:
    MT_ERROR_EXIT("Unknown preconditioner type.");
    break;

  case PreconditionerType::IDENTITY:
    break;

  case PreconditionerType::DIAGONAL:
    return diagonal_precon(A, B);

  case PreconditionerType::MAX_ABS_ROW:
    return maxabsrow_precon(A, B);
  }

  return std::make_pair(A, B);
}

//
// Solve linear system of equations.
// This is meant for the solution of small linear systems of equations
// typically found in constitutive updates.
// Right now the implementation is very inefficient (but accurate)
// as it just Gauss-Jordan elimination. It is intended to be used in
// conjunction with Kokkos to take advantage of thread parallelism.
//
template <typename T, Index N, typename RHS>
RHS solve(Tensor<T, N> const &A, RHS const &b, PreconditionerType const pt) {
  Tensor<T, N>
  PA;

  RHS
  Pb;

  std::tie(PA, Pb) = precon(pt, A, b);

  return solve_full_pivot(PA, Pb);
}

} // namespace minitensor

#endif // MiniTensor_LinearAlgebra_t_h
