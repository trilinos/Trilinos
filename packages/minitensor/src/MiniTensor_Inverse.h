// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Inverse_h)
#define MiniTensor_Inverse_h

// Tensor inversion and linear solves.
#include "MiniTensor_Norms.h"

namespace minitensor {

/// \addtogroup minitensor_inverse
/// @{

///
/// 2nd-order tensor inverse
/// \param A nonsingular tensor
/// \return \f$ A^{-1} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
inverse(Tensor<T, N> const & A);

///
/// 2nd-order tensor inverse using analitical expression for 2 and 3 dimensions
/// \param A nonsingular tensor
/// \return \f$ A^{-1} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
inverse_fast23(Tensor<T, N> const & A);

///
/// 2nd-order tensor inverse using full pivoting, very accurate
/// \param A nonsingular tensor
/// \return \f$ A^{-1} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
inverse_full_pivot(Tensor<T, N> const & A);

///
/// Swap row. Echange rows i and j in place
/// \param A tensor
/// \param i index
/// \param j index
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
swap_row(Tensor<T, N> & A, Index const i, Index const j);

///
/// Swap column. Echange columns i and j in place
/// \param A tensor
/// \param i index
/// \param j index
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
swap_col(Tensor<T, N> & A, Index const i, Index const j);

///
/// Apply rank-one update on the left in place
/// \f$ A = (I - beta v v^T) A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
rank_one_left(T const & beta, Vector<T, N> const & v, Tensor<T, N> & A);

///
/// Apply rank-one update on the right in place
/// \f$ A = A (I - beta v v^T) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
rank_one_right(T const & beta, Vector<T, N> const & v, Tensor<T, N> & A);

///
/// Preconditioner types
///
enum class PreconditionerType
{
  UNDEFINED = 0,
  IDENTITY = 1,
  DIAGONAL = 2,
  MAX_ABS_ROW = 3,
};

///
/// Compute a preconditioner for improving the conditioning of a
/// linear system.
///
template <typename T, Index N, typename RHS>
std::pair<Tensor<T, N>, RHS> precon(PreconditionerType const pt,
                                    Tensor<T, N> const &A, RHS const &B);

///
/// Solve linear system of equations.
/// This is meant for the solution of small linear systems of equations
/// typically found in constitutive updates.
/// Right now the implementation is very inefficient (but accurate)
/// as it just uses the inverse function. It is intended to be used in
/// conjunction with Kokkos to take advantage of thread parallelism.
/// \param A assumed non-singular tensor
/// \param b rhs of the system Ax=b
/// \param pt preconditioner type
/// \return x solution(s) to the system Ax=b
///
template <typename T, Index N, typename RHS>
RHS solve(Tensor<T, N> const &A, RHS const &b,
          PreconditionerType const pt = PreconditionerType::IDENTITY);

template<typename T, Index N, typename RHS>
KOKKOS_INLINE_FUNCTION
RHS
solve_full_pivot(Tensor<T, N> const & A, RHS const & b);

//
// Swap row. Exchange rows i and j in place
// \param A tensor
// \param i index
// \param j index
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
swap_row(Tensor<T, N> & A, Index const i, Index const j)
{
  Index const
  dimension = A.get_dimension();

  if (i != j) {
    for (Index k = 0; k < dimension; ++k) {
      std::swap(A(i, k), A(j, k));
    }
  }
  return;
}

//
// Swap column. Exchange columns i and j in place
// \param A tensor
// \param i index
// \param j index
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
swap_col(Tensor<T, N> & A, Index const i, Index const j)
{
  Index const
  dimension = A.get_dimension();

  if (i != j) {
    for (Index k = 0; k < dimension; ++k) {
      std::swap(A(k, i), A(k, j));
    }
  }
  return;
}

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
///
/// Solve linear system A x = b via Gauss-Jordan elimination with full
/// pivoting.
///
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

    // The search finds no pivot if all remaining entries are zero or NaN
    // (NaN compares false against the running maximum). Indexing with
    // pivot_row == dimension would write out of bounds, so fail explicitly.
    if (pivot_row >= dimension || pivot_col >= dimension) {
      MT_ERROR_EXIT("Full-pivot solve failed: singular or non-finite matrix.");
    }

    // Gauss-Jordan elimination
    T const
    t{S(pivot_row, pivot_col)};

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

// Auxiliary functions for precondioners.
namespace impl {

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

} // namespace impl

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
    return impl::diagonal_precon(A, B);

  case PreconditionerType::MAX_ABS_ROW:
    return impl::maxabsrow_precon(A, B);
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

/// @}
} // namespace minitensor

#endif // MiniTensor_Inverse_h
