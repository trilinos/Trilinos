// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Norms_h)
#define MiniTensor_Norms_h

// Scalar measures of tensors: norms, invariants, argmax queries.
#include "MiniTensor_Tensor.h"

namespace minitensor {

/// \addtogroup minitensor_norms
/// @{

///
/// Tensor Frobenius norm
/// \return \f$ \sqrt{A:A} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm(Tensor<T, N> const & A);

///
/// Tensor 1-norm
/// \return \f$ \max_{j \in {0,1,2}}\Sigma_{i=0}^2 |A_{ij}| \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_1(Tensor<T, N> const & A);

///
/// Tensor infinity-norm
/// \return \f$ \max_{i \in {0,1,2}}\Sigma_{j=0}^2 |A_{ij}| \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_infinity(Tensor<T, N> const & A);

///
/// Subtensor
/// \param A tensor
/// \param i index
/// \param j index
/// \return Subtensor with i-row and j-col deleted.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
subtensor(Tensor<T, N> const & A, Index const i, Index const j);

///
/// Determinant
/// \return \f$ \det A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
det(Tensor<T, N> const & A);

///
/// Trace
/// \return \f$ A:I \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
trace(Tensor<T, N> const & A);

///
/// First invariant, trace
/// \return \f$ I_A = A:I \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
I1(Tensor<T, N> const & A);

///
/// Second invariant
/// \return \f$ II_A = \frac{1}{2}((I_A)^2-I_{A^2}) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
I2(Tensor<T, N> const & A);

///
/// Third invariant
/// \return \f$ III_A = \det A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
I3(Tensor<T, N> const & A);

///
/// Off-diagonal norm. Useful for SVD and other algorithms
/// that rely on Jacobi-type procedures.
/// \return \f$ \sqrt(\sum_i \sum_{j, j\neq i} a_{ij}^2) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_off_diagonal(Tensor<T, N> const & A);

///
/// Arg max abs. Useful for inverse and other algorithms
/// that rely on Jacobi-type procedures.
/// \return \f$ (p,q) = arg max_{i,j} |a_{ij}| \f$
///
template <typename T, Index N>
std::pair<Index, Index> arg_max_abs(Tensor<T, N> const &A);

///
/// Arg max off-diagonal. Useful for SVD and other algorithms
/// that rely on Jacobi-type procedures.
/// \return \f$ (p,q) = arg max_{i \neq j} |a_{ij}| \f$
///
template <typename T, Index N>
std::pair<Index, Index> arg_max_off_diagonal(Tensor<T, N> const &A);

//
// R^N tensor Frobenius norm
// \return \f$ \sqrt{A:A} \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

    default:
      s = norm_f_square(A);
      break;

    case 3:
      s+= A(0,0)*A(0,0) + A(0,1)*A(0,1) + A(0,2)*A(0,2);
      s+= A(1,0)*A(1,0) + A(1,1)*A(1,1) + A(1,2)*A(1,2);
      s+= A(2,0)*A(2,0) + A(2,1)*A(2,1) + A(2,2)*A(2,2);
      break;

    case 2:
      s+= A(0,0)*A(0,0) + A(0,1)*A(0,1);
      s+= A(1,0)*A(1,0) + A(1,1)*A(1,1);
      break;

    case 1:
      s+= A(0,0)*A(0,0);
      break;

  }

  if (s > 0.0) return std::sqrt(s);
  return 0.0;

}

//
// R^N tensor 1-norm
// \return \f$ \max_{j \in {0,\cdots,N}}\Sigma_{i=0}^N |A_{ij}| \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_1(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Vector<T, N>
  v(dimension);

  T
  s = 0.0;

  switch (dimension) {

  default:

    for (Index i = 0; i < dimension; ++i) {
      T t = 0.0;
      for (Index j = 0; j < dimension; ++j) {
        t += minitensor::abs(A(j, i));
      }
      v(i) = t;
    }

    for (Index i = 0; i < dimension; ++i) {
      s = max(s, v(i));
    }
    break;

  case 3:
    v(0) = minitensor::abs(A(0, 0)) + minitensor::abs(A(1, 0)) + minitensor::abs(A(2, 0));
    v(1) = minitensor::abs(A(0, 1)) + minitensor::abs(A(1, 1)) + minitensor::abs(A(2, 1));
    v(2) = minitensor::abs(A(0, 2)) + minitensor::abs(A(1, 2)) + minitensor::abs(A(2, 2));

    s = max(max(v(0), v(1)), v(2));
    break;

  case 2:
    v(0) = minitensor::abs(A(0, 0)) + minitensor::abs(A(1, 0));
    v(1) = minitensor::abs(A(0, 1)) + minitensor::abs(A(1, 1));

    s = max(v(0), v(1));
    break;

  case 1:
    s = minitensor::abs(A(0, 0));
    break;
  }

  return s;
}

//
// R^N tensor infinity-norm
// \return \f$ \max_{i \in {0,\cdots,N}}\Sigma_{j=0}^N |A_{ij}| \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_infinity(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Vector<T, N>
  v(dimension);

  T s = 0.0;

  switch (dimension) {

    default:
      for (Index i = 0; i < dimension; ++i) {
        T t = 0.0;
        for (Index j = 0; j < dimension; ++j) {
          t += minitensor::abs(A(i, j));
        }
        v(i) = t;
      }

      for (Index i = 0; i < dimension; ++i) {
        s = max(s, v(i));
      }
      break;

    case 3:
      v(0) = minitensor::abs(A(0, 0)) + minitensor::abs(A(0, 1)) + minitensor::abs(A(0, 2));
      v(1) = minitensor::abs(A(1, 0)) + minitensor::abs(A(1, 1)) + minitensor::abs(A(1, 2));
      v(2) = minitensor::abs(A(2, 0)) + minitensor::abs(A(2, 1)) + minitensor::abs(A(2, 2));

      s = max(max(v(0), v(1)), v(2));
      break;

    case 2:
      v(0) = minitensor::abs(A(0, 0)) + minitensor::abs(A(0, 1));
      v(1) = minitensor::abs(A(1, 0)) + minitensor::abs(A(1, 1));

      s = max(v(0), v(1));
      break;

    case 1:
      s = minitensor::abs(A(0, 0));
      break;

  }

  return s;
}

//
// R^N determinant
// Laplace expansion. Warning: no pivoting.
// Casual use only. Use Teuchos LAPACK interface for
// more efficient and robust techniques.
// \param A tensor
// \return \f$ \det A \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
det(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

    default:
    {
      int sign = 1;
      for (Index i = 0; i < dimension; ++i) {
        const T d = det(subtensor(A, i, 1));
        s += sign * d * A(i, 1);
        sign *= -1;
      }
    }
    break;

    case 3:
      s = -A(0,2)*A(1,1)*A(2,0) + A(0,1)*A(1,2)*A(2,0) +
           A(0,2)*A(1,0)*A(2,1) - A(0,0)*A(1,2)*A(2,1) -
           A(0,1)*A(1,0)*A(2,2) + A(0,0)*A(1,1)*A(2,2);
      break;

    case 2:
      s = A(0,0) * A(1,1) - A(1,0) * A(0,1);
      break;

    case 1:
      s = A(0,0);
      break;

  }

  return s;
}

//
// R^N trace
// \param A tensor
// \return \f$ A:I \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
trace(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  T s = 0.0;

  switch (dimension) {

    default:
      for (Index i = 0; i < dimension; ++i) {
        s += A(i,i);
      }
      break;

    case 3:
      s = A(0,0) + A(1,1) + A(2,2);
      break;

    case 2:
      s = A(0,0) + A(1,1);
      break;

    case 1:
      s = A(0,0);
      break;

 }

  return s;
}

//
// R^N first invariant, trace
// \param A tensor
// \return \f$ I_A = A:I \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
I1(Tensor<T, N> const & A)
{
  return trace(A);
}

//
// R^N second invariant
// \param A tensor
// \return \f$ II_A \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
I2(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  T
  s = 0.0;

  T const
  trA = trace(A);

  switch (dimension) {

    default:
#ifdef KOKKOS_ENABLE_CUDA
      Kokkos::abort("I2 for N > 3 not implemented.");
      return T();
#else
      std::cerr << "I2 for N > 3 not implemented." << std::endl;
      exit(1);
#endif
      break;

    case 3:
      s = 0.5 * (trA*trA - A(0,0)*A(0,0) - A(1,1)*A(1,1) - A(2,2)*A(2,2)) -
      A(0,1)*A(1,0) - A(0,2)*A(2,0) - A(1,2)*A(2,1);
      break;

    case 2:
      s = - det(A);
      break;

    case 1:
      s = 0.0;
      break;

  }

  return s;
}

//
// R^N third invariant
// \param A tensor
// \return \f$ III_A \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
I3(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

    default:
#ifdef KOKKOS_ENABLE_CUDA
      Kokkos::abort("I3 for N > 3 not implemented.");
      return T();
#else
      std::cerr << "I3 for N > 3 not implemented." << std::endl;
      exit(1);
#endif
      break;

    case 3:
      s = det(A);
      break;

    case 2:
      s = 0.0;
      break;

    case 1:
      s = 0.0;
      break;

  }

  return s;
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

/// @}
} // namespace minitensor

#endif // MiniTensor_Norms_h
