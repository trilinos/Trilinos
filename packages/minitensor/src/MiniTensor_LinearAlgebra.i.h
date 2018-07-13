// @HEADER
// ************************************************************************
//
//                           MiniTensor Package
//                 Copyright (2016) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(MiniTensor_LinearAlgebra_i_h)
#define MiniTensor_LinearAlgebra_i_h

namespace minitensor {

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
          t += std::abs(A(j, i));
        }
        v(i) = t;
      }

      for (Index i = 0; i < dimension; ++i) {
        s = std::max(s, v(i));
      }
      break;

    case 3:
      v(0) = std::abs(A(0,0)) + std::abs(A(1,0)) + std::abs(A(2,0));
      v(1) = std::abs(A(0,1)) + std::abs(A(1,1)) + std::abs(A(2,1));
      v(2) = std::abs(A(0,2)) + std::abs(A(1,2)) + std::abs(A(2,2));

      s = std::max(std::max(v(0),v(1)),v(2));
      break;

    case 2:
      v(0) = std::abs(A(0,0)) + std::abs(A(1,0));
      v(1) = std::abs(A(0,1)) + std::abs(A(1,1));

      s = std::max(v(0),v(1));
      break;

    case 1:
      s = std::abs(A(0,0));
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
          t += std::abs(A(i, j));
        }
        v(i) = t;
      }

      for (Index i = 0; i < dimension; ++i) {
        s = std::max(s, v(i));
      }
      break;

    case 3:
      v(0) = std::abs(A(0,0)) + std::abs(A(0,1)) + std::abs(A(0,2));
      v(1) = std::abs(A(1,0)) + std::abs(A(1,1)) + std::abs(A(1,2));
      v(2) = std::abs(A(2,0)) + std::abs(A(2,1)) + std::abs(A(2,2));

      s = std::max(std::max(v(0),v(1)),v(2));
      break;

    case 2:
      v(0) = std::abs(A(0,0)) + std::abs(A(0,1));
      v(1) = std::abs(A(1,0)) + std::abs(A(1,1));

      s = std::max(v(0),v(1));
      break;

    case 1:
      s = std::abs(A(0,0));
      break;

  }

  return s;
}

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
#else
      std::cerr << "I2 for N > 3 not implemented." << std::endl;
#endif
      exit(1);
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
#else 
      std::cerr << "I3 for N > 3 not implemented." << std::endl;
#endif
      exit(1);
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
// Condition number.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
cond(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N> const
  S = boost::get<1>(svd(A));

  T const
  k = S(0, 0) / S(dimension - 1, dimension - 1);

  return k;
}

//
// Reciprocal condition number.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
inv_cond(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N> const
  S = boost::get<1>(svd(A));

  T const
  k = S(dimension - 1, dimension - 1) / S(0, 0);

  return k;
}

} // namespace minitensor

#endif // MiniTensor_LinearAlgebra_i_h
