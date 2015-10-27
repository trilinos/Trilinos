// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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

#if !defined(Intrepid2_MiniTensor_Tensor3_i_h)
#define Intrepid2_MiniTensor_Tensor3_i_h

namespace Intrepid2 {

//
// 3rd-order tensor constructor with NaNs
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>::Tensor3() :
TensorBase<T, Store>::TensorBase()
{
  set_dimension(N);
  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>::Tensor3(Index const dimension) :
TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

//
// 3rd-order tensor constructor with a specified value
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>::Tensor3(ComponentValue const value) :
TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>::Tensor3(Index const dimension, ComponentValue const value) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

//
//  Create 3rd-order tensor from array
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>::Tensor3(T const * data_ptr) :
TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>::Tensor3(Index const dimension, T const * data_ptr) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>::Tensor3(Tensor3<T, N, ES> const & A) :
TensorBase<T, Store>::TensorBase(A)
{
  return;
}

//
// 3rd-order tensor simple destructor
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>::~Tensor3()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Index
Tensor3<T, N, ES>::get_dimension() const
{
  return TensorBase<T, Store>::get_dimension();
}

//
// Set dimension
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Tensor3<T, N, ES>::set_dimension(Index const dimension)
{
  if (IS_DYNAMIC == false) {
    assert(dimension <= N);
  }

  TensorBase<T, Store>::set_dimension(dimension, ORDER);

  return;
}

//
// 3rd-order tensor addition
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N,  ES>
operator+(Tensor3<S, N, ES>const & A, Tensor3<T, N, ES> const & B)
{
  Tensor3<typename Promote<S, T>::type, N, ES>
  C(A.get_dimension());

  add(A, B, C);

  return C;
}

//
// 3rd-order tensor subtraction
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
operator-(Tensor3<S, N, ES> const & A, Tensor3<T, N, ES> const & B)
{
  Tensor3<typename Promote<S, T>::type, N, ES>
  C(A.get_dimension());

  subtract(A, B, C);

  return C;
}

//
// 3rd-order tensor minus
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>
operator-(Tensor3<T, N, ES> const & A)
{
  Tensor3<T, N, ES>
  B(A.get_dimension());

  minus(A, B);

  return B;
}

//
// 3rd-order tensor equality
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor3<T, N, ES> const & A, Tensor3<T, N, ES> const & B)
{
  return equal(A, B);
}

//
// 3rd-order tensor inequality
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor3<T, N, ES> const & A, Tensor3<T, N, ES> const & B)
{
  return not_equal(A, B);
}

//
// Scalar 3rd-order tensor product
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N, ES> >::type
operator*(S const & s, Tensor3<T, N, ES> const & A)
{
  Tensor3<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 3rd-order tensor scalar product
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N, ES> >::type
operator*(Tensor3<T, N, ES> const & A, S const & s)
{
  Tensor3<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 3rd-order tensor scalar division
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
operator/(Tensor3<T, N, ES> const & A, S const & s)
{
  Tensor3<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  divide(A, s, B);

  return B;
}

//
// 3rd-order scalar tensor division
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
operator/(S const & s, Tensor3<T, N, ES> const & A)
{
  Tensor3<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  split(A, s, B);

  return B;
}

//
// Indexing for constant 3rd order tensor
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T const &
Tensor3<T, N, ES>::operator()(Index const i, Index const j, Index const k) const
{
  Tensor3<T, N, ES> const &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[(i * dimension + j) * dimension + k];
}

//
// 3rd-order tensor indexing
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T &
Tensor3<T, N, ES>::operator()(Index const i, Index const j, Index const k)
{
  Tensor3<T, N, ES> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[(i * dimension + j) * dimension + k];
}

// Local utility functions
namespace {

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void ones_in_diagonal(Tensor3<T, N, ES> & A)
{
  Index const
  dimension = A.get_dimension();

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      A(i, i, i) = 1.0;
    }
    break;

  case 3:
    A(0, 0, 0) = 1.0;
    A(1, 1, 1) = 1.0;
    A(2, 2, 2) = 1.0;
    break;

  case 2:
    A(0, 0, 0) = 1.0;
    A(1, 1, 1) = 1.0;
    break;

  }

  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void fill_levi_civita(Tensor3<T, N, ES> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        A(i, j, k) = levi_civita<T>(i, j, k);
      }
    }
  }

  return;
}

} // anonymous namespace

//
// Levi-Civita symbol
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
levi_civita_3()
{
  Tensor3<T, N, ES>
  A(N, ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC,ES> const
levi_civita_3(Index const dimension)
{
  Tensor3<T, DYNAMIC,ES>
  A(dimension, ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
levi_civita_3(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor3<T, DYNAMIC,ES>
  A(dimension, ZEROS);

  fill_levi_civita(A);

  return A;
}
// Permutation symbol
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
permutation_3()
{
  return levi_civita_3<T, N, ES>();
}

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC,ES> const
permutation_3(Index const dimension)
{
  return levi_civita_3<T>(dimension);
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
permutation_3(Index const dimension)
{
  return levi_civita_3<T, N, ES>(dimension);
}

//
// Alternating symbol
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
alternator_3()
{
  return levi_civita_3<T, N, ES>();
}

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC,ES> const
alternator_3(Index const dimension)
{
  return levi_civita_3<T>(dimension);
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
alternator_3(Index const dimension)
{
  return levi_civita_3<T, N, ES>(dimension);
}

} // namespace Intrepid

#endif // Intrepid2_MiniTensor_Tensor3_i_h
