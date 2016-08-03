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

#if !defined(Intrepid2_MiniTensor_Tensor4_i_h)
#define Intrepid2_MiniTensor_Tensor4_i_h

namespace Intrepid2
{

//
// 4th-order tensor constructor with NaNs
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>::Tensor4() :
TensorBase<T, Store>::TensorBase()
{
  set_dimension(N);
  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>::Tensor4(Index const dimension) :
TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

//
// 4th-order tensor constructor with a specified value
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>::Tensor4(ComponentValue const value) :
TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>::Tensor4(Index const dimension, ComponentValue const value) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

//
//  Create 4th-order tensor from array
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>::Tensor4(T const * data_ptr) :
TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>::Tensor4(Index const dimension, T const * data_ptr) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>::Tensor4(Tensor4<T, N, ES> const & A) :
TensorBase<T, Store>::TensorBase(A)
{
  return;
}

//
// 4th-order tensor from 2nd-order tensor
//

namespace {

KOKKOS_INLINE_FUNCTION
Index
second_to_fourth_dimension(Index const dimension_2nd)
{
  Index
  dimension_4th = 0;

  switch (dimension_2nd) {

  default:
    MT_ERROR_EXIT("Invalid dimension for 2nd-order tensor.");
    break;

  case 1:
    dimension_4th = 1;
    break;

  case 4:
    dimension_4th = 2;
    break;

  case 9:
    dimension_4th = 3;
    break;

  case 16:
    dimension_4th = 4;
    break;

  }

  return dimension_4th;
}

} //anonymous namespace

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>::Tensor4(Tensor<T, dimension_square<N>::value, ES> const & A)
{
  Index const
  dimension_2nd = A.get_dimension();

  Index const
  dimension_4th = second_to_fourth_dimension(dimension_2nd);

  Tensor4<T, N, ES> &
  self = (*this);

  self.set_dimension(dimension_4th);

  Index const
  number_components = dimension_2nd * dimension_2nd;

  for (Index i = 0; i < number_components; ++i) {
    self[i] = A[i];
  }

  return;
}
//
// 4th-order tensor simple destructor
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>::~Tensor4()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Index
Tensor4<T, N, ES>::get_dimension() const
{
  return TensorBase<T, Store>::get_dimension();
}

//
// Set dimension
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Tensor4<T, N, ES>::set_dimension(Index const dimension)
{
  if (IS_DYNAMIC == false) {
    assert(dimension <= N);
  }
  TensorBase<T, Store>::set_dimension(dimension, ORDER);
  return;
}

//
// 4th-order tensor addition
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N, ES>
operator+(Tensor4<S, N, ES> const & A, Tensor4<T, N, ES> const & B)
{
  Tensor4<typename Promote<S, T>::type, N, ES>
  C(A.get_dimension());

  add(A, B, C);

  return C;
}

//
// 4th-order tensor subtraction
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N, ES>
operator-(Tensor4<S, N, ES> const & A, Tensor4<T, N, ES> const & B)
{
  Tensor4<typename Promote<S, T>::type, N, ES>
  C(A.get_dimension());

  subtract(A, B, C);

  return C;
}

//
// 4th-order tensor minus
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>
operator-(Tensor4<T, N, ES> const & A)
{
  Tensor4<T, N, ES>
  B(A.get_dimension());

  minus(A, B);

  return B;
}

//
// 4th-order equality
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor4<T, N, ES> const & A, Tensor4<T, N, ES> const & B)
{
  return equal(A, B);
}

//
// 4th-order inequality
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor4<T, N, ES> const & A, Tensor4<T, N, ES> const & B)
{
  return not_equal(A, B);
}

//
// Scalar 4th-order tensor product
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T>, N, ES>>::type
operator*(S const & s, Tensor4<T, N, ES> const & A)
{
  Tensor4<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 4th-order tensor scalar product
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T>, N, ES>>::type
operator*(Tensor4<T, N, ES> const & A, S const & s)
{
  Tensor4<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 4th-order tensor scalar division
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N, ES>
operator/(Tensor4<T, N, ES> const & A, S const & s)
{
  Tensor4<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  divide(A, s, B);

  return B;
}

//
// 4th-order scalar tensor division
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N, ES>
operator/(S const & s, Tensor4<T, N, ES> const & A)
{
  Tensor4<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  split(A, s, B);

  return B;
}

//
// Indexing for constant 4th order tensor
// \param i index
// \param j index
// \param k index
// \param l index
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T const &
Tensor4<T, N, ES>::operator()(
    Index const i, Index const j, Index const k, Index const l) const
{
  Tensor4<T, N, ES> const &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[((i * dimension + j) * dimension + k) * dimension + l];
}

//
// 4th-order tensor indexing
// \param i index
// \param j index
// \param k index
// \param l index
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T &
Tensor4<T, N, ES>::operator()(
    Index const i, Index const j, Index const k, Index const l)
{
  Tensor4<T, N, ES> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[((i * dimension + j) * dimension + k) * dimension + l];
}

///
/// 4th-order inverse
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N, ES>
inverse(Tensor4<T, N, ES> const & A)
{
  return Tensor4<T, N, ES>(inverse(Tensor<T, dimension_square<N>::value, ES>(A)));
}

} // namespace Intrepid

#endif // Intrepid2_MiniTensor_Tensor4_i_h
