// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Tensor4_i_h)
#define MiniTensor_Tensor4_i_h

namespace minitensor
{

//
// 4th-order tensor constructor with NaNs
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4() :
TensorBase<T, Store>::TensorBase()
{
  set_dimension(N);
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Index const dimension) :
TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

//
// 4th-order tensor constructor with a specified value
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Filler const value) :
TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Index const dimension, Filler const value) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

//
//  Create 4th-order tensor from array
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(T const * data_ptr) :
TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Index const dimension, T const * data_ptr) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Tensor4<T, N> const & A) :
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

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Tensor<T, dimension_square<N>::value> const & A)
{
  Index const
  dimension_2nd = A.get_dimension();

  Index const
  dimension_4th = second_to_fourth_dimension(dimension_2nd);

  Tensor4<T, N> &
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
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::~Tensor4()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Tensor4<T, N>::get_dimension() const
{
  return TensorBase<T, Store>::get_dimension(ORDER);
}

//
// Set dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
Tensor4<T, N>::set_dimension(Index const dimension)
{
  TensorBase<T, Store>::set_dimension(dimension, ORDER);
  return;
}

//
// 4th-order tensor addition
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator+(Tensor4<S, N> const & A, Tensor4<T, N> const & B)
{
  Tensor4<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  add(A, B, C);

  return C;
}

//
// 4th-order tensor subtraction
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator-(Tensor4<S, N> const & A, Tensor4<T, N> const & B)
{
  Tensor4<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  subtract(A, B, C);

  return C;
}

//
// 4th-order tensor minus
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>
operator-(Tensor4<T, N> const & A)
{
  Tensor4<T, N>
  B(A.get_dimension());

  minus(A, B);

  return B;
}

//
// 4th-order equality
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor4<T, N> const & A, Tensor4<T, N> const & B)
{
  return equal(A, B);
}

//
// 4th-order inequality
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor4<T, N> const & A, Tensor4<T, N> const & B)
{
  return not_equal(A, B);
}

//
// Scalar 4th-order tensor product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T>, N>>::type
operator*(S const & s, Tensor4<T, N> const & A)
{
  Tensor4<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 4th-order tensor scalar product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T>, N>>::type
operator*(Tensor4<T, N> const & A, S const & s)
{
  Tensor4<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 4th-order tensor scalar division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator/(Tensor4<T, N> const & A, S const & s)
{
  Tensor4<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  divide(A, s, B);

  return B;
}

//
// 4th-order scalar tensor division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator/(S const & s, Tensor4<T, N> const & A)
{
  Tensor4<typename Promote<S, T>::type, N>
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
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T const &
Tensor4<T, N>::operator()(
    Index const i, Index const j, Index const k, Index const l) const
{
  Tensor4<T, N> const &
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
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T &
Tensor4<T, N>::operator()(
    Index const i, Index const j, Index const k, Index const l)
{
  Tensor4<T, N> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[((i * dimension + j) * dimension + k) * dimension + l];
}

//
// 4th-order inverse
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>
inverse(Tensor4<T, N> const & A)
{
  return Tensor4<T, N>(inverse(Tensor<T, dimension_square<N>::value>(A)));
}

} // namespace minitensor

#endif // MiniTensor_Tensor4_i_h
