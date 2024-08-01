// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Tensor3_i_h)
#define MiniTensor_Tensor3_i_h

namespace minitensor {

//
// 3rd-order tensor constructor with NaNs
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3() :
TensorBase<T, Store>::TensorBase()
{
  set_dimension(N);
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Index const dimension) :
TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

//
// 3rd-order tensor constructor with a specified value
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Filler const value) :
TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Index const dimension, Filler const value) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

//
//  Create 3rd-order tensor from array
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(T const * data_ptr) :
TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Index const dimension, T const * data_ptr) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Tensor3<T, N> const & A) :
TensorBase<T, Store>::TensorBase(A)
{
  return;
}

//
// 3rd-order tensor simple destructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::~Tensor3()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Tensor3<T, N>::get_dimension() const
{
  return TensorBase<T, Store>::get_dimension(ORDER);
}

//
// Set dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
Tensor3<T, N>::set_dimension(Index const dimension)
{
  TensorBase<T, Store>::set_dimension(dimension, ORDER);
  return;
}

//
// 3rd-order tensor addition
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
operator+(Tensor3<S, N>const & A, Tensor3<T, N> const & B)
{
  Tensor3<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  add(A, B, C);

  return C;
}

//
// 3rd-order tensor subtraction
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
operator-(Tensor3<S, N> const & A, Tensor3<T, N> const & B)
{
  Tensor3<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  subtract(A, B, C);

  return C;
}

//
// 3rd-order tensor minus
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>
operator-(Tensor3<T, N> const & A)
{
  Tensor3<T, N>
  B(A.get_dimension());

  minus(A, B);

  return B;
}

//
// 3rd-order tensor equality
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor3<T, N> const & A, Tensor3<T, N> const & B)
{
  return equal(A, B);
}

//
// 3rd-order tensor inequality
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor3<T, N> const & A, Tensor3<T, N> const & B)
{
  return not_equal(A, B);
}

//
// Scalar 3rd-order tensor product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N>>::type
operator*(S const & s, Tensor3<T, N> const & A)
{
  Tensor3<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 3rd-order tensor scalar product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N>>::type
operator*(Tensor3<T, N> const & A, S const & s)
{
  Tensor3<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 3rd-order tensor scalar division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
operator/(Tensor3<T, N> const & A, S const & s)
{
  Tensor3<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  divide(A, s, B);

  return B;
}

//
// 3rd-order scalar tensor division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
operator/(S const & s, Tensor3<T, N> const & A)
{
  Tensor3<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  split(A, s, B);

  return B;
}

//
// Indexing for constant 3rd order tensor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T const &
Tensor3<T, N>::operator()(Index const i, Index const j, Index const k) const
{
  Tensor3<T, N> const &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[(i * dimension + j) * dimension + k];
}

//
// 3rd-order tensor indexing
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T &
Tensor3<T, N>::operator()(Index const i, Index const j, Index const k)
{
  Tensor3<T, N> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[(i * dimension + j) * dimension + k];
}

// Local utility functions
namespace {

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void ones_in_diagonal(Tensor3<T, N> & A)
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

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void fill_levi_civita(Tensor3<T, N> & A)
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
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
levi_civita_3()
{
  Tensor3<T, N>
  A(N, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC> const
levi_civita_3(Index const dimension)
{
  Tensor3<T, DYNAMIC>
  A(dimension, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
levi_civita_3(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor3<T, DYNAMIC>
  A(dimension, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}
// Permutation symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
permutation_3()
{
  return levi_civita_3<T, N>();
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC> const
permutation_3(Index const dimension)
{
  return levi_civita_3<T>(dimension);
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
permutation_3(Index const dimension)
{
  return levi_civita_3<T, N>(dimension);
}

//
// Alternating symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
alternator_3()
{
  return levi_civita_3<T, N>();
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC> const
alternator_3(Index const dimension)
{
  return levi_civita_3<T>(dimension);
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
alternator_3(Index const dimension)
{
  return levi_civita_3<T, N>(dimension);
}

} // namespace minitensor

#endif // MiniTensor_Tensor3_i_h
