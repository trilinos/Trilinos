// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_ArithTraits.hpp"
#if !defined(MiniTensor_Vector_i_h)
#define MiniTensor_Vector_i_h

namespace minitensor
{

//
// Default constructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector() :
    TensorBase<T, Store>::TensorBase()
{
  set_dimension(N);
  return;
}

//
// Constructor that initializes to NaNs
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(Index const dimension) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

//
// Create vector from a specified value
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(Index const dimension, Filler const value) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(Filler const value) :
    TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

//
// Create vector from array
//
template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    ArrayT & data,
    Index index1) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    ArrayT & data,
    Index index1,
    Index index2) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1, index2)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1, index2, index3)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4) :
    TensorBase<T, Store>::TensorBase(
        N,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4,
    Index index5) :
    TensorBase<T, Store>::TensorBase(
        N,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4,
        index5)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4,
    Index index5,
    Index index6) :
    TensorBase<T, Store>::TensorBase(
        N,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4,
        index5,
        index6)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    Index const dimension,
    ArrayT & data,
    Index index1) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data, index1)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    Index const dimension,
    ArrayT & data,
    Index index1,
    Index index2) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data, index1, index2)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    Index const dimension,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3) :
    TensorBase<T, Store>::TensorBase(
        dimension,
        ORDER,
        data,
        index1,
        index2,
        index3)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    Index const dimension,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4) :
    TensorBase<T, Store>::TensorBase(
        dimension,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    Index const dimension,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4,
    Index index5) :
    TensorBase<T, Store>::TensorBase(
        dimension,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4,
        index5)
{
  return;
}

template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(
    Source const,
    Index const dimension,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4,
    Index index5,
    Index index6) :
    TensorBase<T, Store>::TensorBase(
        dimension,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4,
        index5,
        index6)
{
  return;
}


template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(Index const dimension, T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(Vector<T, N> const & v) :
    TensorBase<T, Store>::TensorBase(v)
{
  return;
}

//
// Create vector specifying components
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(T const & s0, T const & s1)
{
  Vector <T, N> &
  self = (*this);

  self.set_dimension(2);

  self[0] = s0;
  self[1] = s1;

  return;
}

//
// Create vector specifying components
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::Vector(T const & s0, T const & s1, T const & s2)
{
  Vector <T, N> &
  self = (*this);

  self.set_dimension(3);

  self[0] = s0;
  self[1] = s1;
  self[2] = s2;

  return;
}

//
// Simple destructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>::~Vector()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Vector<T, N>::get_dimension() const
{
  return TensorBase<T, Store>::get_dimension(ORDER);
}

//
// Get number rows
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Vector<T, N>::get_num_rows() const
{
  return get_dimension();
}

//
// Get number columns
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Vector<T, N>::get_num_cols() const
{
  return static_cast<Index>(1);
}

//
// Set dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
Vector<T, N>::set_dimension(Index const dimension)
{
  TensorBase<T, Store>::set_dimension(dimension, ORDER);
  return;
}

//
// Indexing for constant vector
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T const &
Vector<T, N>::operator()(Index const i) const
{
  return (*this)[i];
}

//
// Indexing for vector
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T &
Vector<T, N>::operator()(Index const i)
{
  return (*this)[i];
}

//
// Indexing for constant vector (convenience)
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T const &
Vector<T, N>::operator()(Index const i, Index const) const
{
  return (*this)[i];
}

//
// Indexing for vector (convenience)
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T &
Vector<T, N>::operator()(Index const i, Index const)
{
  return (*this)[i];
}

//
// Vector addition
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator+(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Vector<typename Promote<S, T>::type, N>
  w(u.get_dimension());

  add(u, v, w);

  return w;
}

//
// Vector subtraction
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator-(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Vector<typename Promote<S, T>::type, N>
  w(u.get_dimension());

  subtract(u, v, w);

  return w;
}

//
// Vector minus
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
operator-(Vector<T, N> const & u)
{
  Vector<T, N>
  v(u.get_dimension());

  minus(u, v);

  return v;
}

//
// Vector dot product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
operator*(Vector<S, N> const & u, Vector<T, N> const & v)
{
  return dot(u, v);
}

//
// Vector equality tested by components
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION

bool
operator==(Vector<T, N> const & u, Vector<T, N> const & v)
{
  return equal(u, v);
}

//
// Vector inequality tested by components
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION

bool
operator!=(Vector<T, N> const & u, Vector<T, N> const & v)
{
  return not_equal(u, v);
}

//
// Scalar vector product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename
lazy_disable_if<order_1234<S>, apply_vector<Promote<S, T>, N>>::type
operator*(S const & s, Vector<T, N> const & u)
{
  Vector<typename Promote<S, T>::type, N>
  v(u.get_dimension());

  scale(u, s, v);

  return v;
}

//
// Vector scalar product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename
lazy_disable_if<order_1234<S>, apply_vector<Promote<S, T>, N>>::type
operator*(Vector<T, N> const & u, S const & s)
{
  Vector<typename Promote<S, T>::type, N>
  v(u.get_dimension());

  scale(u, s, v);

  return v;
}

//
// Vector scalar division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator/(Vector<T, N> const & u, S const & s)
{
  Vector<typename Promote<S, T>::type, N>
  v(u.get_dimension());

  divide(u, s, v);

  return v;
}

//
// Scalar vector division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator/(S const & s, Vector<T, N> const & u)
{
  Vector<typename Promote<S, T>::type, N>
  v(u.get_dimension());

  split(u, s, v);

  return v;
}

//
// Vector dot product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dot(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Index const
  dimension = u.get_dimension();

  assert(v.get_dimension() == dimension);

  typename Promote<S, T>::type
  s = 0.0;

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      s += u(i) * v(i);
    }
    break;

  case 3:
    s = u(0) * v(0) + u(1) * v(1) + u(2) * v(2);
    break;

  case 2:
    s = u(0) * v(0) + u(1) * v(1);
    break;

  }

  return s;
}

//
// Cross product only valid for R^3.
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
cross(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Index const
  dimension = u.get_dimension();

  assert(v.get_dimension() == dimension);

  Vector<typename Promote<S, T>::type, N>
  w(dimension);

  switch (dimension) {

  case 3:
    w(0) = u(1) * v(2) - u(2) * v(1);
    w(1) = u(2) * v(0) - u(0) * v(2);
    w(2) = u(0) * v(1) - u(1) * v(0);
    break;

  default:
    MT_ERROR_EXIT("Cross product defined for 3D only.");
    break;

  }

  return w;
}

//
// R^N vector 2-norm
// \return \f$ \sqrt{u \cdot u} \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm(Vector<T, N> const & u)
{
  T const
  s = norm_square(u);

  if (s > 0.0) return std::sqrt(s);

  return 0.0;
}

//
// R^N vector 2-norm square for fast distance calculations.
// \return \f$ u \cdot u \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_square(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

  default:
    s = dot(u, u);
    break;

  case 3:
    s = u(0) * u(0) + u(1) * u(1) + u(2) * u(2);
    break;

  case 2:
    s = u(0) * u(0) + u(1) * u(1);
    break;

  }

  return s;
}

//
// R^N vector 1-norm
// \return \f$ \sum_i |u_i| \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_1(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {
 
  default:
    for (Index i = 0; i < dimension; ++i) {
      s += std::abs(u(i));
    }
    break;

  case 3:
    s = std::abs(u(0)) + std::abs(u(1)) + std::abs(u(2));
    break;

  case 2:
    s = std::abs(u(0)) + std::abs(u(1));
    break;
  }

  return s;
}

//
// R^N vector infinity-norm
// \return \f$ \max(|u_0|,...|u_i|,...|u_N|) \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_infinity(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      s = max(minitensor::abs(u(i)), s);
    }
    break;

  case 3:
    s = max(max(minitensor::abs(u(0)), minitensor::abs(u(1))), minitensor::abs(u(2)));
    break;

  case 2:
    s = max(minitensor::abs(u(0)), minitensor::abs(u(1)));
    break;
  }

  return s;
}

//
// \return u / |u|, fails for |u| = 0
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
unit(Vector<T, N> const & u)
{
  return u / norm(u);
}

//
// Compute Householder vector
//
template <typename T, Index N>
std::pair<Vector<T, N>, T> house(Vector<T, N> const &x) {
  Vector<T, N>
  v = x;

  v[0] = 1.0;

  Index const
  dimension = x.get_dimension();

  T
  sigma = 0.0;

  for (Index i = 1; i < dimension; ++i) {
    sigma = v[i] * v[i];
  }

  T
  beta = 0.0;

  if (sigma > 0.0) {
    T const
    mu = std::sqrt(x[0] * x[0] + sigma);
    v[0] = x[0] > 0.0 ? -sigma / (x[0] + mu) : x[0] - mu;
    beta = 2.0 * v[0] * v[0] / (sigma + v[0] * v[0]);
    v = v / v[0];
  }

  return std::make_pair(v, beta);
}

} // namespace minitensor

#endif
 // MiniTensor_Vector_i_h
