// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Vector_h)
#define MiniTensor_Vector_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include <Kokkos_Core.hpp>
#include "MiniTensor_TensorBase.h"

namespace minitensor {

template<typename T, Index N>
using vector_store = Storage<T, dimension_power<N, 1>::value>;

///
/// Vector class.
///
template<typename T, Index N = DYNAMIC>
class Vector: public TensorBase<T, vector_store<T, N>>
{
public:

  ///
  /// Order
  ///
  static constexpr
  Index
  ORDER = 1;

  ///
  /// Static or dynamic
  ///
  static constexpr
  bool
  IS_DYNAMIC = N == DYNAMIC;

  ///
  /// Storage type
  ///
  using Store = vector_store<T, N>;

  ///
  /// Vector order
  ///
  KOKKOS_INLINE_FUNCTION
  static constexpr
  Index
  get_order()
  {
    return ORDER;
  }

  ///
  /// Construction that initializes to NaNs
  /// \param dimension the space dimension
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Vector();

  explicit
  KOKKOS_INLINE_FUNCTION
  Vector(Index const dimension);

  ///
  /// Create vector from a specified value
  /// \param dimension the space dimension
  /// \param value all components are set equal to this
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Vector(Filler const value);

  explicit
  KOKKOS_INLINE_FUNCTION
  Vector(Index const dimension, Filler const value);

  ///
  /// Create vector from array.
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      ArrayT & data,
      Index index1);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5,
      Index index6);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION  
  Vector(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5,
      Index index6);

  KOKKOS_INLINE_FUNCTION
  Vector(T const * data_ptr);

  KOKKOS_INLINE_FUNCTION
  Vector(Index const dimension, T const * data_ptr);

  ///
  /// Copy constructor
  ///
  KOKKOS_INLINE_FUNCTION
  Vector(Vector<T, N> const & v);

  ///
  /// Create vector specifying components
  /// \param s0 s1 are the vector components in the R^2 canonical basis
  ///
  KOKKOS_INLINE_FUNCTION
  Vector(T const & s0, T const & s1);

  ///
  /// Create vector specifying components
  /// the vector components in the R^3 canonical basis
  /// \param s0 s1 s2 are the vector components in the R^3 canonical basis
  ///
  KOKKOS_INLINE_FUNCTION
  Vector(T const & s0, T const & s1, T const & s2);

  ///
  /// Simple destructor
  ///
  virtual
  KOKKOS_INLINE_FUNCTION
  ~Vector();

  ///
  /// Indexing for constant vector
  /// \param i the index
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  operator()(Index const i) const;

  ///
  /// Vector indexing
  /// \param i the index
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  operator()(Index const i);

  ///
  /// const Vector indexing convenience operator for some algorithms
  /// \param i the index
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  operator()(Index const i, Index const) const;

  ///
  /// Vector indexing convenience operator for some algorithms
  /// \param i the index
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  operator()(Index const i, Index const);

  ///
  /// \return dimension
  ///
  KOKKOS_INLINE_FUNCTION
  Index
  get_dimension() const;

  ///
  /// \return number rows
  ///
  KOKKOS_INLINE_FUNCTION
  Index
  get_num_rows() const;

  ///
  /// \return number rows
  ///
  KOKKOS_INLINE_FUNCTION
  Index
  get_num_cols() const;

  ///
  /// \param dimension of vector
  ///
  KOKKOS_INLINE_FUNCTION
  void
  set_dimension(Index const dimension);

};

///
/// Vector addition
/// \return \f$ u + v \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator+(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Vector substraction
/// \return \f$ u - v \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator-(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Vector minus
/// \return \f$ -u \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
operator-(Vector<T, N> const & u);

///
/// Vector dot product
/// \return \f$ u \cdot v \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
operator*(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Vector equality tested by components
/// \return \f$ u \equiv v \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator==(Vector<T, N> const & u, Vector<T, N> const & v);

///
/// Vector inequality tested by components
/// \return \f$ u \neq v \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Vector<T, N> const & u, Vector<T, N> const & v);

///
/// Scalar vector product
/// \param s scalar factor
/// \param u vector factor
/// \return \f$ s u \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_vector<Promote<S, T>, N>>::type
operator*(S const & s, Vector<T, N> const & u);

///
/// Vector scalar product
/// \param u vector factor
/// \param s scalar factor
/// \return \f$ s u \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_vector<Promote<S, T>, N>>::type
operator*(Vector<T, N> const & u, S const & s);

///
/// Vector scalar division
/// \param u vector
/// \param s scalar that divides each component of vector
/// \return \f$ u / s \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator/(Vector<T, N> const & u, S const & s);

///
/// Scalar vector division
/// \param s scalar
/// \param u vector that divides scalar with each component
/// \return \f$ s / u \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator/(S const & s, Vector<T, N> const & u);

///
/// Vector dot product
/// \return \f$ u \cdot v \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dot(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Cross product only valid for R^3.
/// R^N with N != 3 will produce an error.
/// \return \f$ u \times v \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
cross(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Vector 2-norm
/// \return \f$ \sqrt{u \cdot u} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm(Vector<T, N> const & u);

///
/// Vector 2-norm square. Used for fast distance calculation.
/// \return \f$ u \cdot u \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_square(Vector<T, N> const & u);

///
/// Vector 1-norm
/// \return \f$ |u_0|+|u_1|+|u_2| \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_1(Vector<T, N> const & u);

///
/// Vector infinity-norm
/// \return \f$ \max(|u_0|,|u_1|,|u_2|) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_infinity(Vector<T, N> const & u);

///
/// \return u / |u|, fails for |u| = 0
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
unit(Vector<T, N> const & u);

///
/// Compute Householder vector
/// See algorithm 5.1.1, Matrix Computations 3rd ed, Golub - Van Loan
/// \f$ Px = |X|e_1, P := I - \beta v v^T\f$
/// \return v, beta
///
template<typename T, Index N>
std::pair<Vector<T, N>, T>
house(Vector<T, N> const & x);

///
/// Vector input
/// \param u vector
/// \param is input stream
/// \return is input stream
///
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Vector<T, N> & u);

///
/// Vector output
/// \param u vector
/// \param os output stream
/// \return os output stream
///
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Vector<T, N> const & u);

} // namespace minitensor

#include "MiniTensor_Vector.i.h"
#include "MiniTensor_Vector.t.h"

#endif //MiniTensor_Vector_h
