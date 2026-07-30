// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Tensor_h)
#define MiniTensor_Tensor_h

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Kokkos_Core.hpp>
#include "MiniTensor_Vector.h"

namespace minitensor {

template<typename T, Index N>
using tensor_store = Storage<T, dimension_power<N, 2>::value>;

///
/// Second order tensor.
///
template<typename T, Index N = DYNAMIC>
class Tensor: public TensorBase<T, tensor_store<T, N>>
{
public:

  ///
  /// Order
  ///
  static constexpr
  Index
  ORDER = 2;

  ///
  /// Static or dynamic
  ///
  static constexpr
  bool
  IS_DYNAMIC = N == DYNAMIC;

  ///
  /// Storage type
  ///
  using Store = tensor_store<T, N>;

  ///
  /// Tensor order
  ///
  KOKKOS_INLINE_FUNCTION
  static constexpr
  Index
  get_order()
  {
    return ORDER;
  }

  ///
  /// Constructor that initializes to NaNs
  /// \param dimension the space dimension
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor();

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(Index const dimension);

  ///
  /// Create tensor from a specified value
  /// \param dimension the space dimension
  /// \param value all components are set equal to this
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(Filler const value);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(Index const dimension, Filler const value);

  ///
  /// Create tensor from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      ArrayT & data,
      Index index1);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5,
      Index index6);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5);

  template<class ArrayT>
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Source const source,
      Index const dimension,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5,
      Index index6);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(T const * data_ptr);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(Index const dimension, T const * data_ptr);

  ///
  /// Copy constructor
  ///
  KOKKOS_INLINE_FUNCTION
  Tensor(Tensor<T, N> const & A);

  ///
  /// 2nd-order tensor from 4th-order tensor
  ///
  KOKKOS_INLINE_FUNCTION
  Tensor(Tensor4<T, dimension_sqrt<N>::value> const & A);

  ///
  /// Create tensor specifying components
  /// \param  s00 s01 ... components in the R^2 canonical basis
  ///
  //
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(T const & s00, T const & s01, T const & s10, T const & s11);

  ///
  /// Create tensor specifying components
  /// \param  s00 s01 ... components in the R^3 canonical basis
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      T const & s00, T const & s01, T const & s02,
      T const & s10, T const & s11, T const & s12,
      T const & s20, T const & s21, T const & s22);

  ///
  /// Create tensor from array
  /// \param data_ptr pointer into the array
  /// \param component_order component convention (3D only)
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(T const * data_ptr, ComponentOrder const component_order);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Index const dimension,
      T const * data_ptr,
      ComponentOrder const component_order);

  ///
  /// Simple destructor
  ///
  virtual
  KOKKOS_INLINE_FUNCTION
  ~Tensor();

  ///
  /// Indexing for constant tensor
  /// \param i index
  /// \param j index
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  operator()(Index const i, Index const j) const;

  ///
  /// Tensor indexing
  /// \param i index
  /// \param j index
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  operator()(Index const i, Index const j);

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
/// Tensor addition
/// \return \f$ A + B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator+(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor subtraction
/// \return \f$ A - B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator-(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor minus
/// \return \f$ -A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
operator-(Tensor<T, N> const & A);

///
/// Tensor equality
/// Tested by components
/// \return \f$ A \equiv B \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor<T, N> const & A, Tensor<T, N> const & B);

///
/// Tensor inequality
/// Tested by components
/// \return \f$ A \neq B \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor<T, N> const & A, Tensor<T, N> const & B);

///
/// Tensor vector product v = A u
/// \param A tensor
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator*(Tensor<T, N> const & A, Vector<S, N> const & u);

///
/// Vector tensor product v = u A
/// \param A tensor
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator*(Vector<S, N> const & u, Tensor<T, N> const & A);

///
/// Tensor dot product C = A B
/// \return \f$ A \cdot B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator*(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Scalar tensor product
/// \param s scalar
/// \param A tensor
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N>>::type
operator*(S const & s, Tensor<T, N> const & A);

///
/// Tensor scalar product
/// \param A tensor
/// \param s scalar
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N>>::type
operator*(Tensor<T, N> const & A, S const & s);

///
/// Tensor scalar division
/// \param A tensor
/// \param s scalar
/// \return \f$ A / s \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator/(Tensor<T, N> const & A, S const & s);

///
/// Scalar tensor division
/// \param s scalar
/// \param A tensor that divides scalar with each component
/// \return \f$ s / A \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator/(S const & s, Tensor<T, N> const & A);

///
/// Tensor input
/// \param A tensor
/// \param is input stream
/// \return is input stream
///
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Tensor<T, N> & A);

///
/// Tensor output
/// \param A tensor
/// \param os output stream
/// \return os output stream
///
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Tensor<T, N> const & A);

///
/// Extract a row as a vector
/// \param A tensor
/// \param i index of row
/// \return \f$ v = A(i,:) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
row(Tensor<T, N> const & A, Index const i);

///
/// Extract a column as a vector
/// \param A tensor
/// \param j index of column
/// \return \f$ v = A(:,j) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
col(Tensor<T, N> const & A, Index const j);

///
/// Tensor vector product v = A u
/// \param A tensor
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
dot(Tensor<T, N> const & A, Vector<S, N> const & u);

///
/// Vector tensor product v = u A
/// \param A tensor
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
dot(Vector<S, N> const & u, Tensor<T, N> const & A);

///
/// Tensor tensor product C = A B
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A \cdot B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor tensor product C = A^T B
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A^T \cdot B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
t_dot(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor matrix product C = A^T B
/// \param A tensor
/// \param B matrix
/// \return a matrix \f$ A^T \cdot B \f$
///
template<typename S, typename T, Index N, Index P>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, N, P>
t_dot(Tensor<S, N> const & A, Matrix<T, N, P> const & B);

///
/// Tensor vector product C = A^T B
/// \param A tensor
/// \param B vector
/// \return a vector \f$ A^T \cdot B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
t_dot(Tensor<S, N> const & A, Vector<T, N> const & B);

///
/// Tensor tensor product C = A B^T
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A \cdot B^T \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot_t(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor tensor product C = A^T B^T
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A^T \cdot B^T \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
t_dot_t(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor tensor double dot product (contraction)
/// \param A tensor
/// \param B tensor
/// \return a scalar \f$ A : B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dotdot(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Dyad
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dyad(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Bun operator, just for Jay, and now Reese too.
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
bun(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Tensor product
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
tensor(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Diagonal tensor from vector
/// \param v vector
/// \return A = diag(v)
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
diag(Vector<T, N> const & v);

///
/// Diagonal of tensor in a vector
/// \param A tensor
/// \return v = diag(A)
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
diag(Tensor<T, N> const & A);

///
/// Zero 2nd-order tensor
/// All components are zero
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION 
Tensor<T, N> const
zero();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
zero(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
zero(Index const dimension);

///
/// 2nd-order identity tensor
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
identity();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
identity(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
identity(Index const dimension);

///
/// 2nd-order identity tensor, à la Matlab
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
eye();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
eye(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
eye(Index const dimension);

///
/// Levi-Civita symbol
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
levi_civita_2();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
levi_civita_2(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
levi_civita_2(Index const dimension);

///
/// Permutation symbol
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
permutation_2();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
permutation_2(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
permutation_2(Index const dimension);

///
/// Alternating symbol
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
alternator_2();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
alternator_2(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
alternator_2(Index const dimension);

///
/// 2nd-order tensor transpose
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
transpose(Tensor<T, N> const & A);

///
/// C^N 2nd-order tensor adjoint
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
adjoint(Tensor<T, N> const & A);

///
/// Symmetric part of 2nd-order tensor
/// \return \f$ \frac{1}{2}(A + A^T) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
sym(Tensor<T, N> const & A);

///
/// Skew symmetric part of 2nd-order tensor
/// \return \f$ \frac{1}{2}(A - A^T) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
skew(Tensor<T, N> const & A);

///
/// Skew symmetric 2nd-order tensor from vector valid for R^3 only.
/// R^N with N != 3 will produce an error
/// \param u vector
/// \return \f$ {{0, -u_2, u_1}, {u_2, 0, -u_0}, {-u_1, u+0, 0}} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
skew(Vector<T, N> const & u);

} // namespace minitensor

namespace minitensor {

//
// Constructor that initializes to NaNs
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor() :
    TensorBase<T, Store>::TensorBase()
{
  set_dimension(N);
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(Index const dimension) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

//
// Create tensor from a specified value
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(Filler const value) :
    TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(Index const dimension, Filler const value) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

//
//  Create tensor from array
//
//
template<typename T, Index N>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(
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
Tensor<T, N>::Tensor(T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(Index const dimension, T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}
//
// Copy constructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(Tensor<T, N> const & A) :
    TensorBase<T, Store>::TensorBase(A)
{
  return;
}

//
// Create tensor specifying components
// \param  s00, s01, ... components in the R^2 canonical basis
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(
    T const & s00, T const & s01,
    T const & s10, T const & s11)
{
  Tensor<T, N> &
  self = (*this);

  self.set_dimension(2);

  self[0] = s00;
  self[1] = s01;

  self[2] = s10;
  self[3] = s11;

  return;
}

//
// Create tensor specifying components
// \param  s00, s01, ... components in the R^3 canonical basis
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(
    T const & s00, T const & s01, T const & s02,
    T const & s10, T const & s11, T const & s12,
    T const & s20, T const & s21, T const & s22)
{
  Tensor<T, N> &
  self = (*this);

  self.set_dimension(3);

  self[0] = s00;
  self[1] = s01;
  self[2] = s02;

  self[3] = s10;
  self[4] = s11;
  self[5] = s12;

  self[6] = s20;
  self[7] = s21;
  self[8] = s22;

  return;
}

//
//  Create tensor from array with component order
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(
    T const * data_ptr,
    ComponentOrder const component_order)
{
  assert(data_ptr != NULL);

  fill(data_ptr, component_order);

  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(
    Index const dimension,
    T const * data_ptr,
    ComponentOrder const component_order)
{
  assert(data_ptr != NULL);

  Tensor<T, N> &
  self = (*this);

  self.set_dimension(dimension);

  fill(data_ptr, component_order);

  return;
}

//
// 2nd-order tensor from 4th-order tensor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::Tensor(Tensor4<T, dimension_sqrt<N>::value> const & A)
{
  Index const
  dimension_4th = A.get_dimension();

  Index const
  dimension_2nd = dimension_4th * dimension_4th;

  Tensor<T, N> &
  self = (*this);

  self.set_dimension(dimension_2nd);

  Index const
  number_components = dimension_2nd * dimension_2nd;

  for (Index i = 0; i < number_components; ++i) {
    self[i] = A[i];
  }

  return;
}

//
// Simple destructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>::~Tensor()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Tensor<T, N>::get_dimension() const
{
  return TensorBase<T, Store>::get_dimension(ORDER);
}

//
// Get number rows
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Tensor<T, N>::get_num_rows() const
{
  return get_dimension();
}

//
// Get number columns
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Tensor<T, N>::get_num_cols() const
{
  return get_dimension();
}

//
// Set dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N>::set_dimension(Index const dimension)
{
  TensorBase<T, Store>::set_dimension(dimension, ORDER);
  return;
}

//
// Indexing for constant tensor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T const &
Tensor<T, N>::operator()(Index const i, Index const j) const
    {
  Tensor<T, N> const &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[i * dimension + j];
}

//
//Tensor indexing
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T &
Tensor<T, N>::operator()(Index const i, Index const j)
{
  Tensor<T, N> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[i * dimension + j];
}

namespace {

template<typename S>
KOKKOS_INLINE_FUNCTION
bool
greater_than(S const & a, S const & b)
{
  return a.first > b.first;
}

} // anonymous namespace

//
// Extract a row as a vector
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
row(Tensor<T, N> const & A, Index const i)
{
  Index const
  dimension = A.get_dimension();

  Vector<T, N>
  v(dimension);

  switch (dimension) {
  default:
    for (Index j = 0; j < dimension; ++j) {
      v(j) = A(i, j);
    }
    break;

  case 2:
    v(0) = A(i, 0);
    v(1) = A(i, 1);
    break;

  case 3:
    v(0) = A(i, 0);
    v(1) = A(i, 1);
    v(2) = A(i, 2);
    break;
  }

  return v;
}

//
// Extract a column as a vector
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
col(Tensor<T, N> const & A, Index const j)
{
  Index const
  dimension = A.get_dimension();

  Vector<T, N>
  v(dimension);

  switch (dimension) {
  default:
    for (Index i = 0; i < dimension; ++i) {
      v(i) = A(i, j);
    }
    break;

  case 2:
    v(0) = A(0, j);
    v(1) = A(1, j);
    break;

  case 3:
    v(0) = A(0, j);
    v(1) = A(1, j);
    v(2) = A(2, j);
    break;
  }

  return v;
}

//
// Tensor addition
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator+(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Tensor<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  add(A, B, C);

  return C;
}

//
// Tensor subtraction
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator-(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Tensor<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  subtract(A, B, C);

  return C;
}

//
// Tensor minus
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
operator-(Tensor<T, N> const & A)
{
  Tensor<T, N>
  B(A.get_dimension());

  minus(A, B);

  return B;
}

//
// Tensor equality
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor<T, N> const & A, Tensor<T, N> const & B)
{
  return equal(A, B);
}

//
// Tensor inequality
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor<T, N> const & A, Tensor<T, N> const & B)
{
  return not_equal(A, B);
}

//
// Scalar tensor product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename
lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N>>::type
operator*(S const & s, Tensor<T, N> const & A)
{
  Tensor<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// Tensor scalar product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename
lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N>>::type
operator*(Tensor<T, N> const & A, S const & s)
{
  Tensor<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// Tensor scalar division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator/(Tensor<T, N> const & A, S const & s)
{
  Tensor<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  divide(A, s, B);

  return B;
}

//
// Scalar tensor division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator/(S const & s, Tensor<T, N> const & A)
{
  Tensor<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  split(A, s, B);

  return B;
}

//
// Tensor vector product v = A u
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator*(Tensor<T, N> const & A, Vector<S, N> const & u)
{
  return dot(A, u);
}

//
// Vector tensor product v = u A
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator*(Vector<S, N> const & u, Tensor<T, N> const & A)
{
  return dot(u, A);
}

//
// Tensor dot product C = A B
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
operator*(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  return dot(A, B);
}

//
// Tensor vector product v = A u
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
dot(Tensor<T, N> const & A, Vector<S, N> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

   Vector<typename Promote<S, T>::type, N>
  v(dimension);

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(i, p) * u(p);
      }
      v(i) = s;
    }
    break;

  case 3:
    v(0) = A(0, 0) * u(0) + A(0, 1) * u(1) + A(0, 2) * u(2);
    v(1) = A(1, 0) * u(0) + A(1, 1) * u(1) + A(1, 2) * u(2);
    v(2) = A(2, 0) * u(0) + A(2, 1) * u(1) + A(2, 2) * u(2);
    break;

  case 2:
    v(0) = A(0, 0) * u(0) + A(0, 1) * u(1);
    v(1) = A(1, 0) * u(0) + A(1, 1) * u(1);
    break;

  }

  return v;
}

//
// Vector tensor product v = u A
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
dot(Vector<S, N> const & u, Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

   Vector<typename Promote<S, T>::type, N>
  v(dimension);

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(p, i) * u(p);
      }
      v(i) = s;
    }
    break;

  case 3:
    v(0) = A(0, 0) * u(0) + A(1, 0) * u(1) + A(2, 0) * u(2);
    v(1) = A(0, 1) * u(0) + A(1, 1) * u(1) + A(2, 1) * u(2);
    v(2) = A(0, 2) * u(0) + A(1, 2) * u(1) + A(2, 2) * u(2);
    break;

  case 2:
    v(0) = A(0, 0) * u(0) + A(1, 0) * u(1);
    v(1) = A(0, 1) * u(0) + A(1, 1) * u(1);
    break;

  }

  return v;
}

//
// Tensor tensor product C = A B
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  C(dimension);

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i, p) * B(p, j);
        }
        C(i, j) = s;
      }
    }
    break;

  case 3:
    C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
    C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1);
    C(0, 2) = A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2);

    C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0);
    C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1);
    C(1, 2) = A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2);

    C(2, 0) = A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0);
    C(2, 1) = A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1);
    C(2, 2) = A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2);
    break;

  case 2:
    C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
    C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);

    C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
    C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);
    break;

  }

  return C;
}

//
// Tensor tensor product C = A^T B
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
t_dot(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  C(dimension);

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(p, i) * B(p, j);
        }
        C(i, j) = s;
      }
    }
    break;

  case 3:
    C(0, 0) = A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0);
    C(0, 1) = A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1);
    C(0, 2) = A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2);

    C(1, 0) = A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0);
    C(1, 1) = A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1);
    C(1, 2) = A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2);

    C(2, 0) = A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0);
    C(2, 1) = A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1);
    C(2, 2) = A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2);
    break;

  case 2:
    C(0, 0) = A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0);
    C(0, 1) = A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1);

    C(1, 0) = A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0);
    C(1, 1) = A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1);
    break;

  }

  return C;
}

//
// Tensor matrix product C = A^T B
//
template<typename S, typename T, Index N, Index P>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, N, P>
t_dot(Tensor<S, N> const & A, Matrix<T, N, P> const & B)
{
  Index const
  dimension{A.get_dimension()};

  Index const
  num_cols{B.get_num_cols()};

  assert(B.get_num_rows() == dimension);

  Matrix<typename Promote<S, T>::type, N, P>
  C(dimension, num_cols);

  for (Index i{0}; i < dimension; ++i) {
    for (Index j{0}; j < num_cols; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p{0}; p < dimension; ++p) {
        s += A(p, i) * B(p, j);
      }
      C(i, j) = s;
    }
  }

  return C;
}

//
// Tensor vector product C = A^T B
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
t_dot(Tensor<S, N> const & A, Vector<T, N> const & B)
{
  return dot(B, A);
}

//
// Tensor tensor product C = A B^T
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot_t(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  C(dimension);

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i, p) * B(j, p);
        }
        C(i, j) = s;
      }
    }
    break;

  case 3:
    C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2);
    C(0, 1) = A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2);
    C(0, 2) = A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2);

    C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2);
    C(1, 1) = A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2);
    C(1, 2) = A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2);

    C(2, 0) = A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2);
    C(2, 1) = A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2);
    C(2, 2) = A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2);
    break;

  case 2:
    C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1);
    C(0, 1) = A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1);

    C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1);
    C(1, 1) = A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1);
    break;

  }

  return C;
}

//
// Tensor tensor product C = A^T B^T
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
t_dot_t(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  C(dimension);

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(p, i) * B(j, p);
        }
        C(i, j) = s;
      }
    }
    break;

  case 3:
    C(0, 0) = A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1) + A(2, 0) * B(0, 2);
    C(0, 1) = A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1) + A(2, 0) * B(1, 2);
    C(0, 2) = A(0, 0) * B(2, 0) + A(1, 0) * B(2, 1) + A(2, 0) * B(2, 2);

    C(1, 0) = A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1) + A(2, 1) * B(0, 2);
    C(1, 1) = A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1) + A(2, 1) * B(1, 2);
    C(1, 2) = A(0, 1) * B(2, 0) + A(1, 1) * B(2, 1) + A(2, 1) * B(2, 2);

    C(2, 0) = A(0, 2) * B(0, 0) + A(1, 2) * B(0, 1) + A(2, 2) * B(0, 2);
    C(2, 1) = A(0, 2) * B(1, 0) + A(1, 2) * B(1, 1) + A(2, 2) * B(1, 2);
    C(2, 2) = A(0, 2) * B(2, 0) + A(1, 2) * B(2, 1) + A(2, 2) * B(2, 2);
    break;

  case 2:
    C(0, 0) = A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1);
    C(0, 1) = A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1);

    C(1, 0) = A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1);
    C(1, 1) = A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1);
    break;

  }

  return C;
}

//
// Tensor tensor double dot product (contraction)
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dotdot(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  typename Promote<S, T>::type
  s = 0.0;

  switch (dimension) {

  default:
    for (Index p = 0; p < dimension; ++p) {
      for (Index q = 0; q < dimension; ++q) {
        s += A(p, q) * B(p, q);
      }
    }
    break;

  case 3:
    s += A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2);
    s += A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2);
    s += A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2);
    break;

  case 2:
    s += A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1);
    s += A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1);
    break;

  }

  return s;
}

//
// dyad
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dyad(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Index const
  dimension = u.get_dimension();

  assert(v.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  A(dimension);

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {

      typename Promote<S, T>::type const
      s = u(i);

      for (Index j = 0; j < dimension; ++j) {
        A(i, j) = s * v(j);
      }
    }
    break;

  case 3:
    A(0, 0) = u(0) * v(0);
    A(0, 1) = u(0) * v(1);
    A(0, 2) = u(0) * v(2);

    A(1, 0) = u(1) * v(0);
    A(1, 1) = u(1) * v(1);
    A(1, 2) = u(1) * v(2);

    A(2, 0) = u(2) * v(0);
    A(2, 1) = u(2) * v(1);
    A(2, 2) = u(2) * v(2);
    break;

  case 2:
    A(0, 0) = u(0) * v(0);
    A(0, 1) = u(0) * v(1);

    A(1, 0) = u(1) * v(0);
    A(1, 1) = u(1) * v(1);
    break;

  }

  return A;
}

//
// bun operator, just for Jay, and now Reese too.
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
bun(Vector<S, N> const & u, Vector<T, N> const & v)
{
  return dyad(u, v);
}

//
// tensor product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
tensor(Vector<S, N> const & u, Vector<T, N> const & v)
{
  return dyad(u, v);
}

//
// diagonal tensor from vector
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
diag(Vector<T, N> const & v)
{
  Index const
  dimension = v.get_dimension();

  Tensor<T, N>
  A = zero<T, N>(dimension);

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      A(i, i) = v(i);
    }
    break;

  case 3:
    A(0, 0) = v(0);
    A(1, 1) = v(1);
    A(2, 2) = v(2);
    break;

  case 2:
    A(0, 0) = v(0);
    A(1, 1) = v(1);
    break;

  }

  return A;
}

//
// diagonal of tensor in a vector
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
diag(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Vector<T, N>
  v(dimension);

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      v(i) = A(i, i);
    }
    break;

  case 3:
    v(0) = A(0, 0);
    v(1) = A(1, 1);
    v(2) = A(2, 2);
    break;

  case 2:
    v(0) = A(0, 0);
    v(1) = A(1, 1);
    break;

  }

  return v;
}

//
// Zero 2nd-order tensor
// All components are zero
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
zero()
{
  return Tensor<T, N>(N, Filler::ZEROS);
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
zero(Index const dimension)
{
  return Tensor<T, DYNAMIC>(dimension, Filler::ZEROS);
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
zero(Index const dimension)
{
  return Tensor<T, N>(dimension, Filler::ZEROS);
}

// Local utility functions
namespace {

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void ones_in_diagonal(Tensor<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      A(i, i) = 1.0;
    }
    break;

  case 3:
    A(0, 0) = 1.0;
    A(1, 1) = 1.0;
    A(2, 2) = 1.0;
    break;

  case 2:
    A(0, 0) = 1.0;
    A(1, 1) = 1.0;
    break;

  }

  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void fill_levi_civita(Tensor<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      A(i, j) = levi_civita<T>(i, j);
    }
  }

  return;
}

} // anonymous namespace

//
// 2nd-order identity tensor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
identity()
{
  Tensor<T, N>
  A(N, Filler::ZEROS);

  ones_in_diagonal(A);

  return A;
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
identity(Index const dimension)
{
  Tensor<T, DYNAMIC>
  A(dimension, Filler::ZEROS);

  ones_in_diagonal(A);

  return A;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
identity(Index const dimension)
{

  Tensor<T, N>
  A(dimension, Filler::ZEROS);

  ones_in_diagonal(A);

  return A;
}

//
// 2nd-order identity tensor, à la Matlab
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
eye()
{
  return identity<T, N>();
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
eye(Index const dimension)
{
  return identity<T, DYNAMIC>(dimension);
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
eye(Index const dimension)
{
  return identity<T, N>(dimension);
}

//
// Levi-Civita symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
levi_civita_2()
{
  Tensor<T, N>
  A(N, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
levi_civita_2(Index const dimension)
{
  Tensor<T, DYNAMIC>
  A(dimension, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
levi_civita_2(Index const dimension)
{

  Tensor<T, DYNAMIC>
  A(dimension, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

//
// Permutation symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
permutation_2()
{
  return levi_civita_2<T, N>();
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
permutation_2(Index const dimension)
{
  return levi_civita_2<T, DYNAMIC>(dimension);
}


template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
permutation_2(Index const dimension)
{
  return levi_civita_2<T, N>(dimension);
}

//
// Alternating symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
alternator_2()
{
  return levi_civita_2<T, N>();
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC> const
alternator_2(Index const dimension)
{
  return levi_civita_2<T, DYNAMIC>(dimension);
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N> const
alternator_2(Index const dimension)
{
  return levi_civita_2<T, N>(dimension);
}

//
// 2nd-order tensor transpose
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
transpose(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  B = A;

  switch (dimension) {
  default:
    for (Index i = 0; i < dimension; ++i) {
      for (Index j = i + 1; j < dimension; ++j) {
        minitensor::swap(B(i, j), B(j, i));
      }
    }
    break;

  case 3:
    minitensor::swap(B(0, 1), B(1, 0));
    minitensor::swap(B(0, 2), B(2, 0));
    minitensor::swap(B(1, 2), B(2, 1));
    break;

  case 2:
    minitensor::swap(B(0, 1), B(1, 0));
    break;
  }

  return B;
}

//
// symmetric part of 2nd-order tensor
// \return \f$ \frac{1}{2}(A + A^T) \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
sym(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  B(dimension);

  switch (dimension) {

  default:
    B = 0.5 * (A + transpose(A));
    break;

  case 3:
    {
      T const & s00 = A(0, 0);
      T const & s11 = A(1, 1);
      T const & s22 = A(2, 2);

      T const s01 = 0.5 * (A(0, 1) + A(1, 0));
      T const s02 = 0.5 * (A(0, 2) + A(2, 0));
      T const s12 = 0.5 * (A(1, 2) + A(2, 1));

      B(0, 0) = s00;
      B(0, 1) = s01;
      B(0, 2) = s02;

      B(1, 0) = s01;
      B(1, 1) = s11;
      B(1, 2) = s12;

      B(2, 0) = s02;
      B(2, 1) = s12;
      B(2, 2) = s22;
    }
    break;

  case 2:
    {
      T const & s00 = A(0, 0);
      T const & s11 = A(1, 1);

      T const s01 = 0.5 * (A(0, 1) + A(1, 0));

      B(0, 0) = s00;
      B(0, 1) = s01;

      B(1, 0) = s01;
      B(1, 1) = s11;
    }
    break;

  }

  return B;
}

//
// skew symmetric part of 2nd-order tensor
// \return \f$ \frac{1}{2}(A - A^T) \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
skew(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N>
  B(dimension);

  switch (dimension) {

  default:
    B = 0.5 * (A - transpose(A));
    break;

  case 3:
    {
      T const s01 = 0.5 * (A(0, 1) - A(1, 0));
      T const s02 = 0.5 * (A(0, 2) - A(2, 0));
      T const s12 = 0.5 * (A(1, 2) - A(2, 1));

      B(0, 0) = 0.0;
      B(0, 1) = s01;
      B(0, 2) = s02;

      B(1, 0) = -s01;
      B(1, 1) = 0.0;
      B(1, 2) = s12;

      B(2, 0) = -s02;
      B(2, 1) = -s12;
      B(2, 2) = 0.0;
    }
    break;

  case 2:
    {
      T const s01 = 0.5 * (A(0, 1) - A(1, 0));

      B(0, 0) = 0.0;
      B(0, 1) = s01;

      B(1, 0) = -s01;
      B(1, 1) = 0.0;
    }
    break;

  }

  return B;
}

//
// skew symmetric 2nd-order tensor from vector, undefined
// for N!=3.
// \param u vector
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
skew(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  Tensor<T, N>
  A(dimension);

  switch (dimension) {

  case 3:
    A(0, 0) = 0.0;
    A(0, 1) = -u(2);
    A(0, 2) = u(1);

    A(1, 0) = u(2);
    A(1, 1) = 0.0;
    A(1, 2) = -u(0);

    A(2, 0) = -u(1);
    A(2, 1) = u(0);
    A(2, 2) = 0.0;
    break;

  default:
    MT_ERROR_EXIT("Skew from vector defined for 3D only");
    break;
  }

  return A;
}

} // namespace minitensor
namespace minitensor {

//
// tensor input
//
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Tensor<T, N> & A)
{

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      is >> A(i,j);
    }
  }

  return is;
}

//
// tensor output
//
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  if (dimension == 0) {
    return os;
  }

  os << std::scientific << std::setprecision(17);

  for (Index i = 0; i < dimension; ++i) {

    os << std::setw(24) << A(i,0);

    for (Index j = 1; j < dimension; ++j) {
      os << "," << std::setw(24) << A(i,j);
    }

    os << std::endl;
  }

  return os;
}

} // namespace minitensor

#endif //MiniTensor_Tensor_h
