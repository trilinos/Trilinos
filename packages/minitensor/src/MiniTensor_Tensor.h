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

#include "MiniTensor_Tensor.i.h"
#include "MiniTensor_Tensor.t.h"

#endif //MiniTensor_Tensor_h
