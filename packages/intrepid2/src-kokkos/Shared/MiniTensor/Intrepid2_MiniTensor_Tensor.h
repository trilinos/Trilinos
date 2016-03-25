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

#if !defined(Intrepid2_MiniTensor_Tensor_h)
#define Intrepid2_MiniTensor_Tensor_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include <boost/tuple/tuple.hpp>
#include<Kokkos_Core.hpp>
#include "Intrepid2_MiniTensor_Vector.h"

namespace Intrepid2 {

///
/// Component ordering convention
///
enum ComponentOrder {
  CANONICAL, SIERRA_FULL, SIERRA_SYMMETRIC
};

template<typename T, Index N,  typename ES>
using tensor_store = Storage<T, dimension_power<N, 2>::value, ES>;

///
/// Second order tensor.
///
template<typename T, Index N = DYNAMIC,  typename ES=NOKOKKOS>
class Tensor: public TensorBase<T, tensor_store<T, N,ES>>
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
  using Store = tensor_store<T, N, ES>;

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
  Tensor(ComponentValue const value);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor(Index const dimension, ComponentValue const value);

  ///
  /// Create tensor from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(ArrayT & data, iType index1);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(ArrayT & data, iType index1, iType index2, iType index3, iType index4);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(Index const dimension, ArrayT & data, iType index1);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Index const dimension,
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Tensor(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

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
  Tensor(Tensor<T, N, ES> const & A);

  ///
  /// 2nd-order tensor from 4th-order tensor
  ///
  KOKKOS_INLINE_FUNCTION
  Tensor(Tensor4<T, dimension_sqrt<N>::value, ES> const & A);

  ///
  /// Create tensor specifying components
  /// \param  s00 s01 ... components in the R^2 canonical basis
  ///
  //
  KOKKOS_INLINE_FUNCTION
  Tensor(T const & s00, T const & s01, T const & s10, T const & s11);

  ///
  /// Create tensor specifying components
  /// \param  s00 s01 ... components in the R^3 canonical basis
  ///
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
  Tensor(Index const dimension, T const * data_ptr,
      ComponentOrder const component_order);

  ///
  /// Simple destructor
  ///
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
  /// \param dimension of vector
  ///
  KOKKOS_INLINE_FUNCTION
  void
  set_dimension(Index const dimension);

  ///
  /// Fill components with value specification
  /// \param value all components are set equal to this specification
  ///
  KOKKOS_INLINE_FUNCTION
  void
  fill(ComponentValue const value);

  ///
  /// Fill components with value as parameter
  /// \param value all components are set equal to this parameter
  ///
  KOKKOS_INLINE_FUNCTION
  void
  fill(T const & s);

  ///
  /// Fill components from array defined by pointer.
  /// \param data_ptr pointer into array for filling components
  ///
  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  typename Kokkos::Impl::enable_if<
  !Kokkos::Impl::is_same<ArrayT, T*>::value, void>::type
  fill(ArrayT & data, iType index1);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void
  fill(ArrayT & data, iType index1, iType index2);

  template<class ArrayT, typename iType1, typename iType2, typename iType3>
  KOKKOS_INLINE_FUNCTION
  void
  fill(ArrayT & data, iType1 index1, iType2 index2, iType3 index3);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void
  fill(ArrayT & data, iType index1, iType index2, iType index3, iType index4);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void
  fill(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void
  fill(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  KOKKOS_INLINE_FUNCTION
  void
  fill(T const * data_ptr);

  ///
  /// Fill components from array defined by pointer.
  /// \param data_ptr pointer into array for filling components
  /// \param component_order component convention (3D only)
  ///
  KOKKOS_INLINE_FUNCTION
  void
  fill(T const * data_ptr, ComponentOrder const component_order);

};

///
/// Tensor addition
/// \return \f$ A + B \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator+(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Tensor subtraction
/// \return \f$ A - B \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator-(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Tensor minus
/// \return \f$ -A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
operator-(Tensor<T, N, ES> const & A);

///
/// Tensor equality
/// Tested by components
/// \return \f$ A \equiv B \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor<T, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Tensor inequality
/// Tested by components
/// \return \f$ A \neq B \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor<T, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Tensor vector product v = A u
/// \param A tensor
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator*(Tensor<T, N, ES> const & A, Vector<S, N, ES> const & u);

///
/// Vector tensor product v = u A
/// \param A tensor
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator*(Vector<S, N, ES> const & u, Tensor<T, N, ES> const & A);

///
/// Tensor dot product C = A B
/// \return \f$ A \cdot B \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator*(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Scalar tensor product
/// \param s scalar
/// \param A tensor
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N, class  ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N, ES> >::type
operator*(S const & s, Tensor<T, N, ES> const & A);

///
/// Tensor scalar product
/// \param A tensor
/// \param s scalar
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N, class  ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N, ES> >::type
operator*(Tensor<T, N, ES> const & A, S const & s);

///
/// Tensor scalar division
/// \param A tensor
/// \param s scalar
/// \return \f$ A / s \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator/(Tensor<T, N, ES> const & A, S const & s);

///
/// Scalar tensor division
/// \param s scalar
/// \param A tensor that divides scalar with each component
/// \return \f$ s / A \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator/(S const & s, Tensor<T, N, ES> const & A);

///
/// Tensor input
/// \param A tensor
/// \param is input stream
/// \return is input stream
///
template<typename T, Index N, class  ES>
std::istream &
operator>>(std::istream & is, Tensor<T, N, ES> & A);

///
/// Tensor output
/// \param A tensor
/// \param os output stream
/// \return os output stream
///
template<typename T, Index N,  typename ES>
std::ostream &
operator<<(std::ostream & os, Tensor<T, N, ES> const & A);

///
/// Extract a row as a vector
/// \param A tensor
/// \param i index of row
/// \return \f$ v = A(i,:) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
row(Tensor<T, N, ES> const & A, Index const i);

///
/// Extract a column as a vector
/// \param A tensor
/// \param j index of column
/// \return \f$ v = A(:,j) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
col(Tensor<T, N, ES> const & A, Index const j);

///
/// Tensor vector product v = A u
/// \param A tensor
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
dot(Tensor<T, N, ES> const & A, Vector<S, N, ES> const & u);

///
/// Vector tensor product v = u A
/// \param A tensor
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index N, class  ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
dot(Vector<S, N, ES> const & u, Tensor<T, N, ES> const & A);

///
/// Tensor tensor product C = A B
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A \cdot B \f$
///
template<typename S, typename T, Index N, class  ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dot(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Tensor tensor product C = A^T B
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A^T \cdot B \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
t_dot(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Tensor tensor product C = A B^T
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A \cdot B^T \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dot_t(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Tensor tensor product C = A^T B^T
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A^T \cdot B^T \f$
///
template<typename S, typename T, Index N, class  ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
t_dot_t(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Tensor tensor double dot product (contraction)
/// \param A tensor
/// \param B tensor
/// \return a scalar \f$ A : B \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dotdot(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B);

///
/// Dyad
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index N, class  ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dyad(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Bun operator, just for Jay, and now Reese too.
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
bun(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Tensor product
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index N, class  ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
tensor(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Diagonal tensor from vector
/// \param v vector
/// \return A = diag(v)
///
template<typename T, Index N, class  ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
diag(Vector<T, N, ES> const & v);

///
/// Diagonal of tensor in a vector
/// \param A tensor
/// \return v = diag(A)
///
template<typename T, Index N, class  ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
diag(Tensor<T, N, ES> const & A);

///
/// Zero 2nd-order tensor
/// All components are zero
///
template<typename T, Index N, class  ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION 
Tensor<T, N, ES> const
zero();

template<typename T,  typename ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
zero(Index const dimension);

template<typename T, Index N, class  ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
zero(Index const dimension);

///
/// 2nd-order identity tensor
///
template<typename T, Index N, class  ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
identity();

template<typename T, class  ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
identity(Index const dimension);

template<typename T, Index N,  typename ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
identity(Index const dimension);

///
/// 2nd-order identity tensor, à la Matlab
///
template<typename T, Index N, class  ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
eye();

template<typename T,  typename ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
eye(Index const dimension);

template<typename T, Index N,  typename ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
eye(Index const dimension);

///
/// Levi-Civita symbol
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
levi_civita_2();

template<typename T, class  ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
levi_civita_2(Index const dimension);

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
levi_civita_2(Index const dimension);

///
/// Permutation symbol
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
permutation_2();

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
permutation_2(Index const dimension);

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
permutation_2(Index const dimension);

///
/// Alternating symbol
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
alternator_2();

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
alternator_2(Index const dimension);

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
alternator_2(Index const dimension);

///
/// 2nd-order tensor transpose
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
transpose(Tensor<T, N, ES> const & A);

///
/// C^N 2nd-order tensor adjoint
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
adjoint(Tensor<T, N, ES> const & A);

///
/// Symmetric part of 2nd-order tensor
/// \return \f$ \frac{1}{2}(A + A^T) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
sym(Tensor<T, N, ES> const & A);

///
/// Skew symmetric part of 2nd-order tensor
/// \return \f$ \frac{1}{2}(A - A^T) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
skew(Tensor<T, N, ES> const & A);

///
/// Skew symmetric 2nd-order tensor from vector valid for R^3 only.
/// R^N with N != 3 will produce an error
/// \param u vector
/// \return \f$ {{0, -u_2, u_1}, {u_2, 0, -u_0}, {-u_1, u+0, 0}} \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
skew(Vector<T, N, ES> const & u);

} // namespace Intrepid

#include "Intrepid2_MiniTensor_Tensor.i.h"
#include "Intrepid2_MiniTensor_Tensor.t.h"

#endif //Intrepid2_MiniTensor_Tensor_h
