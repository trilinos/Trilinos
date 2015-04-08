// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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

#if !defined(Intrepid_MiniTensor_Tensor_h)
#define Intrepid_MiniTensor_Tensor_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include <boost/tuple/tuple.hpp>
#ifdef HAVE_INTREPID_KOKKOSCORE
#include<Kokkos_Core.hpp>
#endif
#include "Intrepid_MiniTensor_Vector.h"

namespace Intrepid {

///
/// Component ordering convention
///
enum ComponentOrder {
  CANONICAL, SIERRA_FULL, SIERRA_SYMMETRIC
};

template<typename T, Index N>
struct tensor_store
{
  typedef Storage<T, dimension_power<N, 2>::value> type;
};

///
/// Second order tensor.
///
template<typename T, Index N = DYNAMIC>
class Tensor: public TensorBase<T, typename tensor_store<T, N>::type>
{
public:

  ///
  /// Order
  ///
  static
  Index const
  ORDER = 2;

  ///
  /// Static or dynamic
  ///
  static
  bool const
  IS_DYNAMIC = N == DYNAMIC;

  ///
  /// Storage type
  ///
  typedef typename tensor_store<T, N>::type
  Store;

  ///
  /// Tensor order
  ///
  static
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
  Tensor();

  explicit
  Tensor(Index const dimension);

  ///
  /// Create tensor from a specified value
  /// \param dimension the space dimension
  /// \param value all components are set equal to this
  ///
  explicit
  Tensor(ComponentValue const value);

  explicit
  Tensor(Index const dimension, ComponentValue const value);
  ///
  /// Create tensor from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
#ifdef HAVE_INTREPID_KOKKOSCORE
  template<class ArrayT, typename iType>
  Tensor(ArrayT & data, iType index1);

  template<class ArrayT, typename iType>
  Tensor(
      typename apply_diff<ArrayT, Index>::type & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  Tensor(
      typename apply_diff<ArrayT, Index>::type & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  Tensor(ArrayT & data, iType index1, iType index2, iType index3, iType index4);

  template<class ArrayT, typename iType>
  Tensor(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  Tensor(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  template<class ArrayT, typename iType>
  Tensor(Index const dimension, ArrayT & data, iType index1);

  template<class ArrayT, typename iType>
  Tensor(
      Index const dimension,
      typename apply_diff<ArrayT, Index>::type & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  Tensor(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  Tensor(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4);

  template<class ArrayT, typename iType>
  Tensor(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  Tensor(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);
  #endif

  explicit
  Tensor(T const * data_ptr);

  explicit
  Tensor(Index const dimension, T const * data_ptr);
  ///
  /// Copy constructor
  ///
  Tensor(Tensor<T, N> const & A);

  ///
  /// 2nd-order tensor from 4th-order tensor
  ///
  Tensor(Tensor4<T, dimension_sqrt<N>::value> const & A);

  ///
  /// Create tensor specifying components
  /// \param  s00 s01 ... components in the R^2 canonical basis
  ///
  Tensor(T const & s00, T const & s01, T const & s10, T const & s11);

  ///
  /// Create tensor specifying components
  /// \param  s00 s01 ... components in the R^3 canonical basis
  ///
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
  Tensor(T const * data_ptr, ComponentOrder const component_order);

  explicit
  Tensor(Index const dimension, T const * data_ptr,
      ComponentOrder const component_order);
  ///
  /// Simple destructor
  ///
  ~Tensor();

  ///
  /// Indexing for constant tensor
  /// \param i index
  /// \param j index
  ///
  T const &
  operator()(Index const i, Index const j) const;

  ///
  /// Tensor indexing
  /// \param i index
  /// \param j index
  ///
  T &
  operator()(Index const i, Index const j);

  ///
  /// \return dimension
  ///
  Index
  get_dimension() const;

  ///
  /// \param dimension of vector
  ///
  void
  set_dimension(Index const dimension);

  ///
  /// Fill components with value specification
  /// \param value all components are set equal to this specification
  ///
  void
  fill(ComponentValue const value);

  ///
  /// Fill components with value as parameter
  /// \param value all components are set equal to this parameter
  ///
  void
  fill(T const & s);

  ///
  /// Fill components from array defined by pointer.
  /// \param data_ptr pointer into array for filling components
  ///
  template<class ArrayT, typename iType>
  typename if_diff<ArrayT, T*, void>::type
  fill(ArrayT & data, iType index1);

  template<class ArrayT, typename iType>
  void
  fill(ArrayT & data, iType index1, iType index2);

  template<class ArrayT, typename iType1, typename iType2, typename iType3>
  void
  fill(ArrayT & data, iType1 index1, iType2 index2, iType3 index3);

  template<class ArrayT, typename iType>
  void
  fill(ArrayT & data, iType index1, iType index2, iType index3, iType index4);

  template<class ArrayT, typename iType>
  void
  fill(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  void
  fill(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  void
  fill(T const * data_ptr);

  ///
  /// Fill components from array defined by pointer.
  /// \param data_ptr pointer into array for filling components
  /// \param component_order component convention (3D only)
  ///
  void
  fill(T const * data_ptr, ComponentOrder const component_order);

};

///
/// Tensor addition
/// \return \f$ A + B \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
operator+(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor subtraction
/// \return \f$ A - B \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
operator-(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor minus
/// \return \f$ -A \f$
///
template<typename T, Index N>
Tensor<T, N>
operator-(Tensor<T, N> const & A);

///
/// Tensor equality
/// Tested by components
/// \return \f$ A \equiv B \f$
///
template<typename T, Index N>
bool
operator==(Tensor<T, N> const & A, Tensor<T, N> const & B);

///
/// Tensor inequality
/// Tested by components
/// \return \f$ A \neq B \f$
///
template<typename T, Index N>
bool
operator!=(Tensor<T, N> const & A, Tensor<T, N> const & B);

///
/// Tensor vector product v = A u
/// \param A tensor
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index N>
Vector<typename Promote<S, T>::type, N>
operator*(Tensor<T, N> const & A, Vector<S, N> const & u);

///
/// Vector tensor product v = u A
/// \param A tensor
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index N>
Vector<typename Promote<S, T>::type, N>
operator*(Vector<S, N> const & u, Tensor<T, N> const & A);

///
/// Tensor dot product C = A B
/// \return \f$ A \cdot B \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
operator*(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Scalar tensor product
/// \param s scalar
/// \param A tensor
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N>
typename lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N> >::type
operator*(S const & s, Tensor<T, N> const & A);

///
/// Tensor scalar product
/// \param A tensor
/// \param s scalar
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N>
typename lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N> >::type
operator*(Tensor<T, N> const & A, S const & s);

///
/// Tensor scalar division
/// \param A tensor
/// \param s scalar
/// \return \f$ A / s \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
operator/(Tensor<T, N> const & A, S const & s);

///
/// Scalar tensor division
/// \param s scalar
/// \param A tensor that divides scalar with each component
/// \return \f$ s / A \f$
///
template<typename S, typename T, Index N>
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
Vector<T, N>
row(Tensor<T, N> const & A, Index const i);

///
/// Extract a column as a vector
/// \param A tensor
/// \param j index of column
/// \return \f$ v = A(:,j) \f$
///
template<typename T, Index N>
Vector<T, N>
col(Tensor<T, N> const & A, Index const j);

///
/// Tensor vector product v = A u
/// \param A tensor
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index N>
Vector<typename Promote<S, T>::type, N>
dot(Tensor<T, N> const & A, Vector<S, N> const & u);

///
/// Vector tensor product v = u A
/// \param A tensor
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index N>
Vector<typename Promote<S, T>::type, N>
dot(Vector<S, N> const & u, Tensor<T, N> const & A);

///
/// Tensor tensor product C = A B
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A \cdot B \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
dot(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor tensor product C = A^T B
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A^T \cdot B \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
t_dot(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor tensor product C = A B^T
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A \cdot B^T \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
dot_t(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor tensor product C = A^T B^T
/// \param A tensor
/// \param B tensor
/// \return a tensor \f$ A^T \cdot B^T \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
t_dot_t(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Tensor tensor double dot product (contraction)
/// \param A tensor
/// \param B tensor
/// \return a scalar \f$ A : B \f$
///
template<typename S, typename T, Index N>
typename Promote<S, T>::type
dotdot(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// Dyad
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
dyad(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Bun operator, just for Jay, and now Reese too.
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
bun(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Tensor product
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
tensor(Vector<S, N> const & u, Vector<T, N> const & v);

///
/// Diagonal tensor from vector
/// \param v vector
/// \return A = diag(v)
///
template<typename T, Index N>
Tensor<T, N>
diag(Vector<T, N> const & v);

///
/// Diagonal of tensor in a vector
/// \param A tensor
/// \return v = diag(A)
///
template<typename T, Index N>
Vector<T, N>
diag(Tensor<T, N> const & A);

///
/// Zero 2nd-order tensor
/// All components are zero
///
template<typename T, Index N>
Tensor<T, N> const
zero();

template<typename T>
Tensor<T, DYNAMIC> const
zero(Index const dimension);

template<typename T, Index N>
Tensor<T, N> const
zero(Index const dimension);

///
/// 2nd-order identity tensor
///
template<typename T, Index N>
Tensor<T, N> const
identity();

template<typename T>
Tensor<T, DYNAMIC> const
identity(Index const dimension);

template<typename T, Index N>
Tensor<T, N> const
identity(Index const dimension);

///
/// 2nd-order identity tensor, à la Matlab
///
template<typename T, Index N>
Tensor<T, N> const
eye();

template<typename T>
Tensor<T, DYNAMIC> const
eye(Index const dimension);

template<typename T, Index N>
Tensor<T, N> const
eye(Index const dimension);

///
/// Levi-Civita symbol
///
template<typename T, Index N>
Tensor<T, N> const
levi_civita_2();

template<typename T>
Tensor<T, DYNAMIC> const
levi_civita_2(Index const dimension);

template<typename T, Index N>
Tensor<T, N> const
levi_civita_2(Index const dimension);

///
/// Permutation symbol
///
template<typename T, Index N>
Tensor<T, N> const
permutation_2();

template<typename T>
Tensor<T, DYNAMIC> const
permutation_2(Index const dimension);

template<typename T, Index N>
Tensor<T, N> const
permutation_2(Index const dimension);

///
/// Alternating symbol
///
template<typename T, Index N>
Tensor<T, N> const
alternator_2();

template<typename T>
Tensor<T, DYNAMIC> const
alternator_2(Index const dimension);

template<typename T, Index N>
Tensor<T, N> const
alternator_2(Index const dimension);

///
/// R^N 2nd-order tensor transpose
///
template<typename T, Index N>
Tensor<T, N>
transpose(Tensor<T, N> const & A);

///
/// C^N 2nd-order tensor adjoint
///
template<typename T, Index N>
Tensor<T, N>
adjoint(Tensor<T, N> const & A);

///
/// Symmetric part of 2nd-order tensor
/// \return \f$ \frac{1}{2}(A + A^T) \f$
///
template<typename T, Index N>
Tensor<T, N>
sym(Tensor<T, N> const & A);

///
/// Skew symmetric part of 2nd-order tensor
/// \return \f$ \frac{1}{2}(A - A^T) \f$
///
template<typename T, Index N>
Tensor<T, N>
skew(Tensor<T, N> const & A);

///
/// Skew symmetric 2nd-order tensor from vector valid for R^3 only.
/// R^N with N != 3 will produce an error
/// \param u vector
/// \return \f$ {{0, -u_2, u_1}, {u_2, 0, -u_0}, {-u_1, u+0, 0}} \f$
///
template<typename T, Index N>
Tensor<T, N>
skew(Vector<T, N> const & u);

} // namespace Intrepid

#include "Intrepid_MiniTensor_Tensor.i.h"
#include "Intrepid_MiniTensor_Tensor.t.h"

#endif //Intrepid_MiniTensor_Tensor_h
