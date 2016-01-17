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

#if !defined(Intrepid2_MiniTensor_Vector_h)
#define Intrepid2_MiniTensor_Vector_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include<Kokkos_Core.hpp>
#include "Intrepid2_MiniTensor_TensorBase.h"

namespace Intrepid2 {

template<typename T, Index N,  typename ES>
using vector_store = Storage<T, dimension_power<N, 1>::value, ES>;

///
/// Vector class.
///
template<typename T, Index N = DYNAMIC,  typename ES=NOKOKKOS>
class Vector: public TensorBase<T, vector_store<T, N, ES>>
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
  using Store = vector_store<T, N,ES>;

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
  Vector(ComponentValue const value);

  explicit
  KOKKOS_INLINE_FUNCTION
  Vector(Index const dimension, ComponentValue const value);

  ///
  /// Create vector from array.
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Index const dimension,
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Index const dimension,
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Vector(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>

  KOKKOS_INLINE_FUNCTION  
  Vector(
      Index const dimension,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  KOKKOS_INLINE_FUNCTION
  Vector(T const * data_ptr);

  KOKKOS_INLINE_FUNCTION
  Vector(Index const dimension, T const * data_ptr);

  ///
  /// Copy constructor
  ///
  KOKKOS_INLINE_FUNCTION
  Vector(Vector<T, N, ES> const & v);

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

};

///
/// Vector addition
/// \return \f$ u + v \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator+(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Vector substraction
/// \return \f$ u - v \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator-(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Vector minus
/// \return \f$ -u \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
operator-(Vector<T, N, ES> const & u);

///
/// Vector dot product
/// \return \f$ u \cdot v \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
operator*(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Vector equality tested by components
/// \return \f$ u \equiv v \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator==(Vector<T, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Vector inequality tested by components
/// \return \f$ u \neq v \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Vector<T, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Scalar vector product
/// \param s scalar factor
/// \param u vector factor
/// \return \f$ s u \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_vector<Promote<S, T>, N, ES> >::type
operator*(S const & s, Vector<T, N, ES> const & u);

///
/// Vector scalar product
/// \param u vector factor
/// \param s scalar factor
/// \return \f$ s u \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_vector<Promote<S, T>, N, ES> >::type
operator*(Vector<T, N, ES> const & u, S const & s);

///
/// Vector scalar division
/// \param u vector
/// \param s scalar that divides each component of vector
/// \return \f$ u / s \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator/(Vector<T, N, ES> const & u, S const & s);

///
/// Scalar vector division
/// \param s scalar
/// \param u vector that divides scalar with each component
/// \return \f$ s / u \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator/(S const & s, Vector<T, N, ES> const & u);

///
/// Vector dot product
/// \return \f$ u \cdot v \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dot(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Cross product only valid for R^3.
/// R^N with N != 3 will produce an error.
/// \return \f$ u \times v \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
cross(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v);

///
/// Vector 2-norm
/// \return \f$ \sqrt{u \cdot u} \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
norm(Vector<T, N, ES> const & u);

///
/// Vector 2-norm square. Used for fast distance calculation.
/// \return \f$ u \cdot u \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
norm_square(Vector<T, N, ES> const & u);

///
/// Vector 1-norm
/// \return \f$ |u_0|+|u_1|+|u_2| \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
norm_1(Vector<T, N, ES> const & u);

///
/// Vector infinity-norm
/// \return \f$ \max(|u_0|,|u_1|,|u_2|) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
norm_infinity(Vector<T, N, ES> const & u);

///
/// \return u / |u|, fails for |u| = 0
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
unit(Vector<T, N, ES> const & u);

///
/// Compute Householder vector
/// See algorithm 5.1.1, Matrix Computations 3rd ed, Golub - Van Loan
/// \f$ Px = |X|e_1, P := I - \beta v v^T\f$
/// \return v, beta
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Vector<T, N, ES>, T>
house(Vector<T, N, ES> const & x);

///
/// Vector input
/// \param u vector
/// \param is input stream
/// \return is input stream
///
template<typename T, Index N,  typename ES>
std::istream &
operator>>(std::istream & is, Vector<T, N, ES> & u);

///
/// Vector output
/// \param u vector
/// \param os output stream
/// \return os output stream
///
template<typename T, Index N,  typename ES>
std::ostream &
operator<<(std::ostream & os, Vector<T, N, ES> const & u);

} // namespace Intrepid

#include "Intrepid2_MiniTensor_Vector.i.h"
#include "Intrepid2_MiniTensor_Vector.t.h"

#endif //Intrepid2_MiniTensor_Vector_h
