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

#if !defined(Intrepid2_MiniTensor_Matrix_h)
#define Intrepid2_MiniTensor_Matrix_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include <boost/tuple/tuple.hpp>
#include<Kokkos_Core.hpp>
#include "Intrepid2_MiniTensor_Vector.h"

namespace Intrepid2 {

template<typename T, Index M, Index N,  typename ES>
using matrix_store = Storage<T, dimension_product<M, N>::value, ES>;

///
/// General rectangular matrix.
///
template<typename T, Index M = DYNAMIC, Index N = DYNAMIC,  typename ES=NOKOKKOS>
class Matrix: public TensorBase<T, matrix_store<T, M, N, ES>>
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
  IS_DYNAMIC = (M == DYNAMIC || N == DYNAMIC);

  ///
  /// Storage type
  ///
  using Store = matrix_store<T, M, N, ES>;

  ///
  /// Matrix order
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
  /// \param rows Number of rows
  /// \param cols Number of columns
  ///
  KOKKOS_INLINE_FUNCTION
  explicit
  Matrix();

  KOKKOS_INLINE_FUNCTION
  explicit
  Matrix(Index const rows, Index const cols);

  ///
  /// Create tensor from a specified value
  /// \param rows Number of rows
  /// \param cols Number of columns
  /// \param value all components are set equal to this
  ///
  KOKKOS_INLINE_FUNCTION
  explicit
  Matrix(ComponentValue const value);

  KOKKOS_INLINE_FUNCTION
  explicit
  Matrix(Index const rows, Index cols, ComponentValue const value);

  ///
  /// Create tensor from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(ArrayT & data, iType index1, iType index2, iType index3, iType index4);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(Index const rows, Index const cols, ArrayT & data, iType index1);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Index const rows, Index const cols,
      typename Kokkos::Impl::enable_if<
      !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Index const rows, Index const cols,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Index const rows, Index const cols,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Index const rows, Index const cols,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Index const rows, Index const cols,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  KOKKOS_INLINE_FUNCTION
  explicit
  Matrix(T const * data_ptr);

  KOKKOS_INLINE_FUNCTION
  explicit
  Matrix(Index const rows, Index cols, T const * data_ptr);

  ///
  /// Copy constructor
  ///
  KOKKOS_INLINE_FUNCTION
  Matrix(Matrix<T, M, N, ES> const & A);

  ///
  /// Simple destructor
  ///
  KOKKOS_INLINE_FUNCTION
  ~Matrix();

  ///
  /// Indexing for constant tensor
  /// \param i index
  /// \param j index
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  operator()(Index const i, Index const j) const;

  ///
  /// Matrix indexing
  /// \param i index
  /// \param j index
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  operator()(Index const i, Index const j);

  ///
  /// \return dimensions
  ///
  std::pair<Index, Index>
  get_dimensions() const;

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
  set_dimensions(Index const rows, Index const cols);

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

private:

  Index
  rows_{M};

  Index
  cols_{N};

};

///
/// Matrix addition
/// \return \f$ A + B \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator+(Matrix<S, M, N, ES> const & A, Matrix<T, M, N, ES> const & B);

///
/// Matrix subtraction
/// \return \f$ A - B \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator-(Matrix<S, M, N, ES> const & A, Matrix<T, M, N, ES> const & B);

///
/// Matrix minus
/// \return \f$ -A \f$
///
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>
operator-(Matrix<T, M, N, ES> const & A);

///
/// Matrix equality
/// Tested by components
/// \return \f$ A \equiv B \f$
///
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator==(Matrix<T, M, N, ES> const & A, Matrix<T, M, N, ES> const & B);

///
/// Matrix inequality
/// Tested by components
/// \return \f$ A \neq B \f$
///
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Matrix<T, M, N, ES> const & A, Matrix<T, M, N, ES> const & B);

///
/// Matrix vector product v = A u
/// \param A matrix
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, M, ES>
operator*(Matrix<T, M, N, ES> const & A, Vector<S, N, ES> const & u);

///
/// Vector matrix product v = u A
/// \param A matrix
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator*(Vector<S, M, ES> const & u, Matrix<T, M, N, ES> const & A);

///
/// Matrix dot product C = A B
/// \return \f$ A \cdot B \f$
///
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator*(Matrix<S, M, P, ES> const & A, Matrix<T, P, N, ES> const & B);

///
/// Scalar matrix product
/// \param s scalar
/// \param A matrix
/// \return \f$ s A \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_matrix<Promote<S, T>, M, N,ES>>::type
operator*(S const & s, Matrix<T, M, N, ES> const & A);

///
/// Matrix scalar product
/// \param A matrix
/// \param s scalar
/// \return \f$ s A \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_matrix<Promote<S, T>, M, N, ES>>::type
operator*(Matrix<T, M, N, ES> const & A, S const & s);

///
/// Matrix scalar division
/// \param A matrix
/// \param s scalar
/// \return \f$ A / s \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator/(Matrix<T, M, N, ES> const & A, S const & s);

///
/// Scalar matrix division
/// \param s scalar
/// \param A matrix that divides scalar with each component
/// \return \f$ s / A \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator/(S const & s, Matrix<T, M, N, ES> const & A);

///
/// Matrix input
/// \param A matrix
/// \param is input stream
/// \return is input stream
///
template<typename T, Index M, Index N,  typename ES>
std::istream &
operator>>(std::istream & is, Matrix<T, M, N, ES> & A);

///
/// Matrix output
/// \param A matrix
/// \param os output stream
/// \return os output stream
///
template<typename T, Index M, Index N,  typename ES>
std::ostream &
operator<<(std::ostream & os, Matrix<T, M, N, ES> const & A);

///
/// Extract a row as a vector
/// \param A matrix
/// \param i index of row
/// \return \f$ v = A(i,:) \f$
///
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
row(Matrix<T, M, N, ES> const & A, Index const i);

///
/// Extract a column as a vector
/// \param A matrix
/// \param j index of column
/// \return \f$ v = A(:,j) \f$
///
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, M, ES>
col(Matrix<T, M, N, ES> const & A, Index const j);

///
/// Matrix vector product v = A u
/// \param A matrix
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, M, ES>
dot(Matrix<T, M, N, ES> const & A, Vector<S, N, ES> const & u);

///
/// Vector matrix product v = u A
/// \param A matrix
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
dot(Vector<S, M, ES> const & u, Matrix<T, M, N, ES> const & A);

///
/// Matrix matrix product C = A B
/// \param A matrix
/// \param B matrix
/// \return a matrix \f$ A \cdot B \f$
///
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
dot(Matrix<S, M, P, ES> const & A, Matrix<T, P, N, ES> const & B);

///
/// Matrix matrix product C = A^T B
/// \param A matrix
/// \param B matrix
/// \return a matrix \f$ A^T \cdot B \f$
///
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
t_dot(Matrix<S, P, M, ES> const & A, Matrix<T, P, N, ES> const & B);

///
/// Matrix matrix product C = A B^T
/// \param A matrix
/// \param B matrix
/// \return a matrix \f$ A \cdot B^T \f$
///
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
dot_t(Matrix<S, M, P, ES> const & A, Matrix<T, N, P, ES> const & B);

///
/// Matrix matrix product C = A^T B^T
/// \param A matrix
/// \param B matrix
/// \return a tensor \f$ A^T \cdot B^T \f$
///
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
t_dot_t(Matrix<S, P, M, ES> const & A, Matrix<T, N, P, ES> const & B);

///
/// Matrix matrix double dot product (contraction)
/// \param A matrix
/// \param B matrix
/// \return a scalar \f$ A : B \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dotdot(Matrix<S, M, N, ES> const & A, Matrix<T, M, N, ES> const & B);

///
/// exterior matrix product
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
matrix(Vector<S, M, ES> const & u, Vector<T, N, ES> const & v);

///
/// matrix transpose
///
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, N, M, ES>
transpose(Matrix<T, M, N, ES> const & A);

///
/// matrix adjoint
///
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>
adjoint(Matrix<T, N, M, ES> const & A);

} // namespace Intrepid

#include "Intrepid2_MiniTensor_Matrix.i.h"
#include "Intrepid2_MiniTensor_Matrix.t.h"

#endif //Intrepid2_MiniTensor_Matrix_h
