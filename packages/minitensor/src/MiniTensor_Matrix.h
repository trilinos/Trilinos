// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Matrix_h)
#define MiniTensor_Matrix_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include <Kokkos_Core.hpp>
#include "MiniTensor_Vector.h"

namespace minitensor {

template<typename T, Index M, Index N>
using matrix_store = Storage<T, dimension_product<M, N>::value>;

///
/// General rectangular matrix.
///
template<typename T, Index M = DYNAMIC, Index N = DYNAMIC>
class Matrix: public TensorBase<T, matrix_store<T, M, N>>
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
  using Store = matrix_store<T, M, N>;

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
  Matrix(Filler const value);

  KOKKOS_INLINE_FUNCTION
  explicit
  Matrix(Index const rows, Index cols, Filler const value);

  ///
  /// Create tensor from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      ArrayT & data,
      Index index1);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
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
  Matrix(
      Source const source,
      Index const rows,
      Index const cols,
      ArrayT & data,
      Index index1);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      Index const rows,
      Index const cols,
      ArrayT & data,
      Index index1,
      Index index2);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      Index const rows,
      Index const cols,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      Index const rows,
      Index const cols,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      Index const rows,
      Index const cols,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5);

  template<class ArrayT>
  KOKKOS_INLINE_FUNCTION
  Matrix(
      Source const source,
      Index const rows,
      Index const cols,
      ArrayT & data,
      Index index1,
      Index index2,
      Index index3,
      Index index4,
      Index index5,
      Index index6);

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
  Matrix(Matrix<T, M, N> const & A);

  ///
  /// Simple destructor
  ///
  virtual
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
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
operator+(Matrix<S, M, N> const & A, Matrix<T, M, N> const & B);

///
/// Matrix subtraction
/// \return \f$ A - B \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
operator-(Matrix<S, M, N> const & A, Matrix<T, M, N> const & B);

///
/// Matrix minus
/// \return \f$ -A \f$
///
template<typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N>
operator-(Matrix<T, M, N> const & A);

///
/// Matrix equality
/// Tested by components
/// \return \f$ A \equiv B \f$
///
template<typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator==(Matrix<T, M, N> const & A, Matrix<T, M, N> const & B);

///
/// Matrix inequality
/// Tested by components
/// \return \f$ A \neq B \f$
///
template<typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Matrix<T, M, N> const & A, Matrix<T, M, N> const & B);

///
/// Matrix vector product v = A u
/// \param A matrix
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, M>
operator*(Matrix<T, M, N> const & A, Vector<S, N> const & u);

///
/// Vector matrix product v = u A
/// \param A matrix
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
operator*(Vector<S, M> const & u, Matrix<T, M, N> const & A);

///
/// Matrix tensor product C = A B
/// \param A matrix
/// \param B tensor
/// \return \f$ A B \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
operator*(Matrix<T, M, N> const & A, Tensor<S, N> const & B);

///
/// Tensor matrix product C = A B
/// \param A tensor
/// \param B matrix
/// \return \f$ A B = \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
operator*(Tensor<S, M> const & A, Matrix<T, M, N> const & B);

///
/// Matrix dot product C = A B
/// \return \f$ A \cdot B \f$
///
template<typename S, typename T, Index M, Index P, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
operator*(Matrix<S, M, P> const & A, Matrix<T, P, N> const & B);

///
/// Scalar matrix product
/// \param s scalar
/// \param A matrix
/// \return \f$ s A \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_matrix<Promote<S, T>, M, N>>::type
operator*(S const & s, Matrix<T, M, N> const & A);

///
/// Matrix scalar product
/// \param A matrix
/// \param s scalar
/// \return \f$ s A \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if<order_1234<S>, apply_matrix<Promote<S, T>, M, N>>::type
operator*(Matrix<T, M, N> const & A, S const & s);

///
/// Matrix scalar division
/// \param A matrix
/// \param s scalar
/// \return \f$ A / s \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
operator/(Matrix<T, M, N> const & A, S const & s);

///
/// Scalar matrix division
/// \param s scalar
/// \param A matrix that divides scalar with each component
/// \return \f$ s / A \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
operator/(S const & s, Matrix<T, M, N> const & A);

///
/// Matrix input
/// \param A matrix
/// \param is input stream
/// \return is input stream
///
template<typename T, Index M, Index N>
std::istream &
operator>>(std::istream & is, Matrix<T, M, N> & A);

///
/// Matrix output
/// \param A matrix
/// \param os output stream
/// \return os output stream
///
template<typename T, Index M, Index N>
std::ostream &
operator<<(std::ostream & os, Matrix<T, M, N> const & A);

///
/// Extract a row as a vector
/// \param A matrix
/// \param i index of row
/// \return \f$ v = A(i,:) \f$
///
template<typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
row(Matrix<T, M, N> const & A, Index const i);

///
/// Extract a column as a vector
/// \param A matrix
/// \param j index of column
/// \return \f$ v = A(:,j) \f$
///
template<typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, M>
col(Matrix<T, M, N> const & A, Index const j);

///
/// Matrix vector product v = A u
/// \param A matrix
/// \param u vector
/// \return \f$ A u \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, M>
dot(Matrix<T, M, N> const & A, Vector<S, N> const & u);

///
/// Vector matrix product v = u A
/// \param A matrix
/// \param u vector
/// \return \f$ u A = A^T u \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
dot(Vector<S, M> const & u, Matrix<T, M, N> const & A);

///
/// Matrix tensor product C = A B
/// \param A matrix
/// \param B tensor
/// \return \f$ A B \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
dot(Matrix<T, M, N> const & A, Tensor<S, N> const & B);

///
/// Tensor matrix product C = A B
/// \param A tensor
/// \param B matrix
/// \return \f$ A B = \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
dot(Tensor<S, M> const & A, Matrix<T, M, N> const & B);

///
/// Matrix matrix product C = A B
/// \param A matrix
/// \param B matrix
/// \return a matrix \f$ A \cdot B \f$
///
template<typename S, typename T, Index M, Index P, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
dot(Matrix<S, M, P> const & A, Matrix<T, P, N> const & B);

///
/// Matrix matrix product C = A^T B
/// \param A matrix
/// \param B matrix
/// \return a matrix \f$ A^T \cdot B \f$
///
template<typename S, typename T, Index M, Index P, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
t_dot(Matrix<S, P, M> const & A, Matrix<T, P, N> const & B);

///
/// Matrix matrix product C = A B^T
/// \param A matrix
/// \param B matrix
/// \return a matrix \f$ A \cdot B^T \f$
///
template<typename S, typename T, Index M, Index P, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
dot_t(Matrix<S, M, P> const & A, Matrix<T, N, P> const & B);

///
/// Matrix matrix product C = A^T B^T
/// \param A matrix
/// \param B matrix
/// \return a tensor \f$ A^T \cdot B^T \f$
///
template<typename S, typename T, Index M, Index P, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
t_dot_t(Matrix<S, P, M> const & A, Matrix<T, N, P> const & B);

///
/// Matrix matrix double dot product (contraction)
/// \param A matrix
/// \param B matrix
/// \return a scalar \f$ A : B \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dotdot(Matrix<S, M, N> const & A, Matrix<T, M, N> const & B);

///
/// exterior matrix product
/// \param u vector
/// \param v vector
/// \return \f$ u \otimes v \f$
///
template<typename S, typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N>
matrix(Vector<S, M> const & u, Vector<T, N> const & v);

///
/// matrix transpose
///
template<typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N>
transpose(Matrix<T, N, M> const & A);

///
/// matrix adjoint
///
template<typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N>
adjoint(Matrix<T, N, M> const & A);

} // namespace minitensor

#include "MiniTensor_Matrix.i.h"
#include "MiniTensor_Matrix.t.h"

#endif //MiniTensor_Matrix_h
