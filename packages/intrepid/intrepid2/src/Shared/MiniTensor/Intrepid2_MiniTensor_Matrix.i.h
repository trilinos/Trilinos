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

#if !defined(Intrepid2_MiniTensor_Matrix_i_h)
#define Intrepid2_MiniTensor_Matrix_i_h

namespace Intrepid2 {

//
// Constructor that initializes to NaNs
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix() :
    TensorBase<T, Store>::TensorBase()
{
  set_dimensions(M, N);
  return;
}

template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(Index const rows, Index const cols) :
    TensorBase<T, Store>::TensorBase(rows * cols, ORDER),
    rows_(rows), cols_(cols)
{
  return;
}

///
/// Create matrix from a specified value
///
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N,ES>::Matrix(ComponentValue const value) :
    TensorBase<T, Store>::TensorBase(M * N, ORDER, value)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    Index const rows,
    Index const cols ,
    ComponentValue const value) :
    TensorBase<T, Store>::TensorBase(rows * cols, ORDER, value),
    rows_(rows), cols_(cols)
{
  return;
}

//
//  Create tensor from array
//
//
template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1) :
    TensorBase<T, Store>::TensorBase(M * N, ORDER, data, index1)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1,
    iType index2) :
    TensorBase<T, Store>::TensorBase(M * N, ORDER, data, index1, index2)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1,
    iType index2,
    iType index3) :
    TensorBase<T, Store>::TensorBase(M * N, ORDER, data, index1, index2, index3)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4) :
    TensorBase<T, Store>::TensorBase(
        M * N,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5) :
    TensorBase<T, Store>::TensorBase(
        M * N,
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

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5,
    iType index6) :
    TensorBase<T, Store>::TensorBase(
        M * N,
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

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    Index const rows,
    Index const cols,
    ArrayT & data,
    iType index1) :
    TensorBase<T, Store>::TensorBase(rows * cols, ORDER, data, index1),
    rows_(rows), cols_(cols)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    Index const rows,
    Index const cols,
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1,
    iType index2) :
    TensorBase<T, Store>::TensorBase(rows * cols, ORDER, data, index1, index2),
    rows_(rows), cols_(cols)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    Index const rows,
    Index const cols,
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3) :
    TensorBase<T, Store>::TensorBase(
        rows * cols,
        ORDER,
        data,
        index1,
        index2,
        index3),
        rows_(rows), cols_(cols)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    Index const rows,
    Index const cols,
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4) :
    TensorBase<T, Store>::TensorBase(
        rows * cols,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4),
        rows_(rows), cols_(cols)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    Index const rows,
    Index const cols,
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5) :
    TensorBase<T, Store>::TensorBase(
        rows * cols,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4,
        index5),
        rows_(rows), cols_(cols)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    Index const rows,
    Index const cols,
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5,
    iType index6) :
    TensorBase<T, Store>::TensorBase(
        rows * cols,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4,
        index5,
        index6),
        rows_(rows), cols_(cols)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N,ES>::Matrix(T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(M * N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(
    Index const rows,
    Index const cols,
    T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(rows * cols, ORDER, data_ptr),
    rows_(rows), cols_(cols)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::Matrix(Matrix<T, M, N, ES> const & A) :
    TensorBase<T, Store>::TensorBase(A),
    rows_(A.get_num_rows()), cols_(A.get_num_cols())
{
  return;
}

//
// Simple destructor
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>::~Matrix()
{
  return;
}

//
// Get dimensions
//
template<typename T, Index M, Index N,  typename ES>
inline
std::pair<Index, Index>
Matrix<T, M, N, ES>::get_dimensions() const
{
  return std::make_pair(rows_, cols_);
}

//
// Get number rows
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Index
Matrix<T, M, N, ES>::get_num_rows() const
{
  return rows_;
}

//
// Get number cols
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Index
Matrix<T, M, N, ES>::get_num_cols() const
{
  return cols_;
}

//
// Set dimensions
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Matrix<T, M, N, ES>::set_dimensions(Index const rows, Index const cols)
{
  if (IS_DYNAMIC == false) {
    assert(rows * cols <= M * N);
  }

  TensorBase<T, Store>::set_dimension(rows * cols, ORDER);
  rows_ = rows;
  cols_ = cols;

  return;
}

//
// Indexing for constant matrix
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T const &
Matrix<T, M, N, ES>::operator()(Index const i, Index const j) const
    {
  Matrix<T, M, N, ES> const &
  self = (*this);

  Index const
  num_cols = self.get_num_cols();

  return self[i * num_cols + j];
}

//
// Matrix indexing
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T &
Matrix<T, M, N, ES>::operator()(Index const i, Index const j)
{
  Matrix<T, M, N, ES> &
  self = (*this);

  Index const
  num_cols = self.get_num_cols();

  return self[i * num_cols + j];
}

//
// Fill components with value specification
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Matrix<T, M, N, ES>::fill(ComponentValue const value)
{
  TensorBase<T, Store>::fill(value);
  return;
}

//
// Fill components with value as parameter
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Matrix<T, M, N, ES>::fill(T const & s)
{
  TensorBase<T, Store>::fill(s);
  return;
}

//
// Fill components from array defined by pointer.
//
template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
typename Kokkos::Impl::enable_if<
!Kokkos::Impl::is_same<ArrayT, T*>::value, void>::type
Matrix<T, M, N, ES>::fill(
    ArrayT & data,
    iType index1)
{
  TensorBase<T, Store>::fill(data, index1);
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
void
Matrix<T, M, N, ES>::fill(
    ArrayT & data,
    iType index1,
    iType index2)
{
  TensorBase<T, Store>::fill(data, index1, index2);
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType1, typename iType2, typename iType3>
KOKKOS_INLINE_FUNCTION
void
Matrix<T, M, N, ES>::fill(
    ArrayT & data,
    iType1 index1,
    iType2 index2,
    iType3 index3)
{
  TensorBase<T, Store>::fill(data, index1, index2, index3);
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
void
Matrix<T, M, N, ES>::fill(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4)
{
  TensorBase<T, Store>::fill(data, index1, index2, index3, index4);
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
void
Matrix<T, M, N, ES>::fill(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5)
{
  TensorBase<T, Store>::fill(data, index1, index2, index3, index4, index5);
  return;
}

template<typename T, Index M, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
void
Matrix<T, M, N, ES>::fill(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5,
    iType index6)
{
  TensorBase<T, Store>::fill(
      data,
      index1,
      index2,
      index3,
      index4,
      index5,
      index6);
  return;
}

template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Matrix<T, M, N, ES>::fill(T const * data_ptr)
{
  TensorBase<T, Store>::fill(data_ptr);
  return;
}

//
// Extract a row as a vector
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
row(Matrix<T, M, N, ES> const & A, Index const i)
{
  Index const
  num_cols = A.get_num_cols();

  Vector<T, N, ES>
  v(num_cols);

  switch (num_cols) {
  default:
    for (Index j = 0; j < num_cols; ++j) {
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
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, M, ES>
col(Matrix<T, M, N, ES> const & A, Index const j)
{
  Index const
  num_rows = A.get_num_rows();

  Vector<T, M, ES>
  v(num_rows);

  switch (num_rows) {
  default:
    for (Index i = 0; i < num_rows; ++i) {
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
// Matrix addition
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator+(Matrix<S, M, N, ES> const & A, Matrix<T, M, N, ES> const & B)
{
  Matrix<typename Promote<S, T>::type, M, N, ES>
  C(A.get_num_rows(), A.get_num_cols());

  add(A, B, C);

  return C;
}

//
// Matrix subtraction
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator-(Matrix<S, M, N, ES> const & A, Matrix<T, M, N, ES> const & B)
{
  Matrix<typename Promote<S, T>::type, M, N, ES>
  C(A.get_num_rows(), A.get_num_cols());

  subtract(A, B, C);

  return C;
}

//
// Matrix minus
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>
operator-(Matrix<T, M, N, ES> const & A)
{
  Matrix<T, M, N, ES>
  B(A.get_num_rows(), A.get_num_cols());

  minus(A, B);

  return B;
}

//
// Matrix equality
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator==(Matrix<T, M, N, ES> const & A, Matrix<T, M, N, ES> const & B)
{
  return equal(A, B);
}

//
// Matrix inequality
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Matrix<T, M, N, ES> const & A, Matrix<T, M, N, ES> const & B)
{
  return not_equal(A, B);
}

//
// Scalar matrix product
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename
lazy_disable_if<order_1234<S>, apply_matrix<Promote<S, T>, M, N, ES>>::type
operator*(S const & s, Matrix<T, M, N, ES> const & A)
{
  Matrix<typename Promote<S, T>::type, M, N, ES>
  B(A.get_num_rows(), A.get_num_cols());

  scale(A, s, B);

  return B;
}

//
// Matrix scalar product
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename
lazy_disable_if<order_1234<S>, apply_matrix<Promote<S, T>, M, N, ES>>::type
operator*(Matrix<T, M, N, ES> const & A, S const & s)
{
  Matrix<typename Promote<S, T>::type, M, N, ES>
  B(A.get_num_rows(), A.get_num_cols());

  scale(A, s, B);

  return B;
}

//
// Matrix scalar division
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator/(Matrix<T, M, N, ES> const & A, S const & s)
{
  Matrix<typename Promote<S, T>::type, M, N, ES>
  B(A.get_num_rows(), A.get_num_cols());

  divide(A, s, B);

  return B;
}

//
// Scalar matrix division
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator/(S const & s, Matrix<T, M, N, ES> const & A)
{
  Matrix<typename Promote<S, T>::type, M, N, ES>
  B(A.get_num_rows(), A.get_num_cols());

  split(A, s, B);

  return B;
}

//
// Matrix vector product v = A u
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, M, ES>
operator*(Matrix<T, M, N, ES> const & A, Vector<S, N, ES> const & u)
{
  return dot(A, u);
}

//
// Vector matrix product v = u A
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator*(Vector<S, M, ES> const & u, Matrix<T, M, N, ES> const & A)
{
  return dot(u, A);
}

//
// Matrix dot product C = A B
//
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
operator*(Matrix<S, M, P, ES> const & A, Matrix<T, P, N, ES> const & B)
{
  return dot(A, B);
}

//
// Matrix vector product v = A u
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, M, ES>
dot(Matrix<T, M, N, ES> const & A, Vector<S, N, ES> const & u)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = A.get_num_cols();

  assert(u.get_dimension() == num_cols);

  Vector<typename Promote<S, T>::type, M, ES>
  v(num_rows);

  for (Index i = 0; i < num_rows; ++i) {

    typename Promote<S, T>::type
    s = 0.0;

    for (Index p = 0; p < num_cols; ++p) {
      s += A(i, p) * u(p);
    }
    v(i) = s;
  }

  return v;
}

//
// Matrix vector product v = A^T u
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
dot(Vector<S, M, ES> const & u, Matrix<T, M, N, ES> const & A)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = A.get_num_cols();

  assert(u.get_dimension() == num_rows);

  Vector<typename Promote<S, T>::type, N, ES>
  v(num_cols);

  for (Index i = 0; i < num_cols; ++i) {

    typename Promote<S, T>::type
    s = 0.0;

    for (Index p = 0; p < num_rows; ++p) {
      s += A(p, i) * u(p);
    }
    v(i) = s;
  }

  return v;
}

//
// Matrix product C = A B
//
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
dot(Matrix<S, M, P, ES> const & A, Matrix<T, P, N, ES> const & B)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = B.get_num_cols();

  Index const
  dimension = A.get_num_cols();

  assert(B.get_num_rows() == dimension);

  Matrix<typename Promote<S, T>::type, M, N, ES>
  C(num_rows, num_cols);

  for (Index i = 0; i < num_rows; ++i) {
    for (Index j = 0; j < num_cols; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(i, p) * B(p, j);
      }
      C(i, j) = s;
    }
  }

  return C;
}

//
// Matrix product C = A^T B
//
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
t_dot(Matrix<S, P, M, ES> const & A, Matrix<T, P, N, ES> const & B)
{
  Index const
  num_rows = A.get_num_cols();

  Index const
  num_cols = B.get_num_cols();

  Index const
  dimension = A.get_num_rows();

  assert(B.get_num_rows() == dimension);

  Matrix<typename Promote<S, T>::type, M, N, ES>
  C(num_rows, num_cols);

  for (Index i = 0; i < num_rows; ++i) {
    for (Index j = 0; j < num_cols; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(p, i) * B(p, j);
      }
      C(i, j) = s;
    }
  }

  return C;
}

//
// Matrix product C = A B^T
//
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
dot_t(Matrix<S, M, P, ES> const & A, Matrix<T, N, P, ES> const & B)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = B.get_num_rows();

  Index const
  dimension = A.get_num_cols();

  assert(B.get_num_cols() == dimension);

  Matrix<typename Promote<S, T>::type, M, N, ES>
  C(num_rows, num_cols);

  for (Index i = 0; i < num_rows; ++i) {
    for (Index j = 0; j < num_cols; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(i, p) * B(j, p);
      }
      C(i, j) = s;
    }
  }

  return C;
}

//
// Matrix product C = A^T B^T
// \param A tensor
// \param B tensor
// \return a tensor \f$ A^T \cdot B^T \f$
//
template<typename S, typename T, Index M, Index P, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
t_dot_t(Matrix<S, P, M, ES> const & A, Matrix<T, N, P, ES> const & B)
{
  Index const
  num_rows = A.get_num_cols();

  Index const
  num_cols = B.get_num_rows();

  Index const
  dimension = A.get_num_rows();

  assert(B.get_num_cols() == dimension);

  Matrix<typename Promote<S, T>::type, M, N, ES>
  C(num_rows, num_cols);

  for (Index i = 0; i < num_rows; ++i) {
    for (Index j = 0; j < num_cols; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(p, i) * B(j, p);
      }
      C(i, j) = s;
    }
  }

  return C;
}

//
// Matrix double dot product (contraction)
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dotdot(Matrix<S, M, N, ES> const & A, Matrix<T, M, N, ES> const & B)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = B.get_num_cols();

  assert(B.get_num_rows() == num_rows);
  assert(B.get_num_cols() == num_cols);

  typename Promote<S, T>::type
  s = 0.0;

  for (Index p = 0; p < num_rows; ++p) {
    for (Index q = 0; q < num_cols; ++q) {
      s += A(p, q) * B(p, q);
    }
  }

  return s;
}

//
// exterior matrix product
//
template<typename S, typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<typename Promote<S, T>::type, M, N, ES>
matrix(Vector<S, M, ES> const & u, Vector<T, N, ES> const & v)
{
  Index const
  num_rows = u.get_dimension();

  Index const
  num_cols = v.get_dimension();

  Matrix<typename Promote<S, T>::type,M, N, ES>
  A(num_rows, num_cols);

  for (Index i = 0; i < num_rows; ++i) {

    typename Promote<S, T>::type const
    s = u(i);

    for (Index j = 0; j < num_cols; ++j) {
      A(i, j) = s * v(j);
    }
  }

  return A;
}

//
// matrix transpose
//
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, M, N, ES>
transpose(Matrix<T, N, M, ES> const & A)
{
  Index const
  num_rows = A.get_num_cols();

  Index const
  num_cols = A.get_num_rows();

  Matrix<T, M, N, ES>
  B(num_rows, num_cols);

  for (Index i = 0; i < num_rows; ++i) {
    for (Index j = i + 1; j < num_cols; ++j) {
      B(i, j) = A(j, i);
    }
  }

  return B;
}

} // namespace Intrepid

#endif // Intrepid2_MiniTensor_Matrix_i_h
