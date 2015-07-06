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

#if !defined(Intrepid_MiniTensor_Matrix_i_h)
#define Intrepid_MiniTensor_Matrix_i_h

namespace Intrepid {

//
// Constructor that initializes to NaNs
//
template<typename T, Index M, Index N>
inline
Matrix<T, M, N>::Matrix() :
    TensorBase<T, Store>::TensorBase(), rows_(M), cols_(N)
{
  return;
}

template<typename T, Index M, Index N>
inline
Matrix<T, M, N>::Matrix(Index const rows, Index const cols) :
    TensorBase<T, Store>::TensorBase(rows * cols, ORDER),
    rows_(rows), cols_(cols)
{
  return;
}

///
/// Create tensor from a specified value
///
template<typename T, Index M, Index N>
inline
Matrix<T, M, N>::Matrix(ComponentValue const value) :
    TensorBase<T, Store>::TensorBase(M * N, ORDER, value), rows_(M), cols_(N)
{
  return;
}

template<typename T, Index M, Index N>
inline
Matrix<T, M, N>::Matrix(
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
#if defined(HAVE_INTREPID_KOKKOSCORE)
template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1) :
    TensorBase<T, Store>::TensorBase(M * N, ORDER, data, index1),
    rows_(M), cols_(N)
{
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1,
    iType index2) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1, index2),
    rows_(M), cols_(N)
{
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1,
    iType index2,
    iType index3) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1, index2, index3),
    rows_(M), cols_(N)
{
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4) :
    TensorBase<T, Store>::TensorBase(
        N,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4), rows_(M), cols_(N)
{
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5) :
    TensorBase<T, Store>::TensorBase(
        N,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4,
        index5), rows_(M), cols_(N)
{
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5,
    iType index6) :
    TensorBase<T, Store>::TensorBase(
        N,
        ORDER,
        data,
        index1,
        index2,
        index3,
        index4,
        index5,
        index6), rows_(M), cols_(N)
{
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
    Index const rows,
    Index const cols,
    ArrayT & data,
    iType index1) :
    TensorBase<T, Store>::TensorBase(rows * cols, ORDER, data, index1),
    rows_(rows), cols_(cols)
{
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
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

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
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

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
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

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
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

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
Matrix<T, M, N>::Matrix(
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
#endif

//FIXME: Irina D. add kokkos initializations for tensor (data)
template<typename T, Index M, Index N>
inline
Matrix<T, M, N>::Matrix(T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(M * N, ORDER, data_ptr), rows_(M), cols_(N)
{
  return;
}

template<typename T, Index M, Index N>
inline
Matrix<T, M, N>::Matrix(
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
template<typename T, Index M, Index N>
inline
Matrix<T, M, N>::Matrix(Matrix<T, M, N> const & A) :
    TensorBase<T, Store>::TensorBase(A),
    rows_(A.get_num_rows()), cols_(A.get_num_cols())
{
  return;
}

//
// Simple destructor
//
template<typename T, Index M, Index N>
inline
Matrix<T, M, N>::~Matrix()
{
  return;
}

//
// Get dimensions
//
template<typename T, Index M, Index N>
inline
std::pair<Index, Index>
Matrix<T, M, N>::get_dimensions() const
{
  return IS_DYNAMIC == true ?
      std::make_pair(rows_, cols_) :
      std::make_pair(M, N);
}

//
// Get number rows
//
template<typename T, Index M, Index N>
inline
Index
Matrix<T, M, N>::get_num_rows() const
{
  return IS_DYNAMIC == true ? rows_ : M;
}

//
// Get number cols
//
template<typename T, Index M, Index N>
inline
Index
Matrix<T, M, N>::get_num_cols() const
{
  return IS_DYNAMIC == true ? cols_ : M;
}

//
// Set dimensions
//
template<typename T, Index M, Index N>
inline
void
Matrix<T, M, N>::set_dimensions(Index const rows, Index const cols)
{
  if (IS_DYNAMIC == true) {
    TensorBase<T, Store>::set_dimension(rows * cols, ORDER);
    rows_ = rows;
    cols_ = cols;
  }
  else {
    assert(rows == M);
    assert(cols == N);
  }

  return;
}

//
// Indexing for constant tensor
//
template<typename T, Index M, Index N>
inline T const &
Matrix<T, M, N>::operator()(Index const i, Index const j) const
    {
  Matrix<T, M, N> const &
  self = (*this);

  Index const
  num_cols = self.get_num_cols();

  return self[i * num_cols + j];
}

//
//Matrix indexing
//
template<typename T, Index M, Index N>
inline T &
Matrix<T, M, N>::operator()(Index const i, Index const j)
{
  Matrix<T, M, N> &
  self = (*this);

  Index const
  num_cols = self.get_num_cols();

  return self[i * num_cols + j];
}

//
// Fill components with value specification
//
template<typename T, Index M, Index N>
inline
void
Matrix<T, M, N>::fill(ComponentValue const value)
{
  TensorBase<T, Store>::fill(value);
  return;
}

//
// Fill components with value as parameter
//
template<typename T, Index M, Index N>
inline
void
Matrix<T, M, N>::fill(T const & s)
{
  TensorBase<T, Store>::fill(s);
  return;
}

//
// Fill components from array defined by pointer.
//
#if defined(HAVE_INTREPID_KOKKOSCORE)
template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
typename Kokkos::Impl::enable_if<
!Kokkos::Impl::is_same<ArrayT, T*>::value, void>::type
Matrix<T, M, N>::fill(
    ArrayT & data,
    iType index1)
{
  TensorBase<T, Store>::fill(data, index1);
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
void
Matrix<T, M, N>::fill(
    ArrayT & data,
    iType index1,
    iType index2)
{
  TensorBase<T, Store>::fill(data, index1, index2);
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType1, typename iType2, typename iType3>
inline
void
Matrix<T, M, N>::fill(
    ArrayT & data,
    iType1 index1,
    iType2 index2,
    iType3 index3)
{
  TensorBase<T, Store>::fill(data, index1, index2, index3);
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
void
Matrix<T, M, N>::fill(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4)
{
  TensorBase<T, Store>::fill(data, index1, index2, index3, index4);
  return;
}

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
void
Matrix<T, M, N>::fill(
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

template<typename T, Index M, Index N>
template<class ArrayT, typename iType>
inline
void
Matrix<T, M, N>::fill(
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
#endif

template<typename T, Index M, Index N>
inline
void
Matrix<T, M, N>::fill(T const * data_ptr)
{
  TensorBase<T, Store>::fill(data_ptr);
  return;
}

//
// Extract a row as a vector
//
template<typename T, Index M, Index N>
Vector<T, N>
row(Matrix<T, M, N> const & A, Index const i)
{
  Index const
  num_cols = A.get_num_cols();

  Vector<T, N>
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
template<typename T, Index M, Index N>
Vector<T, M>
col(Matrix<T, M, N> const & A, Index const j)
{
  Index const
  num_rows = A.get_num_rows();

  Vector<T, M>
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
template<typename S, typename T, Index M, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
operator+(Matrix<S, M, N> const & A, Matrix<T, M, N> const & B)
{
  Matrix<typename Promote<S, T>::type, M, N>
  C(A.get_num_rows(), A.get_num_cols());

  add(A, B, C);

  return C;
}

//
// Matrix subtraction
//
template<typename S, typename T, Index M, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
operator-(Matrix<S, M, N> const & A, Matrix<T, M, N> const & B)
{
  Matrix<typename Promote<S, T>::type, M, N>
  C(A.get_num_rows(), A.get_num_cols());

  subtract(A, B, C);

  return C;
}

//
// Matrix minus
//
template<typename T, Index M, Index N>
inline
Matrix<T, M, N>
operator-(Matrix<T, M, N> const & A)
{
  Matrix<T, M, N>
  B(A.get_num_rows(), A.get_num_cols());

  minus(A, B);

  return B;
}

//
// Matrix equality
//
template<typename T, Index M, Index N>
inline
bool
operator==(Matrix<T, M, N> const & A, Matrix<T, M, N> const & B)
{
  return equal(A, B);
}

//
// Matrix inequality
//
template<typename T, Index M, Index N>
inline
bool
operator!=(Matrix<T, M, N> const & A, Matrix<T, M, N> const & B)
{
  return not_equal(A, B);
}

//
// Scalar matrix product
//
template<typename S, typename T, Index M, Index N>
inline typename
lazy_disable_if<order_1234<S>, apply_matrix<Promote<S, T>, M, N>>::type
operator*(S const & s, Matrix<T, M, N> const & A)
{
  Matrix<typename Promote<S, T>::type, M, N>
  B(A.get_num_rows(), A.get_num_cols());

  scale(A, s, B);

  return B;
}

//
// Matrix scalar product
//
template<typename S, typename T, Index M, Index N>
inline typename
lazy_disable_if<order_1234<S>, apply_matrix<Promote<S, T>, M, N>>::type
operator*(Matrix<T, M, N> const & A, S const & s)
{
  Matrix<typename Promote<S, T>::type, M, N>
  B(A.get_num_rows(), A.get_num_cols());

  scale(A, s, B);

  return B;
}

//
// Matrix scalar division
//
template<typename S, typename T, Index M, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
operator/(Matrix<T, M, N> const & A, S const & s)
{
  Matrix<typename Promote<S, T>::type, M, N>
  B(A.get_num_rows(), A.get_num_cols());

  divide(A, s, B);

  return B;
}

//
// Scalar matrix division
//
template<typename S, typename T, Index M, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
operator/(S const & s, Matrix<T, M, N> const & A)
{
  Matrix<typename Promote<S, T>::type, M, N>
  B(A.get_num_rows(), A.get_num_cols());

  split(A, s, B);

  return B;
}

//
// Matrix vector product v = A u
//
template<typename S, typename T, Index M, Index N>
inline Vector<typename Promote<S, T>::type, M>
operator*(Matrix<T, M, N> const & A, Vector<S, N> const & u)
{
  return dot(A, u);
}

//
// Vector matrix product v = u A
//
template<typename S, typename T, Index M, Index N>
inline Vector<typename Promote<S, T>::type, N>
operator*(Vector<S, M> const & u, Matrix<T, M, N> const & A)
{
  return dot(u, A);
}

//
// Matrix dot product C = A B
//
template<typename S, typename T, Index M, Index P, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
operator*(Matrix<S, M, P> const & A, Matrix<T, P, N> const & B)
{
  return dot(A, B);
}

//
// Matrix vector product v = A u
//
template<typename S, typename T, Index M, Index N>
inline Vector<typename Promote<S, T>::type, M>
dot(Matrix<T, M, N> const & A, Vector<S, N> const & u)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = A.get_num_cols();

  assert(u.get_dimension() == num_cols);

  Vector<typename Promote<S, T>::type, M>
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
template<typename S, typename T, Index M, Index N>
inline Vector<typename Promote<S, T>::type, N>
dot(Vector<S, M> const & u, Matrix<T, M, N> const & A)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = A.get_num_cols();

  assert(u.get_dimension() == num_rows);

  Vector<typename Promote<S, T>::type, N>
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
template<typename S, typename T, Index M, Index P, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
dot(Matrix<S, M, P> const & A, Matrix<T, P, N> const & B)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = B.get_num_cols();

  Index const
  dimension = A.get_num_cols();

  assert(B.get_num_rows() == dimension);

  Matrix<typename Promote<S, T>::type, M, N>
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
template<typename S, typename T, Index M, Index P, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
t_dot(Matrix<S, P, M> const & A, Matrix<T, P, N> const & B)
{
  Index const
  num_rows = A.get_num_cols();

  Index const
  num_cols = B.get_num_cols();

  Index const
  dimension = A.get_num_rows();

  assert(B.get_num_rows() == dimension);

  Matrix<typename Promote<S, T>::type, M, N>
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
template<typename S, typename T, Index M, Index P, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
dot_t(Matrix<S, M, P> const & A, Matrix<T, N, P> const & B)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = B.get_num_rows();

  Index const
  dimension = A.get_num_cols();

  assert(B.get_num_cols() == dimension);

  Matrix<typename Promote<S, T>::type, M, N>
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
template<typename S, typename T, Index M, Index P, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
t_dot_t(Matrix<S, P, M> const & A, Matrix<T, N, P> const & B)
{
  Index const
  num_rows = A.get_num_cols();

  Index const
  num_cols = B.get_num_rows();

  Index const
  dimension = A.get_num_rows();

  assert(B.get_num_cols() == dimension);

  Matrix<typename Promote<S, T>::type, M, N>
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
template<typename S, typename T, Index M, Index N>
inline
typename Promote<S, T>::type
dotdot(Matrix<S, M, N> const & A, Matrix<T, M, N> const & B)
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
template<typename S, typename T, Index M, Index N>
inline
Matrix<typename Promote<S, T>::type, M, N>
matrix(Vector<S, M> const & u, Vector<T, N> const & v)
{
  Index const
  num_rows = u.get_dimension();

  Index const
  num_cols = v.get_dimension();

  Matrix<typename Promote<S, T>::type, N>
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
template<typename T, Index M, Index N>
inline
Matrix<T, M, N>
transpose(Matrix<T, N, M> const & A)
{
  Index const
  num_rows = A.get_num_cols();

  Index const
  num_cols = A.get_num_rows();

  Matrix<T, M, N>
  B(num_rows, num_cols);

  for (Index i = 0; i < num_rows; ++i) {
    for (Index j = i + 1; j < num_cols; ++j) {
      B(i, j) = A(j, i);
    }
  }

  return B;
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Matrix_i_h
