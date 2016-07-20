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

#if !defined(Intrepid2_MiniTensor_Tensor_i_h)
#define Intrepid2_MiniTensor_Tensor_i_h

namespace Intrepid2 {

//
// Constructor that initializes to NaNs
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor() :
    TensorBase<T, Store>::TensorBase()
{
  set_dimension(N);
  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(Index const dimension) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

///
/// Create tensor from a specified value
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(ComponentValue const value) :
    TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(Index const dimension, ComponentValue const value) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

//
//  Create tensor from array
//
//
template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(ArrayT & data, iType index1) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1)
{
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1,
    iType index2) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1, index2)
{
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1,
    iType index2,
    iType index3) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1, index2, index3)
{
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
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
        index4)
{
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
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
        index5)
{
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
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
        index6)
{
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(Index const dimension, ArrayT & data, iType index1) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data, index1)
{
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    Index const dimension,
    typename Kokkos::Impl::enable_if<
    !Kokkos::Impl::is_same<ArrayT, Index>::value, ArrayT>::type & data,
    iType index1,
    iType index2) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data, index1, index2)
{
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    Index const dimension,
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3) :
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

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    Index const dimension,
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4) :
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

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    Index const dimension,
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5) :
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

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    Index const dimension,
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4,
    iType index5,
    iType index6) :
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

//FIXME: Irina D. add kokkos initializations for tensor (data)
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(Index const dimension, T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}
//
// Copy constructor
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(Tensor<T, N, ES> const & A) :
    TensorBase<T, Store>::TensorBase(A)
{
  return;
}

//
// Create tensor specifying components
// \param  s00, s01, ... components in the R^2 canonical basis
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    T const & s00, T const & s01,
    T const & s10, T const & s11)
{
  Tensor<T, N, ES> &
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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    T const & s00, T const & s01, T const & s02,
    T const & s10, T const & s11, T const & s12,
    T const & s20, T const & s21, T const & s22)
{
  Tensor<T, N, ES> &
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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(T const * data_ptr, ComponentOrder const component_order)
{
  assert(data_ptr != NULL);

  fill(data_ptr, component_order);

  return;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(
    Index const dimension,
    T const * data_ptr,
    ComponentOrder const component_order)
{
  assert(data_ptr != NULL);

  Tensor<T, N, ES> &
  self = (*this);

  self.set_dimension(dimension);

  fill(data_ptr, component_order);

  return;
}

//
// 2nd-order tensor from 4th-order tensor
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::Tensor(Tensor4<T, dimension_sqrt<N>::value, ES> const & A)
{
  Index const
  dimension_4th = A.get_dimension();

  Index const
  dimension_2nd = dimension_4th * dimension_4th;

  Tensor<T, N, ES> &
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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>::~Tensor()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Index
Tensor<T, N, ES>::get_dimension() const
{
  return TensorBase<T, Store>::get_dimension();
}

//
// Set dimension
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::set_dimension(Index const dimension)
{
  if (IS_DYNAMIC == false) {
    assert(dimension <= N);
  }
  TensorBase<T, Store>::set_dimension(dimension, ORDER);
  return;
}

//
// Indexing for constant tensor
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T const &
Tensor<T, N, ES>::operator()(Index const i, Index const j) const
    {
  Tensor<T, N, ES> const &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[i * dimension + j];
}

//
//Tensor indexing
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T &
Tensor<T, N, ES>::operator()(Index const i, Index const j)
{
  Tensor<T, N, ES> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[i * dimension + j];
}

//
// Fill components with value specification
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::fill(ComponentValue const value)
{
  TensorBase<T, Store>::fill(value);
  return;
}

//
// Fill components with value as parameter
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::fill(T const & s)
{
  TensorBase<T, Store>::fill(s);
  return;
}

//
// Fill components from array defined by pointer.
//
template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
typename Kokkos::Impl::enable_if<
!Kokkos::Impl::is_same<ArrayT, T*>::value, void>::type
Tensor<T, N, ES>::fill(ArrayT & data, 
            iType index1)
{
  TensorBase<T, Store>::fill(data, index1);
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::fill(ArrayT & data, iType index1, iType index2)
{
  TensorBase<T, Store>::fill(data, index1, index2);
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType1, typename iType2, typename iType3>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::fill(ArrayT & data, iType1 index1, iType2 index2, iType3 index3)
{
  TensorBase<T, Store>::fill(data, index1, index2, index3);
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::fill(
    ArrayT & data,
    iType index1,
    iType index2,
    iType index3,
    iType index4)
{
  TensorBase<T, Store>::fill(data, index1, index2, index3, index4);
  return;
}

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::fill(
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

template<typename T, Index N,  typename ES>
template<class ArrayT, typename iType>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::fill(
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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::fill(T const * data_ptr)
{
  TensorBase<T, Store>::fill(data_ptr);
  return;
}

//
// Fill components from array defined by pointer.
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
Tensor<T, N, ES>::fill(T const * data_ptr, ComponentOrder const component_order)
{
  assert(data_ptr != NULL);

  Tensor<T, N, ES> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  switch (dimension) {

  default:
    TensorBase<T, Store>::fill(data_ptr);
    break;

  case 3:

    switch (component_order) {

    case CANONICAL:
      TensorBase<T, Store>::fill(data_ptr);
      break;

    case SIERRA_FULL:
      //  0  1  2  3  4  5  6  7  8
      // XX YY ZZ XY YZ ZX YX ZY XZ
      //  0  4  8  1  5  6  3  7  2
      self[0] = data_ptr[0];
      self[4] = data_ptr[1];
      self[8] = data_ptr[2];

      self[1] = data_ptr[3];
      self[5] = data_ptr[4];
      self[6] = data_ptr[5];

      self[3] = data_ptr[6];
      self[7] = data_ptr[7];
      self[2] = data_ptr[8];
      break;

    case SIERRA_SYMMETRIC:
      self[0] = data_ptr[0];
      self[4] = data_ptr[1];
      self[8] = data_ptr[2];

      self[1] = data_ptr[3];
      self[5] = data_ptr[4];
      self[6] = data_ptr[5];

      self[3] = data_ptr[3];
      self[7] = data_ptr[4];
      self[2] = data_ptr[5];
      break;

    default:
#if defined(KOKKOS_HAVE_CUDA)
     Kokkos::abort("ERROR (Tensor::fill(): Unknown component order.");
#else
      std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
      std::cerr << std::endl;
      std::cerr << "Unknown component order.";
      std::cerr << std::endl;
      exit(1);
#endif
      break;

    }

    break;
  }

  return;
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
// Sort and index in descending order. Useful for ordering singular values
// and eigenvalues and corresponding vectors in the respective decompositions.
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Vector<T, N, ES>, Tensor<T, N, ES>>
sort_permutation(Vector<T, N, ES> const & u)
{

  Index const
  dimension = u.get_dimension();

  std::vector <std::pair<T, Index>>
  s(dimension);

  for (Index i = 0; i < dimension; ++i) {
    s[i].first = u(i);
    s[i].second = i;
  }

  std::sort(s.begin(), s.end(), greater_than<std::pair<T, Index>>);

  Vector<T, N, ES> v(dimension);

  Tensor<T, N, ES>
  P = zero<T, N, ES>(dimension);

  for (Index i = 0; i < dimension; ++i) {
    v(i) = s[i].first;
    P(s[i].second, i) = 1.0;
  }

  return std::make_pair(v, P);

}

//
// Extract a row as a vector
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
row(Tensor<T, N, ES> const & A, Index const i)
{
  Index const
  dimension = A.get_dimension();

  Vector<T, N, ES>
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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
col(Tensor<T, N, ES> const & A, Index const j)
{
  Index const
  dimension = A.get_dimension();

  Vector<T, N, ES>
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
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator+(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B)
{
  Tensor<typename Promote<S, T>::type, N, ES>
  C(A.get_dimension());

  add(A, B, C);

  return C;
}

//
// Tensor subtraction
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator-(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B)
{
  Tensor<typename Promote<S, T>::type, N, ES>
  C(A.get_dimension());

  subtract(A, B, C);

  return C;
}

//
// Tensor minus
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
operator-(Tensor<T, N, ES> const & A)
{
  Tensor<T, N, ES>
  B(A.get_dimension());

  minus(A, B);

  return B;
}

//
// Tensor equality
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor<T, N, ES> const & A, Tensor<T, N, ES> const & B)
{
  return equal(A, B);
}

//
// Tensor inequality
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor<T, N, ES> const & A, Tensor<T, N, ES> const & B)
{
  return not_equal(A, B);
}

//
// Scalar tensor product
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename
lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N, ES>>::type
operator*(S const & s, Tensor<T, N, ES> const & A)
{
  Tensor<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// Tensor scalar product
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename
lazy_disable_if<order_1234<S>, apply_tensor<Promote<S, T>, N, ES>>::type
operator*(Tensor<T, N, ES> const & A, S const & s)
{
  Tensor<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// Tensor scalar division
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator/(Tensor<T, N, ES> const & A, S const & s)
{
  Tensor<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  divide(A, s, B);

  return B;
}

//
// Scalar tensor division
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator/(S const & s, Tensor<T, N, ES> const & A)
{
  Tensor<typename Promote<S, T>::type, N, ES>
  B(A.get_dimension());

  split(A, s, B);

  return B;
}

//
// Tensor vector product v = A u
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator*(Tensor<T, N, ES> const & A, Vector<S, N, ES> const & u)
{
  return dot(A, u);
}

//
// Vector tensor product v = u A
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
operator*(Vector<S, N, ES> const & u, Tensor<T, N, ES> const & A)
{
  return dot(u, A);
}

//
// Tensor dot product C = A B
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
operator*(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B)
{
  return dot(A, B);
}

//
// Tensor vector product v = A u
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
dot(Tensor<T, N, ES> const & A, Vector<S, N, ES> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

   Vector<typename Promote<S, T>::type, N, ES>
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
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
dot(Vector<S, N, ES> const & u, Tensor<T, N, ES> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

   Vector<typename Promote<S, T>::type, N, ES>
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
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dot(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N, ES>
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
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
t_dot(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N, ES>
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
// Tensor tensor product C = A B^T
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dot_t(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N, ES>
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
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
t_dot_t(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N, ES>
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
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename Promote<S, T>::type
dotdot(Tensor<S, N, ES> const & A, Tensor<T, N, ES> const & B)
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
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dyad(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v)
{
  Index const
  dimension = u.get_dimension();

  assert(v.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N, ES>
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
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
bun(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v)
{
  return dyad(u, v);
}

//
// tensor product
//
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
tensor(Vector<S, N, ES> const & u, Vector<T, N, ES> const & v)
{
  return dyad(u, v);
}

//
// diagonal tensor from vector
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
diag(Vector<T, N, ES> const & v)
{
  Index const
  dimension = v.get_dimension();

  Tensor<T, N, ES>
  A = zero<T, N, ES>(dimension);

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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
diag(Tensor<T, N, ES> const & A)
{
  Index const
  dimension = A.get_dimension();

  Vector<T, N, ES>
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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
zero()
{
  return Tensor<T, N, ES>(N, ZEROS);
}

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
zero(Index const dimension)
{
  return Tensor<T, DYNAMIC, ES>(dimension, ZEROS);
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
zero(Index const dimension)
{
  return Tensor<T, N, ES>(dimension, ZEROS);
}

// Local utility functions
namespace {

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void ones_in_diagonal(Tensor<T, N, ES> & A)
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

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void fill_levi_civita(Tensor<T, N, ES> & A)
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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
identity()
{
  Tensor<T, N, ES>
  A(N, ZEROS);

  ones_in_diagonal(A);

  return A;
}

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
identity(Index const dimension)
{
  Tensor<T, DYNAMIC, ES>
  A(dimension, ZEROS);

  ones_in_diagonal(A);

  return A;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
identity(Index const dimension)
{

  Tensor<T, N, ES>
  A(dimension, ZEROS);

  ones_in_diagonal(A);

  return A;
}

//
// 2nd-order identity tensor, à la Matlab
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
eye()
{
  return identity<T, N, ES>();
}

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
eye(Index const dimension)
{
  return identity<T, DYNAMIC,ES>(dimension);
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
eye(Index const dimension)
{
  return identity<T, N, ES>(dimension);
}

//
// Levi-Civita symbol
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
levi_civita_2()
{
  Tensor<T, N, ES>
  A(N, ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
levi_civita_2(Index const dimension)
{
  Tensor<T, DYNAMIC, ES>
  A(dimension, ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
levi_civita_2(Index const dimension)
{

  Tensor<T, DYNAMIC, ES>
  A(dimension, ZEROS);

  fill_levi_civita(A);

  return A;
}

//
// Permutation symbol
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
permutation_2()
{
  return levi_civita_2<T, N, ES>();
}

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
permutation_2(Index const dimension)
{
  return levi_civita_2<T, DYNAMIC,ES>(dimension);
}


template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
permutation_2(Index const dimension)
{
  return levi_civita_2<T, N,ES>(dimension);
}

//
// Alternating symbol
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
alternator_2()
{
  return levi_civita_2<T, N,ES>();
}

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, DYNAMIC, ES> const
alternator_2(Index const dimension)
{
  return levi_civita_2<T, DYNAMIC, ES>(dimension);
}

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES> const
alternator_2(Index const dimension)
{
  return levi_civita_2<T, N,ES>(dimension);
}

//
// 2nd-order tensor transpose
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
transpose(Tensor<T, N, ES> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N, ES>
  B = A;

  switch (dimension) {
  default:
    for (Index i = 0; i < dimension; ++i) {
      for (Index j = i + 1; j < dimension; ++j) {
        std::swap(B(i, j), B(j, i));
      }
    }
    break;

  case 3:
    std::swap(B(0, 1), B(1, 0));
    std::swap(B(0, 2), B(2, 0));

    std::swap(B(1, 2), B(2, 1));

    break;

  case 2:
    std::swap(B(0, 1), B(1, 0));

    break;
  }

  return B;
}

//
// symmetric part of 2nd-order tensor
// \return \f$ \frac{1}{2}(A + A^T) \f$
//
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
sym(Tensor<T, N, ES> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N, ES>
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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
skew(Tensor<T, N, ES> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor<T, N, ES>
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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
skew(Vector<T, N, ES> const & u)
{
  Index const
  dimension = u.get_dimension();

  Tensor<T, N, ES>
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
#if defined(KOKKOS_HAVE_CUDA)
    Kokkos::abort("ERROR(Tensor::skew): Skew from vector undefined for R");
#else
    std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
    std::cerr << std::endl;
    std::cerr << "Skew from vector undefined for R^" << N;
    std::cerr << std::endl;
    exit(1);
    break;
#endif
  }

  return A;
}

} // namespace Intrepid

#endif // Intrepid2_MiniTensor_Tensor_i_h
