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

#if !defined(Intrepid_MiniTensor_Tensor_i_h)
#define Intrepid_MiniTensor_Tensor_i_h

namespace Intrepid {

namespace {

template<typename S>
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
template<typename T>
std::pair<Vector<T>, Tensor<T> >
sort_permutation(Vector<T> const & u)
{

  Index const
  N = u.get_dimension();

  std::vector<std::pair<T, Index > >
  s(N);

  for (Index i = 0; i < N; ++i) {
    s[i].first = u(i);
    s[i].second = i;
  }

  std::sort(s.begin(), s.end(), greater_than< std::pair<T, Index > > );

  Vector<T> v(N);

  Tensor<T>
  P = zero<T>(N);

  for (Index i = 0; i < N; ++i) {
    v(i) = s[i].first;
    P(s[i].second, i) = 1.0;
  }

  return std::make_pair(v, P);

}

//
// Fill components with value.
//
template<typename T>
inline
void
Tensor<T>::fill(ComponentValue value)
{
  TensorBase<T>::fill(value);
  return;
}

//
// Fill components from array defined by pointer.
//
template<typename T>
inline
void
Tensor<T>::fill(T const * data_ptr, ComponentOrder const component_order)
{
  assert(data_ptr != NULL);

  Tensor<T> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  switch (dimension) {

    default:
      TensorBase<T>::fill(data_ptr);
      break;

    case 3:

      switch (component_order) {

        case CANONICAL:
          TensorBase<T>::fill(data_ptr);
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
          std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
          std::cerr << std::endl;
          std::cerr << "Unknown component order.";
          std::cerr << std::endl;
          exit(1);
          break;

      }

      break;
  }

  return;
}

//
// Default constructor
//
template<typename T>
inline
Tensor<T>::Tensor() :
TensorBase<T>::TensorBase()
{
  return;
}

//
// Constructor that initializes to NaNs
//
template<typename T>
inline
Tensor<T>::Tensor(Index const dimension) :
TensorBase<T>::TensorBase(dimension, order)
{
  return;
}

///
/// Create tensor from a specified value
///
template<typename T>
inline
Tensor<T>::Tensor(Index const dimension, ComponentValue value) :
TensorBase<T>::TensorBase(dimension, order, value)
{
  return;
}

//
// Create tensor from a scalar
//
template<typename T>
inline
Tensor<T>::Tensor(Index const dimension, T const & s) :
TensorBase<T>::TensorBase(dimension, order, s)
{
  return;
}

//
// Create tensor specifying components
// \param  s00, s01, ... components in the R^2 canonical basis
//
template<typename T>
inline
Tensor<T>::Tensor(
    T const & s00, T const & s01,
    T const & s10, T const & s11)
{
  Tensor<T> &
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
template<typename T>
inline
Tensor<T>::Tensor(
    T const & s00, T const & s01, T const & s02,
    T const & s10, T const & s11, T const & s12,
    T const & s20, T const & s21, T const & s22)
{
  Tensor<T> &
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
//  Create tensor from array
//
template<typename T>
inline
Tensor<T>::Tensor(Index const dimension, T const * data_ptr) :
TensorBase<T>::TensorBase(dimension, order, data_ptr)
{
  return;
}

//
//  Create tensor from array with component order
//
template<typename T>
inline
Tensor<T>::Tensor(
    Index const dimension,
    T const * data_ptr,
    ComponentOrder const component_order)
{
  assert(data_ptr != NULL);

  Tensor<T> &
  self = (*this);

  self.set_dimension(dimension);

  fill(data_ptr, component_order);

  return;
}

//
// Copy constructor
//
template<typename T>
inline
Tensor<T>::Tensor(Tensor<T> const & A) :
TensorBase<T>::TensorBase(A)
{
  return;
}

//
// 2nd-order tensor from 4th-order tensor
//
template<typename T>
inline
Tensor<T>::Tensor(Tensor4<T> const & A)
{
  Index const
  dimension_4th = A.get_dimension();

  Index const
  dimension_2nd = dimension_4th * dimension_4th;

  Tensor<T> &
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
template<typename T>
inline
Tensor<T>::~Tensor()
{
  return;
}

//
// Indexing for constant tensor
//
template<typename T>
inline
T const &
Tensor<T>::operator()(Index const i, Index const j) const
{
  Tensor<T> const &
  self = (*this);

  Index const
  N = self.get_dimension();

  assert(i < N);
  assert(j < N);

  return self[i * N + j];
}

//
//Tensor indexing
//
template<typename T>
inline
T &
Tensor<T>::operator()(Index const i, Index const j)
{
  Tensor<T> &
  self = (*this);

  Index const
  N = self.get_dimension();

  assert(i < N);
  assert(j < N);

  return self[i * N + j];
}

//
// Tensor addition
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
operator+(Tensor<S> const & A, Tensor<T> const & B)
{
  Tensor<typename Promote<S, T>::type>
  C;

  add(A, B, C);

  return C;
}

//
// Tensor subtraction
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
operator-(Tensor<S> const & A, Tensor<T> const & B)
{
  Tensor<typename Promote<S, T>::type>
  C;

  subtract(A, B, C);

  return C;
}

//
// Tensor minus
//
template<typename T>
inline
Tensor<T>
operator-(Tensor<T> const & A)
{
  Tensor<T>
  B;

  minus(A, B);

  return B;
}

//
// Tensor equality
//
template<typename T>
inline
bool
operator==(Tensor<T> const & A, Tensor<T> const & B)
{
  return equal(A, B);
}

//
// Tensor inequality
//
template<typename T>
inline
bool
operator!=(Tensor<T> const & A, Tensor<T> const & B)
{
  return not_equal(A, B);
}

//
// Scalar tensor product
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_tensor< Promote<S,T> > >::type
operator*(S const & s, Tensor<T> const & A)
{
  Tensor<typename Promote<S, T>::type>
  B;

  scale(A, s, B);

  return B;
}

//
// Tensor scalar product
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_tensor< Promote<S,T> > >::type
operator*(Tensor<T> const & A, S const & s)
{
  Tensor<typename Promote<S, T>::type>
  B;

  scale(A, s, B);

  return B;
}

//
// Tensor scalar division
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
operator/(Tensor<T> const & A, S const & s)
{
  Tensor<typename Promote<S, T>::type>
  B;

  divide(A, s, B);

  return B;
}

//
// Tensor vector product v = A u
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
operator*(Tensor<T> const & A, Vector<S> const & u)
{
  return dot(A, u);
}

//
// Vector tensor product v = u A
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
operator*(Vector<S> const & u, Tensor<T> const & A)
{
  return dot(u, A);
}

//
// Tensor dot product C = A B
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
operator*(Tensor<S> const & A, Tensor<T> const & B)
{
  return dot(A, B);
}

//
// Extract a row as a vector
//
template<typename T>
Vector<T>
row(Tensor<T> const & A, Index const i)
{
  Index const
  N = A.get_dimension();

  Vector<T>
  v(N);

  switch (N) {
    default:
      for (Index j = 0; j < N; ++j) {
        v(j) = A(i,j);
      }
      break;

    case 2:
      v(0) = A(i,0);
      v(1) = A(i,1);
      break;

    case 3:
      v(0) = A(i,0);
      v(1) = A(i,1);
      v(2) = A(i,2);
      break;
  }

  return v;
}

//
// Extract a column as a vector
//
template<typename T>
Vector<T>
col(Tensor<T> const & A, Index const j)
{
  Index const
  N = A.get_dimension();

  Vector<T>
  v(N);

  switch (N) {
    default:
      for (Index i = 0; i < N; ++i) {
        v(i) = A(i,j);
      }
      break;

    case 2:
      v(0) = A(0,j);
      v(1) = A(1,j);
      break;

    case 3:
      v(0) = A(0,j);
      v(1) = A(1,j);
      v(2) = A(2,j);
      break;
  }

  return v;
}


//
// R^N tensor vector product v = A u
// \param A tensor
// \param u vector
// \return \f$ A u \f$
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
dot(Tensor<T> const & A, Vector<S> const & u)
{
  Index const
  N = A.get_dimension();

  assert(u.get_dimension() == N);

  Vector<typename Promote<S, T>::type>
  v(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < N; ++p) {
          s += A(i, p) * u(p);
        }
        v(i) = s;
      }
      break;

    case 3:
      v(0) = A(0,0)*u(0) + A(0,1)*u(1) + A(0,2)*u(2);
      v(1) = A(1,0)*u(0) + A(1,1)*u(1) + A(1,2)*u(2);
      v(2) = A(2,0)*u(0) + A(2,1)*u(1) + A(2,2)*u(2);
      break;

    case 2:
      v(0) = A(0,0)*u(0) + A(0,1)*u(1);
      v(1) = A(1,0)*u(0) + A(1,1)*u(1);
      break;

  }

  return v;
}

//
// R^N vector tensor product v = u A
// \param A tensor
// \param u vector
// \return \f$ u A = A^T u \f$
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
dot(Vector<S> const & u, Tensor<T> const & A)
{
  Index const
  N = A.get_dimension();

  assert(u.get_dimension() == N);

  Vector<typename Promote<S, T>::type>
  v(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < N; ++p) {
          s += A(p, i) * u(p);
        }
        v(i) = s;
      }
      break;

    case 3:
      v(0) = A(0,0)*u(0) + A(1,0)*u(1) + A(2,0)*u(2);
      v(1) = A(0,1)*u(0) + A(1,1)*u(1) + A(2,1)*u(2);
      v(2) = A(0,2)*u(0) + A(1,2)*u(1) + A(2,2)*u(2);
      break;

    case 2:
      v(0) = A(0,0)*u(0) + A(1,0)*u(1);
      v(1) = A(0,1)*u(0) + A(1,1)*u(1);
      break;

  }

  return v;
}

//
// R^N tensor tensor product C = A B
// \param A tensor
// \param B tensor
// \return a tensor \f$ A \cdot B \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
dot(Tensor<S> const & A, Tensor<T> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  Tensor<typename Promote<S, T>::type>
  C(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < N; ++p) {
            s += A(i, p) * B(p, j);
          }
          C(i, j) = s;
        }
      }
      break;

    case 3:
      C(0,0) = A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0);
      C(0,1) = A(0,0)*B(0,1) + A(0,1)*B(1,1) + A(0,2)*B(2,1);
      C(0,2) = A(0,0)*B(0,2) + A(0,1)*B(1,2) + A(0,2)*B(2,2);

      C(1,0) = A(1,0)*B(0,0) + A(1,1)*B(1,0) + A(1,2)*B(2,0);
      C(1,1) = A(1,0)*B(0,1) + A(1,1)*B(1,1) + A(1,2)*B(2,1);
      C(1,2) = A(1,0)*B(0,2) + A(1,1)*B(1,2) + A(1,2)*B(2,2);

      C(2,0) = A(2,0)*B(0,0) + A(2,1)*B(1,0) + A(2,2)*B(2,0);
      C(2,1) = A(2,0)*B(0,1) + A(2,1)*B(1,1) + A(2,2)*B(2,1);
      C(2,2) = A(2,0)*B(0,2) + A(2,1)*B(1,2) + A(2,2)*B(2,2);
      break;

    case 2:
      C(0,0) = A(0,0)*B(0,0) + A(0,1)*B(1,0);
      C(0,1) = A(0,0)*B(0,1) + A(0,1)*B(1,1);

      C(1,0) = A(1,0)*B(0,0) + A(1,1)*B(1,0);
      C(1,1) = A(1,0)*B(0,1) + A(1,1)*B(1,1);
      break;

  }

  return C;
}

//
// R^N tensor tensor product C = A^T B
// \param A tensor
// \param B tensor
// \return a tensor \f$ A^T \cdot B \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
t_dot(Tensor<S> const & A, Tensor<T> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  Tensor<typename Promote<S, T>::type>
  C(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < N; ++p) {
            s += A(p, i) * B(p, j);
          }
          C(i, j) = s;
        }
      }
      break;

    case 3:
      C(0,0) = A(0,0)*B(0,0) + A(1,0)*B(1,0) + A(2,0)*B(2,0);
      C(0,1) = A(0,0)*B(0,1) + A(1,0)*B(1,1) + A(2,0)*B(2,1);
      C(0,2) = A(0,0)*B(0,2) + A(1,0)*B(1,2) + A(2,0)*B(2,2);

      C(1,0) = A(0,1)*B(0,0) + A(1,1)*B(1,0) + A(2,1)*B(2,0);
      C(1,1) = A(0,1)*B(0,1) + A(1,1)*B(1,1) + A(2,1)*B(2,1);
      C(1,2) = A(0,1)*B(0,2) + A(1,1)*B(1,2) + A(2,1)*B(2,2);

      C(2,0) = A(0,2)*B(0,0) + A(1,2)*B(1,0) + A(2,2)*B(2,0);
      C(2,1) = A(0,2)*B(0,1) + A(1,2)*B(1,1) + A(2,2)*B(2,1);
      C(2,2) = A(0,2)*B(0,2) + A(1,2)*B(1,2) + A(2,2)*B(2,2);
      break;

    case 2:
      C(0,0) = A(0,0)*B(0,0) + A(1,0)*B(1,0);
      C(0,1) = A(0,0)*B(0,1) + A(1,0)*B(1,1);

      C(1,0) = A(0,1)*B(0,0) + A(1,1)*B(1,0);
      C(1,1) = A(0,1)*B(0,1) + A(1,1)*B(1,1);
      break;

  }

  return C;
}

//
// R^N tensor tensor product C = A B^T
// \param A tensor
// \param B tensor
// \return a tensor \f$ A \cdot B^T \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
dot_t(Tensor<S> const & A, Tensor<T> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  Tensor<typename Promote<S, T>::type>
  C(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < N; ++p) {
            s += A(i, p) * B(j, p);
          }
          C(i, j) = s;
        }
      }
      break;

    case 3:
      C(0,0) = A(0,0)*B(0,0) + A(0,1)*B(0,1) + A(0,2)*B(0,2);
      C(0,1) = A(0,0)*B(1,0) + A(0,1)*B(1,1) + A(0,2)*B(1,2);
      C(0,2) = A(0,0)*B(2,0) + A(0,1)*B(2,1) + A(0,2)*B(2,2);

      C(1,0) = A(1,0)*B(0,0) + A(1,1)*B(0,1) + A(1,2)*B(0,2);
      C(1,1) = A(1,0)*B(1,0) + A(1,1)*B(1,1) + A(1,2)*B(1,2);
      C(1,2) = A(1,0)*B(2,0) + A(1,1)*B(2,1) + A(1,2)*B(2,2);

      C(2,0) = A(2,0)*B(0,0) + A(2,1)*B(0,1) + A(2,2)*B(0,2);
      C(2,1) = A(2,0)*B(1,0) + A(2,1)*B(1,1) + A(2,2)*B(1,2);
      C(2,2) = A(2,0)*B(2,0) + A(2,1)*B(2,1) + A(2,2)*B(2,2);
      break;

    case 2:
      C(0,0) = A(0,0)*B(0,0) + A(0,1)*B(0,1);
      C(0,1) = A(0,0)*B(1,0) + A(0,1)*B(1,1);

      C(1,0) = A(1,0)*B(0,0) + A(1,1)*B(0,1);
      C(1,1) = A(1,0)*B(1,0) + A(1,1)*B(1,1);
      break;

  }

  return C;
}

//
// R^N tensor tensor product C = A^T B^T
// \param A tensor
// \param B tensor
// \return a tensor \f$ A^T \cdot B^T \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
t_dot_t(Tensor<S> const & A, Tensor<T> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  Tensor<typename Promote<S, T>::type>
  C(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < N; ++p) {
            s += A(p, i) * B(j, p);
          }
          C(i, j) = s;
        }
      }
      break;

    case 3:
      C(0,0) = A(0,0)*B(0,0) + A(1,0)*B(0,1) + A(2,0)*B(0,2);
      C(0,1) = A(0,0)*B(1,0) + A(1,0)*B(1,1) + A(2,0)*B(1,2);
      C(0,2) = A(0,0)*B(2,0) + A(1,0)*B(2,1) + A(2,0)*B(2,2);

      C(1,0) = A(0,1)*B(0,0) + A(1,1)*B(0,1) + A(2,1)*B(0,2);
      C(1,1) = A(0,1)*B(1,0) + A(1,1)*B(1,1) + A(2,1)*B(1,2);
      C(1,2) = A(0,1)*B(2,0) + A(1,1)*B(2,1) + A(2,1)*B(2,2);

      C(2,0) = A(0,2)*B(0,0) + A(1,2)*B(0,1) + A(2,2)*B(0,2);
      C(2,1) = A(0,2)*B(1,0) + A(1,2)*B(1,1) + A(2,2)*B(1,2);
      C(2,2) = A(0,2)*B(2,0) + A(1,2)*B(2,1) + A(2,2)*B(2,2);
      break;

    case 2:
      C(0,0) = A(0,0)*B(0,0) + A(1,0)*B(0,1);
      C(0,1) = A(0,0)*B(1,0) + A(1,0)*B(1,1);

      C(1,0) = A(0,1)*B(0,0) + A(1,1)*B(0,1);
      C(1,1) = A(0,1)*B(1,0) + A(1,1)*B(1,1);
      break;

  }

  return C;
}

//
// R^N tensor tensor double dot product (contraction)
// \param A tensor
// \param B tensor
// \return a scalar \f$ A : B \f$
//
template<typename S, typename T>
inline
typename Promote<S, T>::type
dotdot(Tensor<S> const & A, Tensor<T> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  typename Promote<S, T>::type
  s = 0.0;

  switch (N) {

    default:
      for (Index p = 0; p < N; ++p) {
        for (Index q = 0; q < N; ++q) {
          s += A(p, q) * B(p, q);
        }
      }
      break;

    case 3:
      s+= A(0,0)*B(0,0) + A(0,1)*B(0,1) + A(0,2)*B(0,2);
      s+= A(1,0)*B(1,0) + A(1,1)*B(1,1) + A(1,2)*B(1,2);
      s+= A(2,0)*B(2,0) + A(2,1)*B(2,1) + A(2,2)*B(2,2);
      break;

    case 2:
      s+= A(0,0)*B(0,0) + A(0,1)*B(0,1);
      s+= A(1,0)*B(1,0) + A(1,1)*B(1,1);
      break;

  }

  return s;
}

//
// R^N dyad
// \param u vector
// \param v vector
// \return \f$ u \otimes v \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
dyad(Vector<S> const & u, Vector<T> const & v)
{
  Index const
  N = u.get_dimension();

  assert(v.get_dimension() == N);

  Tensor<typename Promote<S, T>::type>
  A(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {

        typename Promote<S, T>::type const
        s = u(i);

        for (Index j = 0; j < N; ++j) {
          A(i, j) = s * v(j);
        }
      }
      break;

    case 3:
      A(0,0) = u(0) * v(0);
      A(0,1) = u(0) * v(1);
      A(0,2) = u(0) * v(2);

      A(1,0) = u(1) * v(0);
      A(1,1) = u(1) * v(1);
      A(1,2) = u(1) * v(2);

      A(2,0) = u(2) * v(0);
      A(2,1) = u(2) * v(1);
      A(2,2) = u(2) * v(2);
      break;

    case 2:
      A(0,0) = u(0) * v(0);
      A(0,1) = u(0) * v(1);

      A(1,0) = u(1) * v(0);
      A(1,1) = u(1) * v(1);
      break;

  }

  return A;
}

//
// R^N bun operator, just for Jay
// \param u vector
// \param v vector
// \return \f$ u \otimes v \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
bun(Vector<S> const & u, Vector<T> const & v)
{
  return dyad(u, v);
}

//
// R^N tensor product
// \param u vector
// \param v vector
// \return \f$ u \otimes v \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
tensor(Vector<S> const & u, Vector<T> const & v)
{
  return dyad(u, v);
}

//
// R^N diagonal tensor from vector
// \param v vector
// \return A = diag(v)
//
template<typename T>
Tensor<T>
diag(Vector<T> const & v)
{
  Index const
  N = v.get_dimension();

  Tensor<T>
  A = zero<T>(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        A(i, i) = v(i);
      }
      break;

    case 3:
      A(0,0) = v(0);
      A(1,1) = v(1);
      A(2,2) = v(2);
      break;

    case 2:
      A(0,0) = v(0);
      A(1,1) = v(1);
      break;

  }

  return A;
}

//
// R^N diagonal of tensor in a vector
// \param A tensor
// \return v = diag(A)
//
template<typename T>
Vector<T>
diag(Tensor<T> const & A)
{
  Index const
  N = A.get_dimension();

  Vector<T>
  v(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        v(i) = A(i, i);
      }
      break;

    case 3:
      v(0) = A(0,0);
      v(1) = A(1,1);
      v(2) = A(2,2);
      break;

    case 2:
      v(0) = A(0,0);
      v(1) = A(1,1);
      break;

  }

  return v;
}

//
// Zero 2nd-order tensor
// All components are zero
//
template<typename T>
inline
Tensor<T> const
zero(Index const N)
{
  return Tensor<T>(N, ZEROS);
}

//
// R^N 2nd-order identity tensor
//
template<typename T>
inline
Tensor<T> const
identity(Index const N)
{
  Tensor<T>
  A(N, ZEROS);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        A(i, i) = 1.0;
      }
      break;

    case 3:
      A(0,0) = 1.0;
      A(1,1) = 1.0;
      A(2,2) = 1.0;
      break;

    case 2:
      A(0,0) = 1.0;
      A(1,1) = 1.0;
      break;

  }

  return A;
}

//
// R^N 2nd-order identity tensor, à la Matlab
//
template<typename T>
inline
Tensor<T> const
eye(Index const N)
{
  return identity<T>(N);
}

//
// R^N 2nd-order tensor transpose
//
template<typename T>
inline
Tensor<T>
transpose(Tensor<T> const & A)
{
  Index const
  N = A.get_dimension();

  Tensor<T>
  B = A;

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = i + 1; j < N; ++j) {
          std::swap(B(i, j), B(j, i));
        }
      }
      break;

    case 3:
      std::swap(B(0,1), B(1,0));
      std::swap(B(0,2), B(2,0));

      std::swap(B(1,2), B(2,1));

      break;

    case 2:
      std::swap(B(0,1), B(1,0));

      break;

  }

  return B;
}

//
// R^N symmetric part of 2nd-order tensor
// \return \f$ \frac{1}{2}(A + A^T) \f$
//
template<typename T>
inline
Tensor<T>
sym(Tensor<T> const & A)
{
  Index const
  N = A.get_dimension();

  Tensor<T>
  B(N);

  switch (N) {

    default:
      B = 0.5 * (A + transpose(A));
      break;

    case 3:
    {
      T const & s00 = A(0,0);
      T const & s11 = A(1,1);
      T const & s22 = A(2,2);

      T const & s01 = 0.5 * (A(0,1) + A(1,0));
      T const & s02 = 0.5 * (A(0,2) + A(2,0));
      T const & s12 = 0.5 * (A(1,2) + A(2,1));

      B(0,0) = s00;
      B(0,1) = s01;
      B(0,2) = s02;

      B(1,0) = s01;
      B(1,1) = s11;
      B(1,2) = s12;

      B(2,0) = s02;
      B(2,1) = s12;
      B(2,2) = s22;
    }
    break;

    case 2:
    {
      T const & s00 = A(0,0);
      T const & s11 = A(1,1);

      T const & s01 = 0.5 * (A(0,1) + A(1,0));

      B(0,0) = s00;
      B(0,1) = s01;

      B(1,0) = s01;
      B(1,1) = s11;
    }
    break;

  }

  return B;
}

//
// R^N skew symmetric part of 2nd-order tensor
// \return \f$ \frac{1}{2}(A - A^T) \f$
//
template<typename T>
inline
Tensor<T>
skew(Tensor<T> const & A)
{
  Index const
  N = A.get_dimension();

  Tensor<T>
  B(N);

  switch (N) {

    default:
      B = 0.5 * (A - transpose(A));
      break;

    case 3:
    {
      T const & s01 = 0.5*(A(0,1)-A(1,0));
      T const & s02 = 0.5*(A(0,2)-A(2,0));
      T const & s12 = 0.5*(A(1,2)-A(2,1));

      B(0,0) = 0.0;
      B(0,1) = s01;
      B(0,2) = s02;

      B(1,0) = -s01;
      B(1,1) = 0.0;
      B(1,2) = s12;

      B(2,0) = -s02;
      B(2,1) = -s12;
      B(2,2) = 0.0;
    }
    break;

    case 2:
    {
      T const & s01 = 0.5*(A(0,1)-A(1,0));

      B(0,0) = 0.0;
      B(0,1) = s01;

      B(1,0) = -s01;
      B(1,1) = 0.0;
    }
    break;

  }

  return B;
}

//
// R^N skew symmetric 2nd-order tensor from vector, undefined
// for N!=3.
// \param u vector
//
template<typename T>
inline
Tensor<T>
skew(Vector<T> const & u)
{
  Index const
  N = u.get_dimension();

  Tensor<T>
  A(N);

  switch (N) {

    case 3:
      A(0,0) = 0.0;
      A(0,1) = -u(2);
      A(0,2) = u(1);

      A(1,0) = u(2);
      A(1,1) = 0.0;
      A(1,2) = -u(0);

      A(2,0) = -u(1);
      A(2,1) = u(0);
      A(2,2) = 0.0;
      break;

    default:
      std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
      std::cerr << std::endl;
      std::cerr << "Skew from vector undefined for R^" << N;
      std::cerr << std::endl;
      exit(1);
      break;

  }

  return A;
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Tensor_i_h
