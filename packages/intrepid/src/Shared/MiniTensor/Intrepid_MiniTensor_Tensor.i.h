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
// and eigenvalues and corresponding vectors in the
// respective decompositions.
// \param u vector to sort
// \return pair<v, P>
// \return v sorted vector
// \return P permutation matrix such that v = P^T u
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
// get dimension
//
template<typename T>
inline
Index
Tensor<T>::get_dimension() const
{
  return dimension;
}

//
// set dimension
//
template<typename T>
inline
void
Tensor<T>::set_dimension(Index const N)
{
  if (N == get_dimension()) return;

  Index const
  number_components = N * N;

  e.resize(number_components);

  dimension = N;

  return;
}

//
// Fill components from array defined by pointer.
// \param data_ptr pointer into array for filling components
// \param order component convention (3D only)
//
template<typename T>
inline
void
Tensor<T>::fill(T const * data_ptr, ComponentOrder const order)
{
  assert(data_ptr != NULL);

  Index const
  N = get_dimension();

  Index const
  number_components = N * N;

  switch (N) {

    default:
      for (Index i = 0; i < number_components; ++i) {
        e[i] = data_ptr[i];
      }
      break;

    case 3:

      switch (order) {

        default:
          std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
          std::cerr << std::endl;
          std::cerr << "Inverse of singular tensor.";
          std::cerr << std::endl;
          exit(1);
          break;

        case CANONICAL:
          e[0] = data_ptr[0];
          e[1] = data_ptr[1];
          e[2] = data_ptr[2];

          e[3] = data_ptr[3];
          e[4] = data_ptr[4];
          e[5] = data_ptr[5];

          e[6] = data_ptr[6];
          e[7] = data_ptr[7];
          e[8] = data_ptr[8];
          break;

        case SIERRA_FULL:

          //  0  1  2  3  4  5  6  7  8
          // XX YY ZZ XY YZ ZX YX ZY XZ
          //  0  4  8  1  5  6  3  7  2
          e[0] = data_ptr[0];
          e[4] = data_ptr[1];
          e[8] = data_ptr[2];

          e[1] = data_ptr[3];
          e[5] = data_ptr[4];
          e[6] = data_ptr[5];

          e[3] = data_ptr[6];
          e[7] = data_ptr[7];
          e[2] = data_ptr[8];
          break;

        case SIERRA_SYMMETRIC:
          e[0] = data_ptr[0];
          e[4] = data_ptr[1];
          e[8] = data_ptr[2];

          e[1] = data_ptr[3];
          e[5] = data_ptr[4];
          e[6] = data_ptr[5];

          e[3] = data_ptr[3];
          e[7] = data_ptr[4];
          e[2] = data_ptr[5];
          break;
      }

      break;

    case 2:
      e[0] = data_ptr[0];
      e[1] = data_ptr[1];

      e[2] = data_ptr[2];
      e[3] = data_ptr[3];
      break;

  }

  return;
}

//
// default constructor
//
template<typename T>
inline
Tensor<T>::Tensor() :
dimension(0)
{
  return;
}

//
// constructor that initializes to NaNs
//
template<typename T>
inline
Tensor<T>::Tensor(Index const N) :
dimension(0)
{
  set_dimension(N);

  Index const
  number_components = N * N;

  switch (N) {

    default:
      for (Index i = 0; i < number_components; ++i) {
        e[i] = not_a_number<T>();
      }
      break;

    case 3:
      e[0] = not_a_number<T>();
      e[1] = not_a_number<T>();
      e[2] = not_a_number<T>();

      e[3] = not_a_number<T>();
      e[4] = not_a_number<T>();
      e[5] = not_a_number<T>();

      e[6] = not_a_number<T>();
      e[7] = not_a_number<T>();
      e[8] = not_a_number<T>();
      break;

    case 2:
      e[0] = not_a_number<T>();
      e[1] = not_a_number<T>();

      e[2] = not_a_number<T>();
      e[3] = not_a_number<T>();
      break;

  }
  return;
}

//
// Create tensor from a scalar
// \param s all components are set equal to this value
//
template<typename T>
inline
Tensor<T>::Tensor(Index const N, T const & s) :
dimension(0)
{
  set_dimension(N);

  Index const
  number_components = N * N;

  switch (N) {

    default:
      for (Index i = 0; i < number_components; ++i) {
        e[i] = s;
      }
      break;

    case 3:
      e[0] = s;
      e[1] = s;
      e[2] = s;

      e[3] = s;
      e[4] = s;
      e[5] = s;

      e[6] = s;
      e[7] = s;
      e[8] = s;
      break;

    case 2:
      e[0] = s;
      e[1] = s;

      e[2] = s;
      e[3] = s;
      break;

  }

  return;
}

//
// Create tensor specifying components
// \param N dimension
// \param  s00, s01, ... components in the R^2 canonical basis
//
template<typename T>
inline
Tensor<T>::Tensor(
    T const & s00, T const & s01,
    T const & s10, T const & s11) :
dimension(0)
{
  set_dimension(2);

  e[0] = s00;
  e[1] = s01;

  e[2] = s10;
  e[3] = s11;

  return;
}

//
// Create tensor specifying components
// \param N dimension
// \param  s00, s01, ... components in the R^3 canonical basis
//
template<typename T>
inline
Tensor<T>::Tensor(
    T const & s00, T const & s01, T const & s02,
    T const & s10, T const & s11, T const & s12,
    T const & s20, T const & s21, T const & s22) :
    dimension(0)
{
  set_dimension(3);

  e[0] = s00;
  e[1] = s01;
  e[2] = s02;

  e[3] = s10;
  e[4] = s11;
  e[5] = s12;

  e[6] = s20;
  e[7] = s21;
  e[8] = s22;

  return;
}

//
// R^N create tensor from array - const version
// \param data_ptr
// \param order component convention (3D only)
//
template<typename T>
inline
Tensor<T>::Tensor(
    Index const N,
    T const * data_ptr,
    ComponentOrder const order) :
dimension(0)
{
  assert(data_ptr != NULL);

  set_dimension(N);

  fill(data_ptr, order);

  return;
}

//
// R^N copy constructor
// \param A the values of its components are copied to the new tensor
//
template<typename T>
inline
Tensor<T>::Tensor(Tensor<T> const & A) :
dimension(0)
{
  Index const
  N = A.get_dimension();

  set_dimension(N);

  Index const
  number_components = N * N;

  switch (N) {

    default:
      for (Index i = 0; i < number_components; ++i) {
        e[i] = A.e[i];
      }
      break;

    case 3:
      e[0] = A.e[0];
      e[1] = A.e[1];
      e[2] = A.e[2];

      e[3] = A.e[3];
      e[4] = A.e[4];
      e[5] = A.e[5];

      e[6] = A.e[6];
      e[7] = A.e[7];
      e[8] = A.e[8];
      break;

    case 2:
      e[0] = A.e[0];
      e[1] = A.e[1];

      e[2] = A.e[2];
      e[3] = A.e[3];
      break;

  }

  return;
}

//
// 2nd-order tensor from 4th-order tensor
//
template<typename T>
inline
Tensor<T>::Tensor(Tensor4<T> const & A) :
dimension(0)
{
  Index const
  N = A.get_dimension();

  Index const
  N2 = N * N;

  set_dimension(N2);

  Index const
  number_components = N2 * N2;

  for (Index i = 0; i < number_components; ++i) {
    e[i] = A[i];
  }

  return;
}

//
// R^N simple destructor
//
template<typename T>
inline
Tensor<T>::~Tensor()
{
  return;
}

//
// R^N indexing for constant tensor
// \param i index
// \param j index
//
template<typename T>
inline
T const &
Tensor<T>::operator()(Index const i, Index const j) const
{
  assert(i < get_dimension());
  assert(j < get_dimension());
  return e[i * get_dimension() + j];
}

//
// R^N tensor indexing
// \param i index
// \param j index
//
template<typename T>
inline
T &
Tensor<T>::operator()(Index const i, Index const j)
{
  assert(i < get_dimension());
  assert(j < get_dimension());
  return e[i * get_dimension() + j];
}

//
// Linear access to components
// \param i the index
//
template<typename T>
inline
T const &
Tensor<T>::operator[](Index const i) const
{
  assert(i < integer_power(get_dimension(), order()));
  return e[i];
}

//
// Linear access to components
// \param i the index
//
template<typename T>
inline
T &
Tensor<T>::operator[](Index const i)
{
  assert(i < integer_power(get_dimension(), order()));
  return e[i];
}

//
// R^N copy assignment
// \param A the values of its components are copied to this tensor
//
template<typename T>
inline
Tensor<T> &
Tensor<T>::operator=(Tensor<T> const & A)
{
  if (this == &A) return *this;

  Index const
  N = A.get_dimension();

  set_dimension(N);

  Index const
  number_components = N * N;

  switch (N) {

    default:
      for (Index i = 0; i < number_components; ++i) {
        e[i] = A.e[i];
      }
      break;

    case 3:
      e[0] = A.e[0];
      e[1] = A.e[1];
      e[2] = A.e[2];

      e[3] = A.e[3];
      e[4] = A.e[4];
      e[5] = A.e[5];

      e[6] = A.e[6];
      e[7] = A.e[7];
      e[8] = A.e[8];
      break;

    case 2:
      e[0] = A.e[0];
      e[1] = A.e[1];

      e[2] = A.e[2];
      e[3] = A.e[3];
      break;

  }

  return *this;
}

//
// R^N tensor increment
// \param A added to current tensor
//
template<typename T>
inline
Tensor<T> &
Tensor<T>::operator+=(Tensor<T> const & A)
{
  Index const
  N = get_dimension();

  assert(A.get_dimension() == N);

  Index const
  number_components = N * N;

  switch (N) {

    default:
      for (Index i = 0; i < number_components; ++i) {
        e[i] += A.e[i];
      }
      break;

    case 3:
      e[0] += A.e[0];
      e[1] += A.e[1];
      e[2] += A.e[2];

      e[3] += A.e[3];
      e[4] += A.e[4];
      e[5] += A.e[5];

      e[6] += A.e[6];
      e[7] += A.e[7];
      e[8] += A.e[8];
      break;

    case 2:
      e[0] += A.e[0];
      e[1] += A.e[1];

      e[2] += A.e[2];
      e[3] += A.e[3];
      break;

  }

  return *this;
}

//
// R^N tensor decrement
// \param A substracted from current tensor
//
template<typename T>
inline
Tensor<T> &
Tensor<T>::operator-=(Tensor<T> const & A)
{
  Index const
  N = get_dimension();

  assert(A.get_dimension() == N);

  Index const
  number_components = N * N;

  switch (N) {

    default:
      for (Index i = 0; i < number_components; ++i) {
        e[i] -= A.e[i];
      }
      break;

    case 3:
      e[0] -= A.e[0];
      e[1] -= A.e[1];
      e[2] -= A.e[2];

      e[3] -= A.e[3];
      e[4] -= A.e[4];
      e[5] -= A.e[5];

      e[6] -= A.e[6];
      e[7] -= A.e[7];
      e[8] -= A.e[8];
      break;

    case 2:
      e[0] -= A.e[0];
      e[1] -= A.e[1];

      e[2] -= A.e[2];
      e[3] -= A.e[3];
      break;

  }


  return *this;
}

//
// R^N fill with zeros
//
template<typename T>
inline
void
Tensor<T>::clear()
{
  Index const
  N = get_dimension();

  Index const
  number_components = N * N;

  switch (N) {

    default:
      for (Index i = 0; i < number_components; ++i) {
        e[i] = 0.0;
      }
      break;

    case 3:
      e[0] = 0.0;
      e[1] = 0.0;
      e[2] = 0.0;

      e[3] = 0.0;
      e[4] = 0.0;
      e[5] = 0.0;

      e[6] = 0.0;
      e[7] = 0.0;
      e[8] = 0.0;
      break;

    case 2:
      e[0] = 0.0;
      e[1] = 0.0;

      e[2] = 0.0;
      e[3] = 0.0;
      break;

  }


  return;
}

//
// R^N tensor addition
// \return \f$ A + B \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
operator+(Tensor<S> const & A, Tensor<T> const & B)
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
          C(i, j) = A(i, j) + B(i, j);
        }
      }
      break;

    case 3:
      C(0,0) = A(0,0) + B(0,0);
      C(0,1) = A(0,1) + B(0,1);
      C(0,2) = A(0,2) + B(0,2);

      C(1,0) = A(1,0) + B(1,0);
      C(1,1) = A(1,1) + B(1,1);
      C(1,2) = A(1,2) + B(1,2);

      C(2,0) = A(2,0) + B(2,0);
      C(2,1) = A(2,1) + B(2,1);
      C(2,2) = A(2,2) + B(2,2);
      break;

    case 2:
      C(0,0) = A(0,0) + B(0,0);
      C(0,1) = A(0,1) + B(0,1);

      C(1,0) = A(1,0) + B(1,0);
      C(1,1) = A(1,1) + B(1,1);
      break;

  }

  return C;
}

//
// R^N Tensor substraction
// \return \f$ A - B \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
operator-(Tensor<S> const & A, Tensor<T> const & B)
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
          C(i, j) = A(i, j) - B(i, j);
        }
      }
      break;

    case 3:
      C(0,0) = A(0,0) - B(0,0);
      C(0,1) = A(0,1) - B(0,1);
      C(0,2) = A(0,2) - B(0,2);

      C(1,0) = A(1,0) - B(1,0);
      C(1,1) = A(1,1) - B(1,1);
      C(1,2) = A(1,2) - B(1,2);

      C(2,0) = A(2,0) - B(2,0);
      C(2,1) = A(2,1) - B(2,1);
      C(2,2) = A(2,2) - B(2,2);
      break;

    case 2:
      C(0,0) = A(0,0) - B(0,0);
      C(0,1) = A(0,1) - B(0,1);

      C(1,0) = A(1,0) - B(1,0);
      C(1,1) = A(1,1) - B(1,1);
      break;

  }

  return C;
}

//
// R^N tensor minus
// \return \f$ -A \f$
//
template<typename T>
inline
Tensor<T>
operator-(Tensor<T> const & A)
{
  Index const
  N = A.get_dimension();

  Tensor<T> S(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {
          S(i, j) = - A(i, j);
        }
      }
      break;

    case 3:
      S(0,0) = -A(0,0);
      S(0,1) = -A(0,1);
      S(0,2) = -A(0,2);

      S(1,0) = -A(1,0);
      S(1,1) = -A(1,1);
      S(1,2) = -A(1,2);

      S(2,0) = -A(2,0);
      S(2,1) = -A(2,1);
      S(2,2) = -A(2,2);
      break;

    case 2:
      S(0,0) = -A(0,0);
      S(0,1) = -A(0,1);

      S(1,0) = -A(1,0);
      S(1,1) = -A(1,1);
      break;

  }

  return S;
}

//
// R^N tensor equality
// Tested by components
// \return \f$ A \equiv B \f$
//
template<typename T>
inline
bool
operator==(Tensor<T> const & A, Tensor<T> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {
          if (A(i, j) != B(i, j)) {
            return false;
          }
        }
      }
      break;

    case 3:
      return
          A(0,0)==B(0,0) && A(0,1)==B(0,1) && A(0,2)==B(0,2) &&
          A(1,0)==B(1,0) && A(1,1)==B(1,1) && A(1,2)==B(1,2) &&
          A(2,0)==B(2,0) && A(2,1)==B(2,1) && A(2,2)==B(2,2);
      break;

    case 2:
      return
          A(0,0)==B(0,0) && A(0,1)==B(0,1) &&
          A(1,0)==B(1,0) && A(1,1)==B(1,1);
      break;

  }

  return true;
}

//
// R^N tensor inequality
// Tested by components
// \return \f$ A \neq B \f$
//
template<typename T>
inline
bool
operator!=(Tensor<T> const & A, Tensor<T> const & B)
{
  return !(A == B);
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
operator*(Tensor<T> const & A, Vector<S> const & u)
{
  return dot(A, u);
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
operator*(Vector<S> const & u, Tensor<T> const & A)
{
  return dot(u, A);
}

//
// R^N tensor dot product C = A B
// \return \f$ A \cdot B \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
operator*(Tensor<S> const & A, Tensor<T> const & B)
{
  return dot(A, B);
}

//
// R^N scalar tensor product
// \param s scalar
// \param A tensor
// \return \f$ s A \f$
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_tensor< Promote<S,T> > >::type
operator*(S const & s, Tensor<T> const & A)
{
  Index const
  N = A.get_dimension();

  Tensor<typename Promote<S, T>::type>
  B(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {
          B(i, j) = s * A(i, j);
        }
      }
      break;

    case 3:
      B(0,0) = s * A(0,0);
      B(0,1) = s * A(0,1);
      B(0,2) = s * A(0,2);

      B(1,0) = s * A(1,0);
      B(1,1) = s * A(1,1);
      B(1,2) = s * A(1,2);

      B(2,0) = s * A(2,0);
      B(2,1) = s * A(2,1);
      B(2,2) = s * A(2,2);
      break;

    case 2:
      B(0,0) = s * A(0,0);
      B(0,1) = s * A(0,1);

      B(1,0) = s * A(1,0);
      B(1,1) = s * A(1,1);
      break;

  }

  return B;
}

//
// R^N tensor scalar product
// \param A tensor
// \param s scalar
// \return \f$ s A \f$
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_tensor< Promote<S,T> > >::type
operator*(Tensor<T> const & A, S const & s)
{
  return s * A;
}

//
// R^N tensor scalar division
// \param A tensor
// \param s scalar
// \return \f$ A / s \f$
//
template<typename S, typename T>
inline
Tensor<typename Promote<S, T>::type>
operator/(Tensor<T> const & A, S const & s)
{
  Index const
  N = A.get_dimension();

  Tensor<typename Promote<S, T>::type>
  B(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {
          B(i, j) = A(i, j) / s;
        }
      }
      break;

    case 3:
      B(0,0) = A(0,0) / s;
      B(0,1) = A(0,1) / s;
      B(0,2) = A(0,2) / s;

      B(1,0) = A(1,0) / s;
      B(1,1) = A(1,1) / s;
      B(1,2) = A(1,2) / s;

      B(2,0) = A(2,0) / s;
      B(2,1) = A(2,1) / s;
      B(2,2) = A(2,2) / s;
      break;

    case 2:
      B(0,0) = A(0,0) / s;
      B(0,1) = A(0,1) / s;

      B(1,0) = A(1,0) / s;
      B(1,1) = A(1,1) / s;
      break;

  }

  return B;
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
// R^N zero 2nd-order tensor
// All components are zero
//
template<typename T>
inline
Tensor<T> const
zero(Index const N)
{
  return Tensor<T>(N, 0.0);
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
  A(N, 0.0);

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

    default:
      std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
      std::cerr << std::endl;
      std::cerr << "Skew from vector undefined for R^" << N;
      std::cerr << std::endl;
      exit(1);
      break;

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

  }

  return A;
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Tensor_i_h
