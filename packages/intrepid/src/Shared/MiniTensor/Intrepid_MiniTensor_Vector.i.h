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

#if !defined(Intrepid_MiniTensor_Vector_i_h)
#define Intrepid_MiniTensor_Vector_i_h

namespace Intrepid
{

//
// Default constructor
//
template<typename T, Index N>
inline
Vector<T, N>::Vector() :
    TensorBase<T, Store>::TensorBase()
{
  return;
}

//
// Constructor that initializes to NaNs
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(Index const dimension) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

///
/// Create vector from a specified value
///
template<typename T, Index N>
inline
Vector<T, N>::Vector(Index const dimension, ComponentValue const value) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

template<typename T, Index N>
inline
Vector<T, N>::Vector(ComponentValue const value) :
    TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

//
// Create vector from array
//
#ifdef HAVE_INTREPID_KOKKOSCORE
template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
    typename not_index<ArrayT>::type & data,
    iType index1) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1)
{
  return;
}

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
    typename not_index<ArrayT>::type & data,
    iType index1,
    iType index2) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1, index2)
{
  return;
}

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
    typename not_index<ArrayT>::type & data,
    iType index1,
    iType index2,
    iType index3) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data, index1, index2, index3)
{
  return;
}

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
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

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
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

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
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

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
    Index const dimension,
    typename not_index<ArrayT>::type & data,
    iType index1) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data, index1)
{
  return;
}

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
    Index const dimension,
    typename not_index<ArrayT>::type & data,
    iType index1,
    iType index2) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data, index1, index2)
{
  return;
}

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
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

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
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

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
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

template<typename T, Index N>
template<class ArrayT, typename iType>
inline
Vector<T, N>::Vector(
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
#endif

template<typename T, Index N>
inline
Vector<T, N>::Vector(Index const dimension, T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N>
inline
Vector<T, N>::Vector(T const * data_ptr) :
    TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(Vector<T, N> const & v) :
    TensorBase<T, Store>::TensorBase(v)
{
  return;
}

//
// Create vector specifying components
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(T const & s0, T const & s1)
{
  Vector < T, N > &
      self = (*this);

  self.set_dimension(2);

  self[0] = s0;
  self[1] = s1;

  return;
}

//
// Create vector specifying components
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(T const & s0, T const & s1, T const & s2)
{
  Vector < T, N > &
      self = (*this);

  self.set_dimension(3);

  self[0] = s0;
  self[1] = s1;
  self[2] = s2;

  return;
}

//
// Simple destructor
//
template<typename T, Index N>
inline
Vector<T, N>::~Vector()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N>
inline Index
Vector<T, N>::get_dimension() const
{
  return IS_DYNAMIC == true ? TensorBase<T, Store>::get_dimension() : N;
}

//
// Set dimension
//
template<typename T, Index N>
inline
void
Vector<T, N>::set_dimension(Index const dimension)
{
  if (IS_DYNAMIC == true) {
    TensorBase<T, Store>::set_dimension(dimension, ORDER);
  }
  else {
    assert(dimension == N);
  }

  return;
}

//
// Indexing for constant vector
//
template<typename T, Index N>
inline T const &
Vector<T, N>::operator()(Index const i) const
    {
  return (*this)[i];
}

//
// Vector indexing
//
template<typename T, Index N>
inline T &
Vector<T, N>::operator()(Index const i)
{
  return (*this)[i];
}

//
// Vector addition
//
template<typename S, typename T, Index N>
inline Vector<typename Promote<S, T>::type, N>
operator+(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Vector<typename Promote<S, T>::type, N>
  w(u.get_dimension());

  add(u, v, w);

  return w;
}

//
// Vector subtraction
//
template<typename S, typename T, Index N>
inline Vector<typename Promote<S, T>::type, N>
operator-(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Vector<typename Promote<S, T>::type, N>
  w(u.get_dimension());

  subtract(u, v, w);

  return w;
}

//
// Vector minus
//
template<typename T, Index N>
inline Vector<T, N>
operator-(Vector<T, N> const & u)
{
  Vector<T, N>
  v(u.get_dimension());

  minus(u, v);

  return v;
}

//
// Vector dot product
//
template<typename S, typename T, Index N>
inline typename Promote<S, T>::type
operator*(Vector<S, N> const & u, Vector<T, N> const & v)
{
  return dot(u, v);
}

//
// Vector equality tested by components
//
template<typename T, Index N>
inline
bool
operator==(Vector<T, N> const & u, Vector<T, N> const & v)
{
  return equal(u, v);
}

//
// Vector inequality tested by components
//
template<typename T, Index N>
inline
bool
operator!=(Vector<T, N> const & u, Vector<T, N> const & v)
{
  return not_equal(u, v);
}

//
// Scalar vector product
//
template<typename S, typename T, Index N>
inline typename
lazy_disable_if<order_1234<S>, apply_vector<Promote<S, T>, N> >::type
operator*(S const & s, Vector<T, N> const & u)
{
  Vector<typename Promote<S, T>::type, N>
  v(u.get_dimension());

  scale(u, s, v);

  return v;
}

//
// Vector scalar product
//
template<typename S, typename T, Index N>
inline typename
lazy_disable_if<order_1234<S>, apply_vector<Promote<S, T>, N> >::type
operator*(Vector<T, N> const & u, S const & s)
{
  Vector<typename Promote<S, T>::type, N>
  v(u.get_dimension());

  scale(u, s, v);

  return v;
}

//
// Vector scalar division
//
template<typename S, typename T, Index N>
inline Vector<typename Promote<S, T>::type, N>
operator/(Vector<T, N> const & u, S const & s)
{
  Vector<typename Promote<S, T>::type, N>
  v(u.get_dimension());

  divide(u, s, v);

  return v;
}

//
// Scalar vector division
//
template<typename S, typename T, Index N>
inline Vector<typename Promote<S, T>::type, N>
operator/(S const & s, Vector<T, N> const & u)
{
  Vector<typename Promote<S, T>::type, N>
  v(u.get_dimension());

  split(u, s, v);

  return v;
}

//
// Vector dot product
//
template<typename S, typename T, Index N>
inline typename Promote<S, T>::type
dot(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Index const
  dimension = u.get_dimension();

  assert(v.get_dimension() == dimension);

  typename Promote<S, T>::type
  s = 0.0;

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      s += u(i) * v(i);
    }
    break;

  case 3:
    s = u(0) * v(0) + u(1) * v(1) + u(2) * v(2);
    break;

  case 2:
    s = u(0) * v(0) + u(1) * v(1);
    break;

  }

  return s;
}

//
// Cross product only valid for R^3.
//
template<typename S, typename T, Index N>
inline Vector<typename Promote<S, T>::type, N>
cross(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Index const
  dimension = u.get_dimension();

  assert(v.get_dimension() == dimension);

  Vector<typename Promote<S, T>::type, N>
  w(dimension);

  switch (dimension) {

  case 3:
    w(0) = u(1) * v(2) - u(2) * v(1);
    w(1) = u(2) * v(0) - u(0) * v(2);
    w(2) = u(0) * v(1) - u(1) * v(0);
    break;

  default:
    std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
    std::cerr << std::endl;
    std::cerr << "Cross product undefined for R^" << dimension;
    std::cerr << std::endl;
    exit(1);
    break;

  }

  return w;
}

//
// R^N vector 2-norm
// \return \f$ \sqrt{u \cdot u} \f$
//
template<typename T, Index N>
inline T
norm(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

  default:
    s = std::sqrt(dot(u, u));
    break;

  case 3:
    s = std::sqrt(u(0) * u(0) + u(1) * u(1) + u(2) * u(2));
    break;

  case 2:
    s = std::sqrt(u(0) * u(0) + u(1) * u(1));
    break;

  }

  return s;
}

//
// R^N vector 2-norm square for fast distance calculations.
// \return \f$ u \cdot u \f$
//
template<typename T, Index N>
inline T
norm_square(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

  default:
    s = dot(u, u);
    break;

  case 3:
    s = u(0) * u(0) + u(1) * u(1) + u(2) * u(2);
    break;

  case 2:
    s = u(0) * u(0) + u(1) * u(1);
    break;

  }

  return s;
}

//
// R^N vector 1-norm
// \return \f$ \sum_i |u_i| \f$
//
template<typename T, Index N>
inline T
norm_1(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      s += std::abs(u(i));
    }
    break;

  case 3:
    s = std::abs(u(0)) + std::abs(u(1)) + std::abs(u(2));
    break;

  case 2:
    s = std::abs(u(0)) + std::abs(u(1));
    break;

  }

  return s;
}

//
// R^N vector infinity-norm
// \return \f$ \max(|u_0|,...|u_i|,...|u_N|) \f$
//
template<typename T, Index N>
inline T
norm_infinity(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      s = std::max(s, std::abs(u(i)));
    }
    break;

  case 3:
    s = std::max(std::max(std::abs(u(0)), std::abs(u(1))), std::abs(u(2)));
    break;

  case 2:
    s = std::max(std::abs(u(0)), std::abs(u(1)));
    break;

  }

  return s;
}

//
// \return u / |u|, fails for |u| = 0
//
template<typename T, Index N>
inline Vector<T, N>
unit(Vector<T, N> const & u)
{
  return u / norm(u);
}

//
// Compute Householder vector
//
template<typename T, Index N>
std::pair<Vector<T, N>, T>
house(Vector<T, N> const & x)
{
  Vector<T, N>
  v = x;

  v[0] = 1.0;

  Index const
  dimension = x.get_dimension();

  T
  sigma = 0.0;

  for (Index i = 1; i < dimension; ++i) {
    sigma = v[i] * v[i];
  }

  T
  beta = 0.0;

  if (sigma > 0.0) {
    T const
    mu = std::sqrt(x[0] * x[0] + sigma);
    v[0] = x[0] > 0.0 ? -sigma / (x[0] + mu) : x[0] - mu;
    beta = 2.0 * v[0] * v[0] / (sigma + v[0] * v[0]);
    v = v / v[0];
  }

  return std::make_pair(v, beta);
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Vector_i_h
