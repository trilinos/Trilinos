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
template<typename T>
inline
Vector<T>::Vector() :
TensorBase<T>::TensorBase()
{
  return;
}

//
// Constructor that initializes to NaNs
//
template<typename T>
inline
Vector<T>::Vector(Index const dimension) :
TensorBase<T>::TensorBase(dimension, order)
{
  return;
}

//
// Create vector from a scalar
//
template<typename T>
inline
Vector<T>::Vector(Index const dimension, T const & s) :
TensorBase<T>::TensorBase(dimension, order, s)
{
  return;
}

//
// Create vector specifying components
//
template<typename T>
inline
Vector<T>::Vector(T const & s0, T const & s1)
{
  Vector<T> &
  self = (*this);

  self.set_dimension(2);

  self[0] = s0;
  self[1] = s1;

  return;
}

//
// Create vector specifying components
//
template<typename T>
inline
Vector<T>::Vector(T const & s0, T const & s1, T const & s2)
{
  Vector<T> &
  self = (*this);

  self.set_dimension(3);

  self[0] = s0;
  self[1] = s1;
  self[2] = s2;

  return;
}

//
// Create vector from array
//
template<typename T>
inline
Vector<T>::Vector(Index const dimension, T const * data_ptr) :
TensorBase<T>::TensorBase(dimension, order, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T>
inline
Vector<T>::Vector(Vector<T> const & v) :
TensorBase<T>::TensorBase(v)
{
  return;
}

//
// Simple destructor
//
template<typename T>
inline
Vector<T>::~Vector()
{
  return;
}

//
// Indexing for constant vector
//
template<typename T>
inline
T const &
Vector<T>::operator()(Index const i) const
{
  return (*this)[i];
}

//
// Vector indexing
//
template<typename T>
inline
T &
Vector<T>::operator()(Index const i)
{
  return (*this)[i];
}

//
// Vector addition
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
operator+(Vector<S> const & u, Vector<T> const & v)
{
  Vector<typename Promote<S, T>::type>
  w;

  add(u, v, w);

  return w;
}

//
// Vector subtraction
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
operator-(Vector<S> const & u, Vector<T> const & v)
{
  Vector<typename Promote<S, T>::type>
  w;

  subtract(u, v, w);

  return w;
}

//
// Vector minus
//
template<typename T>
inline
Vector<T>
operator-(Vector<T> const & u)
{
  Vector<T>
  v;

  minus(u, v);

  return v;
}

//
// Vector dot product
//
template<typename S, typename T>
inline
typename Promote<S, T>::type
operator*(Vector<S> const & u, Vector<T> const & v)
{
  return dot(u, v);
}

//
// Vector equality tested by components
//
template<typename T>
inline
bool
operator==(Vector<T> const & u, Vector<T> const & v)
{
  return equal(u, v);
}

//
// Vector inequality tested by components
//
template<typename T>
inline
bool
operator!=(Vector<T> const & u, Vector<T> const & v)
{
  return not_equal(u, v);
}

//
// Scalar vector product
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_vector< Promote<S,T> > >::type
operator*(S const & s, Vector<T> const & u)
{
  Vector<typename Promote<S, T>::type>
  v;

  scale(u, s, v);

  return v;
}

//
// Vector scalar product
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_vector< Promote<S,T> > >::type
operator*(Vector<T> const & u, S const & s)
{
  Vector<typename Promote<S, T>::type>
  v;

  scale(u, s, v);

  return v;
}

//
// Vector scalar division
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
operator/(Vector<T> const & u, S const & s)
{
  Vector<typename Promote<S, T>::type>
  v;

  divide(u, s, v);

  return v;
}

//
// Vector dot product
//
template<typename S, typename T>
inline
typename Promote<S, T>::type
dot(Vector<S> const & u, Vector<T> const & v)
{
  Index const
  N = u.get_dimension();

  assert(v.get_dimension() == N);

  typename Promote<S, T>::type
  s = 0.0;

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
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
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
cross(Vector<S> const & u, Vector<T> const & v)
{
  Index const
  N = u.get_dimension();

  assert(v.get_dimension() == N);

  Vector<typename Promote<S, T>::type>
  w(N);

  switch (N) {

    default:
      std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
      std::cerr << std::endl;
      std::cerr << "Cross product undefined for R^" << N;
      std::cerr << std::endl;
      exit(1);
      break;

    case 3:
      w(0) = u(1) * v(2) - u(2) * v(1);
      w(1) = u(2) * v(0) - u(0) * v(2);
      w(2) = u(0) * v(1) - u(1) * v(0);
      break;

  }

  return w;
}

//
// R^N vector 2-norm
// \return \f$ \sqrt{u \cdot u} \f$
//
template<typename T>
inline
T
norm(Vector<T> const & u)
{
  Index const
  N = u.get_dimension();

  T
  s = 0.0;

  switch (N) {

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
template<typename T>
inline
T
norm_square(Vector<T> const & u)
{
  Index const
  N = u.get_dimension();

  T
  s = 0.0;

  switch (N) {

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
template<typename T>
inline
T
norm_1(Vector<T> const & u)
{
  Index const
  N = u.get_dimension();

  T
  s = 0.0;

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
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
template<typename T>
inline
T
norm_infinity(Vector<T> const & u)
{
  Index const
  N = u.get_dimension();

  T
  s = 0.0;

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
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

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Vector_i_h
