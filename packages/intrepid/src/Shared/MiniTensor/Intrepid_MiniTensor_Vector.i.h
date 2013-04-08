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
// return dimension
//
template<typename T>
inline
Index
Vector<T>::get_dimension() const
{
  return dimension;
}

//
// set dimension
//
template<typename T>
inline
void
Vector<T>::set_dimension(Index const N)
{
  if (N == get_dimension()) return;

  e.resize(N);

  dimension = N;

  return;
}

//
// Fill components from array defined by pointer.
// \param data_ptr pointer into array for filling components
//
template<typename T>
inline
void
Vector<T>::fill(T const * data_ptr)
{
  assert(data_ptr != NULL);

  Index const
  N = get_dimension();

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        e[i] = data_ptr[i];
      }
      break;

    case 3:
      e[0] = data_ptr[0];
      e[1] = data_ptr[1];
      e[2] = data_ptr[2];
      break;

    case 2:
      e[0] = data_ptr[0];
      e[1] = data_ptr[1];
      break;

  }

  return;
}

//
// default constructor
//
template<typename T>
inline
Vector<T>::Vector() :
dimension(0),
e(Teuchos::null)
{
  return;
}

//
// constructor that initializes to NaNs
//
template<typename T>
inline
Vector<T>::Vector(Index const N) :
dimension(0),
e(Teuchos::null)
{
  set_dimension(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        e[i] = not_a_number<T>();
      }
      break;

    case 3:
      e[0] = not_a_number<T>();
      e[1] = not_a_number<T>();
      e[2] = not_a_number<T>();
      break;

    case 2:
      e[0] = not_a_number<T>();
      e[1] = not_a_number<T>();
      break;

  }

  return;
}

//
// Create vector from a scalar
// \param s all components are set equal to this value
//
template<typename T>
inline
Vector<T>::Vector(Index const N, T const & s) :
dimension(0),
e(Teuchos::null)
{
  set_dimension(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        e[i] = s;
      }
      break;

    case 3:
      e[0] = s;
      e[1] = s;
      e[2] = s;
      break;

    case 2:
      e[0] = s;
      e[1] = s;
      break;

  }

  return;
}

//
// Create vector specifying components
// \param N dimension
// \param s0, s1 are the vector components in the R^2 canonical basis
//
template<typename T>
inline
Vector<T>::Vector(T const & s0, T const & s1) :
dimension(0),
e(Teuchos::null)
{
  set_dimension(2);

  e[0] = s0;
  e[1] = s1;

  return;
}

//
// Create vector specifying components
// \param N dimension
// \param s0, s1, s2 are the vector components in the R^3 canonical basis
//
template<typename T>
inline
Vector<T>::Vector(T const & s0, T const & s1, T const & s2) :
dimension(0),
e(Teuchos::null)
{
  set_dimension(3);

  e[0] = s0;
  e[1] = s1;
  e[2] = s2;

  return;
}

//
// R^N create vector from array - const version
// \param N dimension
// \param data_ptr
//
template<typename T>
inline
Vector<T>::Vector(Index const N, T const * data_ptr) :
dimension(0),
e(Teuchos::null)
{
  assert(data_ptr != NULL);

  set_dimension(N);

  fill(data_ptr);

  return;
}

//
// R^N copy constructor
// \param v the values of its components are copied to the new vector
//
template<typename T>
inline
Vector<T>::Vector(Vector<T> const & v) :
dimension(0),
e(Teuchos::null)
{
  Index const
  N = v.get_dimension();

  set_dimension(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        e[i] = v.e[i];
      }
      break;

    case 3:
      e[0] = v.e[0];
      e[1] = v.e[1];
      e[2] = v.e[2];
      break;

    case 2:
      e[0] = v.e[0];
      e[1] = v.e[1];
      break;

  }

  return;
}

//
// R^N simple destructor
//
template<typename T>
inline
Vector<T>::~Vector()
{
  return;
}

//
// R^N indexing for constant vector
// \param i the index
//
template<typename T>
inline
T const &
Vector<T>::operator()(Index const i) const
{
  assert(i < get_dimension());
  return e[i];
}

//
// R^N vector indexing
// \param i the index
//
template<typename T>
inline
T &
Vector<T>::operator()(Index const i)
{
  assert(i < get_dimension());
  return e[i];
}

//
// R^N copy assignment
// \param v the values of its components are copied to this vector
//
template<typename T>
inline
Vector<T> &
Vector<T>::operator=(Vector<T> const & v)
{
  if (this == &v) return *this;

  Index const
  N = v.get_dimension();

  set_dimension(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        e[i] = v.e[i];
      }
      break;

    case 3:
      e[0] = v.e[0];
      e[1] = v.e[1];
      e[2] = v.e[2];
      break;

    case 2:
      e[0] = v.e[0];
      e[1] = v.e[1];
      break;

  }

  return *this;
}

//
// R^N vector increment
// \param v added to currrent vector
//
template<typename T>
inline
Vector<T> &
Vector<T>::operator+=(Vector<T> const & v)
{
  Index const
  N = get_dimension();

  assert(v.get_dimension() == N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        e[i] += v.e[i];
      }
      break;

    case 3:
      e[0] += v.e[0];
      e[1] += v.e[1];
      e[2] += v.e[2];
      break;

    case 2:
      e[0] += v.e[0];
      e[1] += v.e[1];
      break;

  }

  return *this;
}

//
// R^N vector decrement
// \param v substracted from current vector
//
template<typename T>
inline
Vector<T> &
Vector<T>::operator-=(Vector<T> const & v)
{
  Index const
  N = get_dimension();

  assert(v.get_dimension() == N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        e[i] -= v.e[i];
      }
      break;

    case 3:
      e[0] -= v.e[0];
      e[1] -= v.e[1];
      e[2] -= v.e[2];
      break;

    case 2:
      e[0] -= v.e[0];
      e[1] -= v.e[1];
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
Vector<T>::clear()
{
  Index const
  N = get_dimension();

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        e[i] = 0.0;
      }
      break;

    case 3:
      e[0] = 0.0;
      e[1] = 0.0;
      e[2] = 0.0;
      break;

    case 2:
      e[0] = 0.0;
      e[1] = 0.0;
      break;

  }

  return;
}

//
// R^N vector addition
// \param u
// \param v the operands
// \return \f$ u + v \f$
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
operator+(Vector<S> const & u, Vector<T> const & v)
{
  Index const
  N = u.get_dimension();

  assert(v.get_dimension() == N);

  Vector<typename Promote<S, T>::type>
  s(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        s(i) = u(i) + v(i);
      }
      break;

    case 3:
      s(0) = u(0) + v(0);
      s(1) = u(1) + v(1);
      s(2) = u(2) + v(2);
      break;

    case 2:
      s(0) = u(0) + v(0);
      s(1) = u(1) + v(1);
      break;

  }

  return s;
}

//
// R^N vector substraction
// \param u
// \param v the operands
// \return \f$ u - v \f$
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
operator-(Vector<S> const & u, Vector<T> const & v)
{
  Index const
  N = u.get_dimension();

  assert(v.get_dimension() == N);

  Vector<typename Promote<S, T>::type>
  s(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        s(i) = u(i) - v(i);
      }
      break;

    case 3:
      s(0) = u(0) - v(0);
      s(1) = u(1) - v(1);
      s(2) = u(2) - v(2);
      break;

    case 2:
      s(0) = u(0) - v(0);
      s(1) = u(1) - v(1);
      break;

  }

  return s;
}

//
// R^N vector minus
// \param u
// \return \f$ -u \f$
//
template<typename T>
inline
Vector<T>
operator-(Vector<T> const & u)
{
  Index const
  N = u.get_dimension();

  Vector<T>
  v(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        v(i) = -u(i);
      }
      break;

    case 3:
      v(0) = -u(0);
      v(1) = -u(1);
      v(2) = -u(2);
      break;

    case 2:
      v(0) = -u(0);
      v(1) = -u(1);
      break;

  }

  return v;
}

//
// R^N vector dot product
// \param u
// \param v the operands
// \return \f$ u \cdot v \f$
//
template<typename S, typename T>
inline
typename Promote<S, T>::type
operator*(Vector<S> const & u, Vector<T> const & v)
{
  return dot(u, v);
}

//
// R^N vector equality tested by components
// \param u
// \param v the operands
// \return \f$ u \equiv v \f$
//
template<typename T>
inline
bool
operator==(Vector<T> const & u, Vector<T> const & v)
{
  Index const
  N = u.get_dimension();

  assert(v.get_dimension() == N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        if (v(i) != u(i)) {
          return false;
        }
      }
      break;

    case 3:
      return u(0) == v(0) && u(1) == v(1) && u(2) == v(2);
      break;

    case 2:
      return u(0) == v(0) && u(1) == v(1);
      break;

  }

  return true;
}

//
// R^N, vector inequality tested by components
// \param u
// \param v the operands
// \return \f$ u \neq v \f$
//
template<typename T>
inline
bool
operator!=(Vector<T> const & u, Vector<T> const & v)
{
  return !(u == v);
}

//
// R^N scalar vector product
// \param s scalar factor
// \param u vector factor
// \return \f$ s u \f$
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_vector< Promote<S,T> > >::type
operator*(S const & s, Vector<T> const & u)
{
  Index const
  N = u.get_dimension();

  Vector<typename Promote<S, T>::type>
  v(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        v(i) = s * u(i);
      }
      break;

    case 3:
      v(0) = s * u(0);
      v(1) = s * u(1);
      v(2) = s * u(2);
      break;

    case 2:
      v(0) = s * u(0);
      v(1) = s * u(1);
      break;

  }

  return v;
}

//
// R^N vector scalar product
// \param u vector factor
// \param s scalar factor
// \return \f$ s u \f$
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_vector< Promote<S,T> > >::type
operator*(Vector<T> const & u, S const & s)
{
  return s * u;
}

//
// R^N vector scalar division
// \param u vector
// \param s scalar that divides each component of vector
// \return \f$ u / s \f$
//
template<typename S, typename T>
inline
Vector<typename Promote<S, T>::type>
operator/(Vector<T> const & u, S const & s)
{
  Index const
  N = u.get_dimension();

  Vector<typename Promote<S, T>::type>
  v(N);

  switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        v(i) = u(i) / s;
      }
      break;

    case 3:
      v(0) = u(0) / s;
      v(1) = u(1) / s;
      v(2) = u(2) / s;
      break;

    case 2:
      v(0) = u(0) / s;
      v(1) = u(1) / s;
      break;

  }

  return v;
}

//
// R^N vector dot product
// \param u
// \param v operands
// \return \f$ u \cdot v \f$
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
// R^N with N != 3 will produce an error.
// \param u
// \param v operands
// \return \f$ u \times v \f$
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
