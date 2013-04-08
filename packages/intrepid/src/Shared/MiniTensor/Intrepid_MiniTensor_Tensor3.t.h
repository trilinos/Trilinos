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

#if !defined(Intrepid_MiniTensor_Tensor3_t_h)
#define Intrepid_MiniTensor_Tensor3_t_h

namespace Intrepid {

//
// set dimension
//
//
template<typename T>
void
Tensor3<T>::set_dimension(Index const N)
{
  if (N == get_dimension()) return;

  Index const
  number_components = N * N * N;

  e.resize(number_components);

  dimension = N;

  return;
}

//
// 3rd-order tensor default constructor
//
template<typename T>
Tensor3<T>::Tensor3() :
dimension(0),
e(Teuchos::null)
{
  return;
}

//
// 3rd-order tensor constructor with NaNs
//
template<typename T>
Tensor3<T>::Tensor3(Index const N) :
dimension(0),
e(Teuchos::null)
{
  set_dimension(N);

  Index const
  number_components = N * N * N;

  for (Index i = 0; i < number_components; ++i) {
    e[i] = not_a_number<T>();
  }

  return;
}

//
// R^N 3rd-order tensor constructor with a scalar
// \param s all components set to this scalar
//
template<typename T>
Tensor3<T>::Tensor3(Index const N, T const & s) :
dimension(0),
e(Teuchos::null)
{
  set_dimension(N);

  Index const
  number_components = N * N * N;

  for (Index i = 0; i < number_components; ++i) {
    e[i] = s;
  }

  return;
}

//
// R^N copy constructor
// 3rd-order tensor constructor from 3rd-order tensor
// \param A from which components are copied
//
template<typename T>
Tensor3<T>::Tensor3(Tensor3<T> const & A) :
dimension(0),
e(Teuchos::null)
{
  Index const
  N = A.get_dimension();

  set_dimension(N);

  Index const
  number_components = N * N * N;

  for (Index i = 0; i < number_components; ++i) {
    e[i] = A.e[i];
  }

  return;
}

//
// R^N 3rd-order tensor simple destructor
//
template<typename T>
Tensor3<T>::~Tensor3()
{
  return;
}

//
// R^N 3rd-order tensor copy assignment
//
template<typename T>
Tensor3<T> &
Tensor3<T>::operator=(Tensor3<T> const & A)
{
  if (this != &A) return *this;

  Index const
  N = A.get_dimension();

  set_dimension(N);

  Index const
  number_components = N * N * N;

  for (Index i = 0; i < number_components; ++i) {
    e[i] = A.e[i];
  }

  return *this;
}

//
// 3rd-order tensor increment
// \param A added to this tensor
//
template<typename T>
Tensor3<T> &
Tensor3<T>::operator+=(Tensor3<T> const & A)
{
  Index const
  N = get_dimension();

  assert(A.get_dimension() == N);

  Index const
  number_components = N * N * N;

  for (Index i = 0; i < number_components; ++i) {
    e[i] += A.e[i];
  }

  return *this;
}

//
// 3rd-order tensor decrement
// \param A substracted from this tensor
//
template<typename T>
Tensor3<T> &
Tensor3<T>::operator-=(Tensor3<T> const & A)
{
  Index const
  N = get_dimension();

  assert(A.get_dimension() == N);

  Index const
  number_components = N * N * N;

  for (Index i = 0; i < number_components; ++i) {
    e[i] -= A.e[i];
  }

  return *this;
}

//
// R^N fill 3rd-order tensor with zeros
//
template<typename T>
void
Tensor3<T>::clear()
{
  Index const
  N = get_dimension();

  Index const
  number_components = N * N * N;

  for (Index i = 0; i < number_components; ++i) {
    e[i] = 0.0;;
  }

  return;
}

//
// 3rd-order tensor addition
// \param A 3rd-order tensor
// \param B 3rd-order tensor
// \return \f$ A + B \f$
//
template<typename S, typename T>
Tensor3<typename Promote<S, T>::type>
operator+(Tensor3<S> const & A, Tensor3<T> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  Tensor3<typename Promote<S, T>::type>
  C(N);

  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {
      for (Index k = 0; k < N; ++k) {
        C(i,j,k) = A(i,j,k) + B(i,j,k);
      }
    }
  }

  return C;
}

//
// 3rd-order tensor substraction
// \param A 3rd-order tensor
// \param B 3rd-order tensor
// \return \f$ A - B \f$
//
template<typename S, typename T>
Tensor3<typename Promote<S, T>::type>
operator-(Tensor3<S> const & A, Tensor3<T> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  Tensor3<typename Promote<S, T>::type>
  C(N);

  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {
      for (Index k = 0; k < N; ++k) {
        C(i,j,k) = A(i,j,k) - B(i,j,k);
      }
    }
  }

  return C;
}

//
// 3rd-order tensor minus
// \return \f$ -A \f$
//
template<typename T>
Tensor3<T>
operator-(Tensor3<T> const & A)
{
  Index const
  N = A.get_dimension();

  Tensor3<T>
  S(N);

  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {
      for (Index k = 0; k < N; ++k) {
        S(i,j,k) = - A(i,j,k);
      }
    }
  }

  return S;
}

//
// 3rd-order tensor equality
// Tested by components
//
template<typename T>
inline bool
operator==(Tensor3<T> const & A, Tensor3<T> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {
      for (Index k = 0; k < N; ++k) {
        if (A(i,j,k) != B(i,j,k)) {
          return false;
        }
      }
    }
  }

  return true;
}

//
// 3rd-order tensor inequality
// Tested by components
//
template<typename T>
inline bool
operator!=(Tensor3<T> const & A, Tensor3<T> const & B)
{
  return !(A==B);
}

//
// Scalar 3rd-order tensor product
// \param s scalar
// \param A 3rd-order tensor
// \return \f$ s A \f$
//
template<typename S, typename T>
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T> > >::type
operator*(S const & s, Tensor3<T> const & A)
{
  Index const
  N = A.get_dimension();

  Tensor3<typename Promote<S, T>::type>
  B(N);

  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {
      for (Index k = 0; k < N; ++k) {
        B(i,j,k) = s * A(i,j,k);
      }
    }
  }

  return B;
}

//
// 3rd-order tensor scalar product
// \param A 3rd-order tensor
// \param s scalar
// \return \f$ s A \f$
//
template<typename S, typename T>
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T> > >::type
operator*(Tensor3<T> const & A, S const & s)
{
  return s * A;
}

//
// 3rd-order tensor scalar division
// \param A 3rd-order tensor
// \param s scalar
// \return \f$ s A \f$
//
template<typename S, typename T>
Tensor3<typename Promote<S, T>::type>
operator/(Tensor3<T> const & A, S const & s)
{
  Index const
  N = A.get_dimension();

  Tensor3<typename Promote<S, T>::type>
  B(N);

  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {
      for (Index k = 0; k < N; ++k) {
        B(i,j,k) = A(i,j,k) / s;
      }
    }
  }

  return B;
}

//
// 3rd-order tensor vector product
// \param A 3rd-order tensor
// \param u vector
// \return \f$ A u \f$
//
template<typename S, typename T>
Tensor<typename Promote<S, T>::type>
dot(Tensor3<T> const & A, Vector<S> const & u)
{
  Index const
  N = A.get_dimension();

  assert(u.get_dimension() == N);

  Tensor<typename Promote<S, T>::type>
  B(N);

  for (Index j = 0; j < N; ++j) {
    for (Index k = 0; k < N; ++k) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index i = 0; i < N; ++i) {
        s += A(i,j,k) * u(i);
      }
      B(j,k) = s;
    }
  }

  return B;
}

//
// vector 3rd-order tensor product
// \param A 3rd-order tensor
// \param u vector
// \return \f$ u A \f$
//
template<typename S, typename T>
Tensor<typename Promote<S, T>::type>
dot(Vector<S> const & u, Tensor3<T> const & A)
{
  Index const
  N = A.get_dimension();

  assert(u.get_dimension() == N);

  Tensor<typename Promote<S, T>::type>
  B(N);

  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index k = 0; k < N; ++k) {
        s += A(i,j,k) * u(k);
      }
      B(i,j) = s;
    }
  }

  return B;
}


//
// 3rd-order tensor vector product2 (contract 2nd index)
// \param A 3rd-order tensor
// \param u vector
// \return \f$ A u \f$
//
template<typename S, typename T>
Tensor<typename Promote<S, T>::type>
dot2(Tensor3<T> const & A, Vector<S> const & u)
{
  Index const
  N = A.get_dimension();

  assert(u.get_dimension() == N);

  Tensor<typename Promote<S, T>::type>
  B(N);

  for (Index i = 0; i < N; ++i) {
    for (Index k = 0; k < N; ++k) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index j = 0; j < N; ++j) {
        s += A(i,j,k) * u(j);
      }
      B(i,k) = s;
    }
  }

  return B;
}

//
// vector 3rd-order tensor product2 (contract 2nd index)
// \param A 3rd-order tensor
// \param u vector
// \return \f$ u A \f$
//
template<typename S, typename T>
Tensor<typename Promote<S, T>::type>
dot2(Vector<S> const & u, Tensor3<T> const & A)
{
  return dot2(A, u);
}

///
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ijp} B_{pk} \f$
///
template<typename S, typename T>
Tensor3<typename Promote<S, T>::type>
dot(Tensor3<T> const & A, Tensor<S> const & B)
{
  Index const
  N = A.get_dimension();

  assert(B.get_dimension() == N);

  Tensor3<typename Promote<S, T>::type>
  C(N);

  for (Index i = 0; i < N; ++i) {
    for (Index k = 0; k < N; ++k) {
      for (Index j = 0; j < N; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < N; ++p) {
          s += A(i,j,p) * B(p,k);
        }
        C(i,j,k) = s;
      }
    }
  }

  return C;
}

///
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{pjk} \f$
///
template<typename S, typename T>
Tensor3<typename Promote<S, T>::type>
dot(Tensor<S> const & A, Tensor3<T> const & B);

///
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ipj} B_{pk} \f$
///
template<typename S, typename T>
Tensor3<typename Promote<S, T>::type>
dot2(Tensor3<T> const & A, Tensor<S> const & B);

///
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{jpk} \f$
///
template<typename S, typename T>
Tensor3<typename Promote<S, T>::type>
dot2(Tensor<S> const & A, Tensor3<T> const & B);

//
// 3rd-order tensor input
// \param A 3rd-order tensor
// \param is input stream
// \return is input stream
//
template<typename T>
std::istream &
operator>>(std::istream & is, Tensor3<T> & A)
{
  Index const
  N = A.get_dimension();

  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {
      for (Index k = 0; k < N; ++k) {
        is >> A(i,j,k);
      }
    }
  }

  return is;
}

//
// 3rd-order tensor output
// \param A 3rd-order tensor
// \param os output stream
// \return os output stream
//
template<typename T>
std::ostream &
operator<<(std::ostream & os, Tensor3<T> const & A)
{
  Index const
  N = A.get_dimension();

  if (N == 0) {
    return os;
  }

  for (Index i = 0; i < N; ++i) {

    for (Index j = 0; j < N; ++j) {

      os << std::scientific << A(i,j,0);

      for (Index k = 1; k < N; ++k) {
        os << std::scientific << "," << A(i,j,k);
      }

      os << std::endl;

    }

    os << std::endl;
    os << std::endl;

  }

  return os;
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Tensor3_t_h
