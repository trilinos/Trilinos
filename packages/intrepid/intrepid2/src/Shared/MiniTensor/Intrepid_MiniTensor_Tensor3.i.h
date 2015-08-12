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

#if !defined(Intrepid_MiniTensor_Tensor3_i_h)
#define Intrepid_MiniTensor_Tensor3_i_h

namespace Intrepid {

//
// 3rd-order tensor constructor with NaNs
//
template<typename T, Index N>
inline
Tensor3<T, N>::Tensor3() :
TensorBase<T, Store>::TensorBase()
{
  return;
}

template<typename T, Index N>
inline
Tensor3<T, N>::Tensor3(Index const dimension) :
TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

//
// 3rd-order tensor constructor with a specified value
//
template<typename T, Index N>
inline
Tensor3<T, N>::Tensor3(ComponentValue const value) :
TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

template<typename T, Index N>
inline
Tensor3<T, N>::Tensor3(Index const dimension, ComponentValue const value) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

//
//  Create 3rd-order tensor from array
//
template<typename T, Index N>
inline
Tensor3<T, N>::Tensor3(T const * data_ptr) :
TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N>
inline
Tensor3<T, N>::Tensor3(Index const dimension, T const * data_ptr) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index N>
inline
Tensor3<T, N>::Tensor3(Tensor3<T, N> const & A) :
TensorBase<T, Store>::TensorBase(A)
{
  return;
}

//
// 3rd-order tensor simple destructor
//
template<typename T, Index N>
inline
Tensor3<T, N>::~Tensor3()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N>
inline
Index
Tensor3<T, N>::get_dimension() const
{
  return IS_DYNAMIC == true ? TensorBase<T, Store>::get_dimension() : N;
}

//
// Set dimension
//
template<typename T, Index N>
inline
void
Tensor3<T, N>::set_dimension(Index const dimension)
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
// 3rd-order tensor addition
//
template<typename S, typename T, Index N>
inline
Tensor3<typename Promote<S, T>::type, N>
operator+(Tensor3<S, N> const & A, Tensor3<T, N> const & B)
{
  Tensor3<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  add(A, B, C);

  return C;
}

//
// 3rd-order tensor subtraction
//
template<typename S, typename T, Index N>
inline
Tensor3<typename Promote<S, T>::type, N>
operator-(Tensor3<S, N> const & A, Tensor3<T, N> const & B)
{
  Tensor3<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  subtract(A, B, C);

  return C;
}

//
// 3rd-order tensor minus
//
template<typename T, Index N>
inline
Tensor3<T, N>
operator-(Tensor3<T, N> const & A)
{
  Tensor3<T, N>
  B(A.get_dimension());

  minus(A, B);

  return B;
}

//
// 3rd-order tensor equality
//
template<typename T, Index N>
inline
bool
operator==(Tensor3<T, N> const & A, Tensor3<T, N> const & B)
{
  return equal(A, B);
}

//
// 3rd-order tensor inequality
//
template<typename T, Index N>
inline
bool
operator!=(Tensor3<T, N> const & A, Tensor3<T, N> const & B)
{
  return not_equal(A, B);
}

//
// Scalar 3rd-order tensor product
//
template<typename S, typename T, Index N>
inline
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N> >::type
operator*(S const & s, Tensor3<T, N> const & A)
{
  Tensor3<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 3rd-order tensor scalar product
//
template<typename S, typename T, Index N>
inline
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N> >::type
operator*(Tensor3<T, N> const & A, S const & s)
{
  Tensor3<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 3rd-order tensor scalar division
//
template<typename S, typename T, Index N>
inline
Tensor3<typename Promote<S, T>::type, N>
operator/(Tensor3<T, N> const & A, S const & s)
{
  Tensor3<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  divide(A, s, B);

  return B;
}

//
// 3rd-order scalar tensor division
//
template<typename S, typename T, Index N>
inline
Tensor3<typename Promote<S, T>::type, N>
operator/(S const & s, Tensor3<T, N> const & A)
{
  Tensor3<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  split(A, s, B);

  return B;
}

//
// Indexing for constant 3rd order tensor
//
template<typename T, Index N>
inline
T const &
Tensor3<T, N>::operator()(Index const i, Index const j, Index const k) const
{
  Tensor3<T, N> const &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[(i * dimension + j) * dimension + k];
}

//
// 3rd-order tensor indexing
//
template<typename T, Index N>
inline
T &
Tensor3<T, N>::operator()(Index const i, Index const j, Index const k)
{
  Tensor3<T, N> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[(i * dimension + j) * dimension + k];
}

// Local utility functions
namespace {

template<typename T, Index N>
inline
void ones_in_diagonal(Tensor3<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  switch (dimension) {

  default:
    for (Index i = 0; i < dimension; ++i) {
      A(i, i, i) = 1.0;
    }
    break;

  case 3:
    A(0, 0, 0) = 1.0;
    A(1, 1, 1) = 1.0;
    A(2, 2, 2) = 1.0;
    break;

  case 2:
    A(0, 0, 0) = 1.0;
    A(1, 1, 1) = 1.0;
    break;

  }

  return;
}

template<typename T, Index N>
inline
void fill_levi_civita(Tensor3<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        A(i, j, k) = levi_civita<T>(i, j, k);
      }
    }
  }

  return;
}

} // anonymous namespace

//
// Levi-Civita symbol
//
template<typename T, Index N>
inline
Tensor3<T, N> const
levi_civita_3()
{
  Tensor3<T, N>
  A(N, ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T>
inline
Tensor3<T, DYNAMIC> const
levi_civita_3(Index const dimension)
{
  Tensor3<T, DYNAMIC>
  A(dimension, ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T, Index N>
inline
Tensor3<T, N> const
levi_civita_3(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor3<T, DYNAMIC>
  A(dimension, ZEROS);

  fill_levi_civita(A);

  return A;
}

//
// Permutation symbol
//
template<typename T, Index N>
inline
Tensor3<T, N> const
permutation_3()
{
  return levi_civita_3<T, N>();
}

template<typename T>
inline
Tensor3<T, DYNAMIC> const
permutation_3(Index const dimension)
{
  return levi_civita_3<T>(dimension);
}

template<typename T, Index N>
inline
Tensor3<T, N> const
permutation_3(Index const dimension)
{
  return levi_civita_3<T, N>(dimension);
}

//
// Alternating symbol
//
template<typename T, Index N>
inline
Tensor3<T, N> const
alternator_3()
{
  return levi_civita_3<T, N>();
}

template<typename T>
inline
Tensor3<T, DYNAMIC> const
alternator_3(Index const dimension)
{
  return levi_civita_3<T>(dimension);
}

template<typename T, Index N>
inline
Tensor3<T, N> const
alternator_3(Index const dimension)
{
  return levi_civita_3<T, N>(dimension);
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Tensor3_i_h
