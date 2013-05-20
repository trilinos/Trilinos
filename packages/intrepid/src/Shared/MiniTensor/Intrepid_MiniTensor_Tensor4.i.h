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

#if !defined(Intrepid_MiniTensor_Tensor4_i_h)
#define Intrepid_MiniTensor_Tensor4_i_h

namespace Intrepid
{

//
// 4th-order tensor default constructor
//
template<typename T>
inline
Tensor4<T>::Tensor4() :
TensorBase<T>::TensorBase()
{
  return;
}

//
// 4th-order tensor constructor with NaNs
//
template<typename T>
inline
Tensor4<T>::Tensor4(Index const dimension) :
TensorBase<T>::TensorBase(dimension, order)
{
  return;
}

//
// 4th-order tensor constructor with a specified value
//
template<typename T>
inline
Tensor4<T>::Tensor4(Index const dimension, ComponentValue value) :
TensorBase<T>::TensorBase(dimension, order, value)
{
  return;
}

//
// 4th-order tensor constructor with a scalar
//
template<typename T>
inline
Tensor4<T>::Tensor4(Index const dimension, T const & s) :
TensorBase<T>::TensorBase(dimension, order, s)
{
  return;
}

//
//  Create 4th-order tensor from array
//
template<typename T>
inline
Tensor4<T>::Tensor4(Index const dimension, T const * data_ptr) :
TensorBase<T>::TensorBase(dimension, order, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T>
inline
Tensor4<T>::Tensor4(Tensor4<T> const & A) :
TensorBase<T>::TensorBase(A)
{
  return;
}

//
// 4th-order tensor simple destructor
//
template<typename T>
inline
Tensor4<T>::~Tensor4()
{
  return;
}

//
// 4th-order tensor addition
//
template<typename S, typename T>
inline
Tensor4<typename Promote<S, T>::type>
operator+(Tensor4<S> const & A, Tensor4<T> const & B)
{
  Tensor4<typename Promote<S, T>::type>
  C;

  add(A, B, C);

  return C;
}

//
// 4th-order tensor subtraction
//
template<typename S, typename T>
inline
Tensor4<typename Promote<S, T>::type>
operator-(Tensor4<S> const & A, Tensor4<T> const & B)
{
  Tensor4<typename Promote<S, T>::type>
  C;

  subtract(A, B, C);

  return C;
}

//
// 4th-order tensor minus
//
template<typename T>
inline
Tensor4<T>
operator-(Tensor4<T> const & A)
{
  Tensor4<T>
  B;

  minus(A, B);

  return B;
}

//
// 4th-order equality
//
template<typename T>
inline
bool
operator==(Tensor4<T> const & A, Tensor4<T> const & B)
{
  return equal(A, B);
}

//
// 4th-order inequality
//
template<typename T>
inline
bool
operator!=(Tensor4<T> const & A, Tensor4<T> const & B)
{
  return not_equal(A, B);
}

//
// Scalar 4th-order tensor product
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T> > >::type
operator*(S const & s, Tensor4<T> const & A)
{
  Tensor4<typename Promote<S, T>::type>
  B;

  scale(A, s, B);

  return B;
}

//
// 4th-order tensor scalar product
//
template<typename S, typename T>
inline
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T> > >::type
operator*(Tensor4<T> const & A, S const & s)
{
  Tensor4<typename Promote<S, T>::type>
  B;

  scale(A, s, B);

  return B;
}

//
// 4th-order tensor scalar division
//
template<typename S, typename T>
inline
Tensor4<typename Promote<S, T>::type>
operator/(Tensor4<T> const & A, S const & s)
{
  Tensor4<typename Promote<S, T>::type>
  B;

  divide(A, s, B);

  return B;
}

//
// Indexing for constant 4th order tensor
// \param i index
// \param j index
// \param k index
// \param l index
//
template<typename T>
inline
T const &
Tensor4<T>::operator()(
    Index const i, Index const j, Index const k, Index const l) const
{
  Tensor4<T> const &
  self = (*this);

  Index const
  N = self.get_dimension();

  assert(i < N);
  assert(j < N);
  assert(k < N);
  assert(l < N);

  return self[((i * N + j) * N + k) * N + l];
}

//
// 4th-order tensor indexing
// \param i index
// \param j index
// \param k index
// \param l index
//
template<typename T>
inline
T &
Tensor4<T>::operator()(
    Index const i, Index const j, Index const k, Index const l)
{
  Tensor4<T> &
  self = (*this);

  Index const
  N = self.get_dimension();

  assert(i < N);
  assert(j < N);
  assert(k < N);
  assert(l < N);

  return self[((i * N + j) * N + k) * N + l];
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Tensor4_i_h
