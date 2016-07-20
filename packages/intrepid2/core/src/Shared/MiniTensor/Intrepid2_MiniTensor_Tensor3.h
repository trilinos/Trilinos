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

#if !defined(Intrepid2_MiniTensor_Tensor3_h)
#define Intrepid2_MiniTensor_Tensor3_h

#include "Intrepid2_MiniTensor_Tensor.h"

namespace Intrepid2 {

template<typename T, Index N,  typename ES>
using tensor3_store = Storage<T, dimension_power<N, 3>::value, ES>;

///
/// Third-order tensor.
///
template<typename T, Index N = DYNAMIC,  typename ES=NOKOKKOS>
class Tensor3 : public TensorBase<T, tensor3_store<T, N, ES>>
{
public:

  ///
  /// Order
  ///
  static constexpr
  Index
  ORDER = 3;

  ///
  /// Static or dynamic
  ///
  static constexpr
  bool
  IS_DYNAMIC = N == DYNAMIC;

  ///
  /// Storage type
  ///
  using Store = tensor3_store<T, N, ES>;

  ///
  /// Tensor order
  ///
  KOKKOS_INLINE_FUNCTION
  static constexpr
  Index
  get_order()
  {
    return ORDER;
  }

  ///
  /// 3rd-order tensor constructor with NaNs
  /// \param dimension the space dimension
  ///
  KOKKOS_INLINE_FUNCTION
  explicit
  Tensor3();

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor3(Index const dimension);

  ///
  /// Create 3rd-order tensor from a specified value
  /// \param dimension the space dimension
  /// \param value all components are set equal to this
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor3(ComponentValue const value);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor3(Index const dimension, ComponentValue const value);

  ///
  /// Create 3rd-order tensor from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor3(T const * data_ptr);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor3(Index const dimension, T const * data_ptr);

  ///
  /// Copy constructor
  /// 3rd-order tensor constructor from 3rd-order tensor
  ///
  KOKKOS_INLINE_FUNCTION
  Tensor3(Tensor3<T, N, ES> const & A);

  ///
  /// 3rd-order tensor simple destructor
  ///
  KOKKOS_INLINE_FUNCTION
  ~Tensor3();

  ///
  /// Indexing for constant 3rd-order tensor
  /// \param i index
  /// \param j index
  /// \param k index
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  operator()(Index const i, Index const j, Index const k) const;

  ///
  /// 3rd-order tensor indexing
  /// \param i index
  /// \param j index
  /// \param k index
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  operator()(Index const i, Index const j, Index const k);

  ///
  /// \return dimension
  ///
  KOKKOS_INLINE_FUNCTION
  Index
  get_dimension() const;

  ///
  /// \param dimension of vector
  ///
  KOKKOS_INLINE_FUNCTION
  void
  set_dimension(Index const dimension);

};

///
/// 3rd-order tensor addition
/// \param A 3rd-order tensor
/// \param B 3rd-order tensor
/// \return \f$ A + B \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
operator+(Tensor3<S, N, ES> const & A, Tensor3<T, N, ES> const & B);

///
/// 3rd-order tensor substraction
/// \param A 3rd-order tensor
/// \param B 3rd-order tensor
/// \return \f$ A - B \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
operator-(Tensor3<S, N, ES> const & A, Tensor3<T, N, ES> const & B);

///
/// 3rd-order tensor minus
/// \return \f$ -A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES>
operator-(Tensor3<T, N, ES> const & A);

///
/// 3rd-order tensor equality
/// Tested by components
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor3<T, N, ES> const & A, Tensor3<T, N, ES> const & B);

///
/// 3rd-order tensor inequality
/// Tested by components
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor3<T, N, ES> const & A, Tensor3<T, N, ES> const & B);

///
/// Scalar 3rd-order tensor product
/// \param s scalar
/// \param A 3rd-order tensor
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N, ES>>::type
operator*(S const & s, Tensor3<T, N, ES> const & A);

///
/// 3rd-order tensor scalar product
/// \param A 3rd-order tensor
/// \param s scalar
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N, ES>>::type
operator*(Tensor3<T, N, ES> const & A, S const & s);

///
/// 3rd-order tensor scalar division
/// \param A 3rd-order tensor
/// \param s scalar
/// \return \f$ A / s \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
operator/(Tensor3<T, N, ES> const & A, S const & s);

///
/// 3rd-order scalar tensor division
/// \param s scalar
/// \param A 3rd-order tensor
/// \return \f$ s / A \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
operator/(S const & s, Tensor3<T, N, ES> const & A);

///
/// 3rd-order tensor 2nd-order tensor double dot product
/// \param A 3rd-order tensor
/// \param u 2nd-order tensor
/// \return \f$ B = A : u := B_i = A_{ijk} u_{jk} \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N, ES>
dotdot(Tensor3<T, N, ES> const & A, Tensor<S, N, ES> const & u);

///
/// 3rd-order tensor vector product
/// \param A 3rd-order tensor
/// \param u vector
/// \return \f$ B = A \cdot u := B_{ij} = A_{ijp} u_p \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dot(Tensor3<T, N, ES> const & A, Vector<S, N, ES> const & u);

///
/// vector 3rd-order tensor product
/// \param A 3rd-order tensor
/// \param u vector
/// \return \f$ B = u \cdot A := B_{ij} = u_p A{pij} \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dot(Vector<S, N, ES> const & u, Tensor3<T, N, ES> const & A);

///
/// 3rd-order tensor vector product
/// \param A 3rd-order tensor
/// \param u vector
/// \return \f$ B = A \cdot u := B_{ij} = A_{ipj} u_p \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dot2(Tensor3<T, N, ES> const & A, Vector<S> const & u);

///
/// vector 3rd-order tensor product
/// \param u vector
/// \param A 3rd-order tensor
/// \return \f$ B = u \cdot A := B_{ij} = u_p A_{ipj} \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N, ES>
dot2(Vector<S, N, ES> const & u, Tensor3<T, N, ES> const & A);

///
/// 3rd-order tensor 2nd-order tensor product
/// \param A 3rd-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ijp} B_{pk} \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
dot(Tensor3<T, N, ES> const & A, Tensor<S, N, ES> const & B);

///
/// 2nd-order tensor 3rd-order tensor product
/// \param A 2nd-order tensor
/// \param B 3rd-order tensor
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{pjk} \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
dot(Tensor<S, N, ES> const & A, Tensor3<T, N, ES> const & B);

///
/// 3rd-order tensor 2nd-order tensor product
/// \param A 3rd-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ipj} B_{pk} \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
dot2(Tensor3<T, N, ES> const & A, Tensor<S, N, ES> const & B);

///
/// 2nd-order tensor 3rd-order tensor product
/// \param A 2nd-order tensor
/// \param B 3rd-order tensor
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{jpk} \f$
///
template<typename S, typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N, ES>
dot2(Tensor<S, N, ES> const & A, Tensor3<T, N, ES> const & B);

///
/// Levi-Civita symbol
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
levi_civita_3();

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC,ES> const
levi_civita_3(Index const dimension);

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
levi_civita_3(Index const dimension);

///
/// Permutation symbol
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
permutation_3();

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC,ES> const
permutation_3(Index const dimension);

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
permutation_3(Index const dimension);

///
/// Alternating symbol
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
alternator_3();

template<typename T,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC,ES> const
alternator_3(Index const dimension);

template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N, ES> const
alternator_3(Index const dimension);

///
/// 3rd-order tensor input
/// \param A 3rd-order tensor
/// \param is input stream
/// \return is input stream
///
template<typename T, Index N,  typename ES>
std::istream &
operator>>(std::istream & is, Tensor3<T, N, ES> & A);

///
/// 3rd-order tensor output
/// \param A 3rd-order tensor
/// \param os output stream
/// \return os output stream
///
template<typename T, Index N,  typename ES>
std::ostream &
operator<<(std::ostream & os, Tensor3<T, N, ES> const & A);

} // namespace Intrepid

#include "Intrepid2_MiniTensor_Tensor3.i.h"
#include "Intrepid2_MiniTensor_Tensor3.t.h"

#endif //Intrepid2_MiniTensor_Tensor3_h
