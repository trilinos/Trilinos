// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Tensor4_h)
#define MiniTensor_Tensor4_h

#include "MiniTensor_Tensor3.h"

namespace minitensor {

template<typename T, Index N>
using tensor4_store = Storage<T, dimension_power<N, 4>::value>;

///
/// Fourth-order tensor.
///
template<typename T, Index N = DYNAMIC>
class Tensor4 : public TensorBase<T, tensor4_store<T, N>>
{
public:

  ///
  /// Order
  ///
  static constexpr
  Index
  ORDER = 4;

  ///
  /// Static or dynamic
  ///
  static constexpr
  bool
  IS_DYNAMIC = N == DYNAMIC;

  ///
  /// Storage type
  ///
  using Store = tensor4_store<T, N>;

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
  /// 4th-order tensor constructor with NaNs
  /// \param dimension the space dimension
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4();

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4(Index const dimension);

  ///
  /// Create 4th-order tensor from a specified value
  /// \param dimension the space dimension
  /// \param value all components are set equal to this
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4(Filler const value);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4(Index const dimension, Filler const value);

  ///
  /// Create 4th-order tensor from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4(T const * data_ptr);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4(Index const dimension, T const * data_ptr);

  ///
  /// Copy constructor
  /// 4th-order tensor constructor with 4th-order tensor
  ///
  KOKKOS_INLINE_FUNCTION
  Tensor4(Tensor4<T, N> const & A);

  ///
  /// 4th-order tensor from 2nd-order tensor
  ///
  KOKKOS_INLINE_FUNCTION
  Tensor4(Tensor<T, dimension_square<N>::value> const & A);

  ///
  /// 4th-order tensor simple destructor
  ///
  virtual
  KOKKOS_INLINE_FUNCTION
  ~Tensor4();

  ///
  /// Indexing for constant 4th-order tensor
  /// \param i index
  /// \param j index
  /// \param k index
  /// \param l index
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  operator()(
      Index const i,
      Index const j,
      Index const k,
      Index const l) const;

  ///
  /// 4th-order tensor indexing
  /// \param i index
  /// \param j index
  /// \param k index
  /// \param l index
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  operator()(
      Index const i,
      Index const j,
      Index const k,
      Index const l);

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
/// 4th-order tensor addition
/// \param A 4th-order tensor
/// \param B 4th-order tensor
/// \return \f$ A + B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator+(Tensor4<S, N> const & A, Tensor4<T, N> const & B);

///
/// 4th-order tensor substraction
/// \param A 4th-order tensor
/// \param B 4th-order tensor
/// \return \f$ A - B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator-(Tensor4<S, N> const & A, Tensor4<T, N> const & B);

///
/// 4th-order tensor minus
/// \return \f$ -A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>
operator-(Tensor4<T, N> const & A);

///
/// 4th-order equality
/// Tested by components
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor4<T, N> const & A, Tensor4<T, N> const & B);

///
/// 4th-order inequality
/// Tested by components
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor4<T, N> const & A, Tensor4<T, N> const & B);

///
/// Scalar 4th-order tensor product
/// \param s scalar
/// \param A 4th-order tensor
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T>, N>>::type
operator*(S const & s, Tensor4<T, N> const & A);

///
/// 4th-order tensor scalar product
/// \param A 4th-order tensor
/// \param s scalar
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T>, N>>::type
operator*(Tensor4<T, N> const & A, S const & s);

///
/// 4th-order tensor scalar division
/// \param A 4th-order tensor
/// \param s scalar
/// \return \f$ A / s \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator/(Tensor4<T, N> const & A, S const & s);

///
/// 4th-order scalar tensor division
/// \param s scalar
/// \param A 4th-order tensor
/// \return \f$ s / A \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator/(S const & s, Tensor4<T, N> const & A);

///
/// 4th-order tensor transpose
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>
transpose(Tensor4<T, N> const & A);

///
/// 4th-order identity I1
/// \return \f$ \delta_{ik} \delta_{jl} \f$ such that \f$ A = I_1 A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_1();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor4<T, DYNAMIC> const
identity_1(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_1(Index const dimension);

///
/// 4th-order identity I2
/// \return \f$ \delta_{il} \delta_{jk} \f$ such that \f$ A^T = I_2 A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_2();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor4<T, DYNAMIC> const
identity_2(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_2(Index const dimension);

///
/// 4th-order identity I3
/// \return \f$ \delta_{ij} \delta_{kl} \f$ such that \f$ I_A I = I_3 A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_3();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor4<T, DYNAMIC> const
identity_3(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_3(Index const dimension);

///
/// Levi-Civita symbol
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
levi_civita_4();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor4<T, DYNAMIC> const
levi_civita_4(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
levi_civita_4(Index const dimension);

///
/// Permutation symbol
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
permutation_4();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor4<T, DYNAMIC> const
permutation_4(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
permutation_4(Index const dimension);

///
/// Alternating symbol
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
alternator_4();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor4<T, DYNAMIC> const
alternator_4(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
alternator_4(Index const dimension);

///
/// 4th-order inverse
/// \return \f$ B such that B : A = A : B = I_1 \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>
inverse(Tensor4<T, N> const & A);

///
/// 4th-order tensor vector dot product
/// \param A 4th-order tensor
/// \param u vector
/// \return 3rd-order tensor \f$ B = A \cdot u := B_{ijk}=A_{ijkp} u_{p} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Tensor4<T, N> const & A, Vector<S, N> const & u);

///
/// vector 4th-order tensor dot product
/// \param A 4th-order tensor
/// \param u vector
/// \return 3rd-order tensor \f$ u dot A \f$ as \f$ B_{ijk}=u_{p} A_{pijk} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Vector<S, N> const & u, Tensor4<T, N> const & A);

///
/// 4th-order tensor vector dot2 product
/// \param A 4th-order tensor
/// \param u vector
/// \return 3rd-order tensor \f$ B = A \cdot u := B_{ijk} = A_{ijpk} u_{p} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Tensor4<T, N> const & A, Vector<S, N> const & u);

///
/// vector 4th-order tensor dot2 product
/// \param A 4th-order tensor
/// \param u vector
/// \return 3rd-order tensor \f$ u dot2 A \f$ as \f$ B_{ijk}=u_{p} A_{ipjk} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Vector<S, N> const & u, Tensor4<T, N> const & A);

///
/// 4th-order tensor 2nd-order tensor double dot product
/// \param A 4th-order tensor
/// \param B 2nd-order tensor
/// \return 2nd-order tensor \f$ C = A : B := C_{ij} = A_{ijpq} B_{pq} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dotdot(Tensor4<T, N> const & A, Tensor<S, N> const & B);

///
/// 2nd-order tensor 4th-order tensor double dot product
/// \param B 2nd-order tensor
/// \param A 4th-order tensor
/// \return 2nd-order tensor \f$ C = B : A := C_{ij} = B_{pq} A_{pqij} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dotdot(Tensor<S, N> const & B, Tensor4<T, N> const & A);

///
/// 4th-order tensor 4th-order tensor double dot product
/// \param A 4th-order tensor
/// \param B 4th-order tensor
/// \return 2nd-order tensor \f$ C = A : B := C_{ij} = A_{ijpq} B_{pq} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dotdot(Tensor4<S, N> const & A, Tensor4<T, N> const & B);

///
/// 2nd-order tensor 2nd-order tensor tensor product
/// \param A 2nd-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C = A \otimes B := C_{ijkl} = A_{ij} B_{kl} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
tensor(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// 2nd-order tensor 2nd-order tensor tensor product
/// \param A 2nd-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C_{ijkl} = A_{ik} B_{jl} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
tensor2(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// 2nd-order tensor 2nd-order tensor tensor product
/// \param A 2nd-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C_{ijkl} = A_{il} B_{kj} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
tensor3(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// 4th-order tensor 2nd-order tensor dot product
/// \param A 4th-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C = A \cdot B := C_{ijkl} = A_{ijkp} B_{pl} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot(Tensor4<T, N> const & A, Tensor<S, N> const & B);

///
/// 4th-order tensor 2nd-order tensor transpose dot product
/// \param A 4th-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C = A \cdot B^T := C_{ijkl} = A_{ijkp} B_{lp} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot_t(Tensor4<T, N> const & A, Tensor<S, N> const & B);

///
/// 2nd-order tensor 4th-order tensor dot product
/// \param A 2nd-order tensor
/// \param B 4th-order tensor
/// \return \f$ C = A \cdot B := C_{ijkl} = A_{ip} B_{pjkl} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot(Tensor<S> const & A, Tensor4<T, N> const & B);

///
/// 2nd-order tensor transpose 4th-order tensor dot product
/// \param A 2nd-order tensor
/// \param B 4th-order tensor
/// \return \f$ C = A^T \cdot B := C_{ijkl} = A_{pi} B_{pjkl} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
t_dot(Tensor<S, N> const & A, Tensor4<T, N> const & B);

///
/// 4th-order tensor 2nd-order tensor dot product
/// \param A 4th-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C = A \cdot B := C_{ijkl} = A_{ijpl} B_{pk} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot2(Tensor4<T, N> const & A, Tensor<S, N> const & B);

///
/// 4th-order tensor 2nd-order tensor transpose dot product
/// \param A 4th-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C = A \cdot B^T := C_{ijkl} = A_{ijpl} B_{kp} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot2_t(Tensor4<T, N> const & A, Tensor<S, N> const & B);

///
/// 2nd-order tensor 4th-order tensor dot product
/// \param A 2nd-order tensor
/// \param B 4th-order tensor
/// \return \f$ C = A \cdot B := C_{ijkl} = A_{jp} B_{ipkl} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot2(Tensor<S, N> const & A, Tensor4<T, N> const & B);

///
/// 2nd-order tensor transpose 4th-order tensor dot product
/// \param A 2nd-order tensor
/// \param B 4th-order tensor
/// \return \f$ C = A^T \cdot B := C_{ijkl} = A_{pj} B_{ipkl} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
t_dot2(Tensor<S, N> const & A, Tensor4<T, N> const & B);

///
/// odot operator useful for \f$ \frac{\partial A^{-1}}{\partial A} \f$
/// see Holzapfel eqn 6.165
/// \param A 2nd-order tensor
/// \param B 2nd-order tensor
/// \return \f$ A \odot B \f$ which is
/// \f$ C_{ijkl} = \frac{1}{2}(A_{ik} B_{jl} + A_{il} B_{jk}) \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
odot(Tensor<S, N> const & A, Tensor<T, N> const & B);

///
/// 4th-order input
/// \param A 4th-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C'_{i'j'k'l'} = A_{i'i} A_{j'j} A_{k'k} A_{l'l} B_{ijkl} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
kronecker(Tensor<S, N> const & A, Tensor4<T, N> const & B);

///
/// 4th-order input
/// \param A 4th-order tensor
/// \param is input stream
/// \return is input stream
///
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Tensor4<T, N> & A);

///
/// 4th-order output
/// \param A 4th-order tensor
/// \param os output stream
/// \return os output stream
///
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Tensor4<T, N> const & A);

} // namespace minitensor

#include "MiniTensor_Tensor4.i.h"
#include "MiniTensor_Tensor4.t.h"

#endif //MiniTensor_Tensor4_h
