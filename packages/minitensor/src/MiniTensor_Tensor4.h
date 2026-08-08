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

#include <iomanip>

#include "MiniTensor_Tensor3.h"

namespace minitensor {

/// \addtogroup minitensor_containers
/// @{

///
/// Storage type alias for Tensor4.
///
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
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4();

  ///
  /// 4th-order tensor constructor with NaNs
  /// \param dimension the space dimension
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4(Index const dimension);

  ///
  /// Create 4th-order tensor from a specified value
  /// \param value all components are set equal to this
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4(Filler const value);

  ///
  /// Create 4th-order tensor from a specified value
  /// \param dimension the space dimension
  /// \param value all components are set equal to this
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4(Index const dimension, Filler const value);

  ///
  /// Create 4th-order tensor from array
  /// \param data_ptr pointer into the array
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor4(T const * data_ptr);

  ///
  /// Create 4th-order tensor from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
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

///
/// 4th-order identity I1: \f$ \delta_{ik} \delta_{jl} \f$ such that \f$ A = I_1 A \f$.
///
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

///
/// 4th-order identity I2: \f$ \delta_{il} \delta_{jk} \f$ such that \f$ A^T = I_2 A \f$.
///
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

///
/// 4th-order identity I3: \f$ \delta_{ij} \delta_{kl} \f$ such that \f$ I_A I = I_3 A \f$.
///
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

namespace minitensor
{

//
// 4th-order tensor constructor with NaNs
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4() :
TensorBase<T, Store>::TensorBase()
{
  set_dimension(N);
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Index const dimension) :
TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

//
// 4th-order tensor constructor with a specified value
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Filler const value) :
TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Index const dimension, Filler const value) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

//
//  Create 4th-order tensor from array
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(T const * data_ptr) :
TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Index const dimension, T const * data_ptr) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Tensor4<T, N> const & A) :
TensorBase<T, Store>::TensorBase(A)
{
  return;
}

//
// 4th-order tensor from 2nd-order tensor
//

namespace {

KOKKOS_INLINE_FUNCTION
Index
second_to_fourth_dimension(Index const dimension_2nd)
{
  Index
  dimension_4th = 0;

  switch (dimension_2nd) {

  default:
    MT_ERROR_EXIT("Invalid dimension for 2nd-order tensor.");
    break;

  case 1:
    dimension_4th = 1;
    break;

  case 4:
    dimension_4th = 2;
    break;

  case 9:
    dimension_4th = 3;
    break;

  case 16:
    dimension_4th = 4;
    break;

  }

  return dimension_4th;
}

} //anonymous namespace

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::Tensor4(Tensor<T, dimension_square<N>::value> const & A)
{
  Index const
  dimension_2nd = A.get_dimension();

  Index const
  dimension_4th = second_to_fourth_dimension(dimension_2nd);

  Tensor4<T, N> &
  self = (*this);

  self.set_dimension(dimension_4th);

  Index const
  number_components = dimension_2nd * dimension_2nd;

  for (Index i = 0; i < number_components; ++i) {
    self[i] = A[i];
  }

  return;
}
//
// 4th-order tensor simple destructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>::~Tensor4()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Tensor4<T, N>::get_dimension() const
{
  return TensorBase<T, Store>::get_dimension(ORDER);
}

//
// Set dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
Tensor4<T, N>::set_dimension(Index const dimension)
{
  TensorBase<T, Store>::set_dimension(dimension, ORDER);
  return;
}

//
// 4th-order tensor addition
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator+(Tensor4<S, N> const & A, Tensor4<T, N> const & B)
{
  Tensor4<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  add(A, B, C);

  return C;
}

//
// 4th-order tensor subtraction
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator-(Tensor4<S, N> const & A, Tensor4<T, N> const & B)
{
  Tensor4<typename Promote<S, T>::type, N>
  C(A.get_dimension());

  subtract(A, B, C);

  return C;
}

//
// 4th-order tensor minus
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>
operator-(Tensor4<T, N> const & A)
{
  Tensor4<T, N>
  B(A.get_dimension());

  minus(A, B);

  return B;
}

//
// 4th-order equality
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor4<T, N> const & A, Tensor4<T, N> const & B)
{
  return equal(A, B);
}

//
// 4th-order inequality
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor4<T, N> const & A, Tensor4<T, N> const & B)
{
  return not_equal(A, B);
}

//
// Scalar 4th-order tensor product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T>, N>>::type
operator*(S const & s, Tensor4<T, N> const & A)
{
  Tensor4<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 4th-order tensor scalar product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T>, N>>::type
operator*(Tensor4<T, N> const & A, S const & s)
{
  Tensor4<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  scale(A, s, B);

  return B;
}

//
// 4th-order tensor scalar division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator/(Tensor4<T, N> const & A, S const & s)
{
  Tensor4<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  divide(A, s, B);

  return B;
}

//
// 4th-order scalar tensor division
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
operator/(S const & s, Tensor4<T, N> const & A)
{
  Tensor4<typename Promote<S, T>::type, N>
  B(A.get_dimension());

  split(A, s, B);

  return B;
}

//
// Indexing for constant 4th order tensor
// \param i index
// \param j index
// \param k index
// \param l index
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T const &
Tensor4<T, N>::operator()(
    Index const i, Index const j, Index const k, Index const l) const
{
  Tensor4<T, N> const &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[((i * dimension + j) * dimension + k) * dimension + l];
}

//
// 4th-order tensor indexing
// \param i index
// \param j index
// \param k index
// \param l index
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T &
Tensor4<T, N>::operator()(
    Index const i, Index const j, Index const k, Index const l)
{
  Tensor4<T, N> &
  self = (*this);

  Index const
  dimension = self.get_dimension();

  return self[((i * dimension + j) * dimension + k) * dimension + l];
}

//
// 4th-order inverse
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>
inverse(Tensor4<T, N> const & A)
{
  return Tensor4<T, N>(inverse(Tensor<T, dimension_square<N>::value>(A)));
}

} // namespace minitensor
namespace minitensor {

namespace {

template< typename T, Index N>
KOKKOS_INLINE_FUNCTION
void ones_in_ikjl(Tensor4<T, N> & A)
{

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          if (i == k && j == l) {
            A(i,j,k,l) = 1;
          }
        }
      }
    }
  }

  return;
}

template< typename T, Index N>
KOKKOS_INLINE_FUNCTION
void ones_in_iljk(Tensor4<T, N> & A)
{

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          if (i == l && j == k) {
            A(i,j,k,l) = 1;
          }
        }
      }
    }
  }

  return;
}

template< typename T, Index N>
KOKKOS_INLINE_FUNCTION
void ones_in_ijkl(Tensor4<T, N> & A)
{

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          if (i == j && k == l) {
            A(i,j,k,l) = 1;
          }
        }
      }
    }
  }

  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void fill_levi_civita(Tensor4<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          A(i, j, k, l) = levi_civita<T>(i, j, k, l);
        }
      }
    }
  }

  return;
}

} // anonymous namespace

//
// 4th-order identity I1
// \return \f$ \delta_{ik} \delta_{jl} \f$ such that \f$ A = I_1 A \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_1()
{
  Tensor4<T, N> I(N, Filler::ZEROS);
  ones_in_ikjl(I);
  return I;
}

template<typename T>
KOKKOS_INLINE_FUNCTION
const Tensor4<T, DYNAMIC>
identity_1(Index const dimension)
{
  Tensor4<T, DYNAMIC> I(dimension, Filler::ZEROS);
  ones_in_ikjl(I);
  return I;
}

///
/// 4th-order identity I1: \f$ \delta_{ik} \delta_{jl} \f$ such that \f$ A = I_1 A \f$.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_1(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor4<T, N> I(dimension, Filler::ZEROS);
  ones_in_ikjl(I);
  return I;
}

//
// 4th-order identity I2
// \return \f$ \delta_{il} \delta_{jk} \f$ such that \f$ A^T = I_2 A \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_2()
{
  Tensor4<T, N> I(N, Filler::ZEROS);
  ones_in_iljk(I);
  return I;
}

template<typename T>
KOKKOS_INLINE_FUNCTION
const Tensor4<T, DYNAMIC>
identity_2(Index const dimension)
{
  Tensor4<T, DYNAMIC> I(dimension, Filler::ZEROS);
  ones_in_iljk(I);
  return I;
}

///
/// 4th-order identity I2: \f$ \delta_{il} \delta_{jk} \f$ such that \f$ A^T = I_2 A \f$.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_2(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor4<T, N> I(dimension, Filler::ZEROS);
  ones_in_iljk(I);
  return I;
}

//
// 4th-order identity I3
// \return \f$ \delta_{ij} \delta_{kl} \f$ such that \f$ I_A I = I_3 A \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_3()
{
  Tensor4<T, N> I(N, Filler::ZEROS);
  ones_in_ijkl(I);
  return I;
}

template<typename T>
KOKKOS_INLINE_FUNCTION
const Tensor4<T, DYNAMIC>
identity_3(Index const dimension)
{
  Tensor4<T, DYNAMIC> I(dimension, Filler::ZEROS);
  ones_in_ijkl(I);
  return I;
}

///
/// 4th-order identity I3: \f$ \delta_{ij} \delta_{kl} \f$ such that \f$ I_A I = I_3 A \f$.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
identity_3(Index const dimension)
{
  Tensor4<T, N> I(dimension, Filler::ZEROS);
  ones_in_ijkl(I);
  return I;
}

//
// Levi-Civita symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
levi_civita_4()
{
  Tensor4<T, N>
  A(N, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

///
/// Levi-Civita symbol as a 4th-order tensor.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor4<T, DYNAMIC> const
levi_civita_4(Index const dimension)
{
  Tensor4<T, DYNAMIC>
  A(dimension, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

///
/// Levi-Civita symbol as a 4th-order tensor.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
levi_civita_4(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor4<T, DYNAMIC>
  A(dimension, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

//
// Permutation symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
permutation_4()
{
  return levi_civita_4<T, N>();
}

///
/// Permutation symbol as a 4th-order tensor; alias of levi_civita_4.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor4<T, DYNAMIC> const
permutation_4(Index const dimension)
{
  return levi_civita_4<T>(dimension);
}

///
/// Permutation symbol as a 4th-order tensor; alias of levi_civita_4.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
permutation_4(Index const dimension)
{
  return levi_civita_4<T, N>(dimension);
}

//
// Alternating symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
alternator_4()
{
  return levi_civita_4<T, N>();
}

///
/// Alternator as a 4th-order tensor; alias of levi_civita_4.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor4<T, DYNAMIC> const
alternator_4(Index const dimension)
{
  return levi_civita_4<T>(dimension);
}


///
/// Alternator as a 4th-order tensor; alias of levi_civita_4.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N> const
alternator_4(Index const dimension)
{
  return levi_civita_4<T, N>(dimension);
}

//
// 4th-order tensor transpose
// per Holzapfel 1.157
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<T, N>
transpose(Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor4<T, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          B(i, j, k, l) = A(k, l, i, j);
        }
      }
    }
  }

  return B;
}

//
// 4th-order tensor vector dot product
// \param A 4th-order tensor
// \param u vector
// \return 3rd-order tensor \f$ B = A \cdot u := B_{ijk}=A_{ijkp} u_{p} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Tensor4<T, N> const & A, Vector<S, N> const & u)
{
  Index const
  dimension= A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,j,k,p) * u(p);
        }
        B(i,j,k) = s;
      }
    }
  }

  return B;
}

//
// vector 4th-order tensor dot product
// \param A 4th-order tensor
// \param u vector
// \return 3rd-order tensor \f$ u dot A \f$ as \f$ B_{ijk}=u_{p} A_{pijk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Vector<S, N> const & u, Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += u(p) * A(p,i,j,k);
        }
        B(i,j,k) = s;
      }
    }
  }

  return B;
}

//
// 4th-order tensor vector dot2 product
// \param A 4th-order tensor
// \param u vector
// \return 3rd-order tensor \f$ B = A \cdot u := B_{ijk} = A_{ijpk} u_{p} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Tensor4<T, N> const & A, Vector<S, N> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,j,p,k) * u(p);
        }
        B(i,j,k) = s;
      }
    }
  }

  return B;
}

//
// vector 4th-order tensor dot2 product
// \param A 4th-order tensor
// \param u vector
// \return 3rd-order tensor \f$ u dot2 A \f$ as \f$ B_{ijk}=u_{p} A_{ipjk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Vector<S, N> const & u, Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += u(p) * A(i,p,j,k);
        }
        B(i,j,k) = s;
      }
    }
  }

  return B;
}

//
// \return 2nd-order tensor \f$ C = A : B := C_{ij} = A_{ijpq} B_{pq} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dotdot(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        for (Index q = 0; q < dimension; ++q) {
          s += A(i,j,p,q) * B(p,q);
        }
      }
      C(i,j) = s;
    }
  }

  return C;
}

//
// \return 2nd-order tensor \f$ C = B : A := C_{ij} = B_{pq} A_{pqij} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dotdot(Tensor<S, N> const & B, Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        for (Index q = 0; q < dimension; ++q) {
          s += B(p,q) * A(p,q,i,j);
        }
      }
      C(i,j) = s;
    }
  }

  return C;
}

//
// \return \f$ C = A : B := C_{ijkl} = A_{ijpq} B{pqkl} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dotdot(Tensor4<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            for (Index q = 0; q < dimension; ++q) {
              s += A(i,j,p,q) * B(p,q,k,l);
            }
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \otimes B := C_{ijkl} = A_{ij} B_{kl} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
tensor(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          C(i,j,k,l) = A(i,j) * B(k,l);
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C_{ijkl} = A_{ik} B_{jl} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
tensor2(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          C(i,j,k,l) = A(i,k) * B(j,l);
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C_{ijkl} = A_{il} B_{kj} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
tensor3(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          C(i,j,k,l) = A(i,l) * B(k,j);
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijkl} = A_{ijkp} B_{pl} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,j,k,p) * B(p,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B^T := C_{ijkl} = A_{ijkp} B_{lp} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot_t(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,j,k,p) * B(l,p);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijkl} = A_{ip} B_{pjkl} \f$
//
///
/// 2nd-order by 4th-order dot product.
/// \return \f$ C_{ijkl} = A_{ip} B_{pjkl} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,p) * B(p,j,k,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A^T \cdot B := C_{ijkl} = A_{pi} B_{pjkl} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
t_dot(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(p,i) * B(p,j,k,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijkl} = A_{ijpl} B_{pk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot2(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,j,p,l) * B(p,k);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B^T := C_{ijkl} = A_{ijpl} B_{kp} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot2_t(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,j,p,l) * B(k,p);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijkl} = A_{jp} B_{ipkl} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
dot2(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(j,p) * B(i,p,k,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A^T \cdot B := C_{ijkl} = A_{pj} B_{ipkl} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
t_dot2(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(p,j) * B(i,p,k,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// odot operator useful for \f$ \frac{\partial A^{-1}}{\partial A} \f$
// see Holzapfel eqn 6.165
// \param A 2nd-order tensor
// \param B 2nd-order tensor
// \return \f$ A \odot B \f$ which is
// \f$ C_{ijkl} = \frac{1}{2}(A_{ik} B_{jl} + A_{il} B_{jk}) \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
odot(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          C(i,j,k,l) = 0.5 * (A(i,k) * B(j,l) + A(i,l) * B(j,k));
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C'_{i'j'k'l'} = Q_{i'i} Q_{j'j} Q_{k'k} Q_{l'l} C_{ijkl}  \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor4<typename Promote<S, T>::type, N>
kronecker(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          // we assume A is the direction cosine matrix
          // have s_{i}\be_{i} = s_{i'}\be_{i'}, want s_{i}
          // s_{i} = A_{i',i}s_{i'}
          for (Index p = 0; p < dimension; ++p) {
            for (Index q = 0; q < dimension; ++q) {
              for (Index m = 0; m < dimension; ++m) {
                for (Index n = 0; n < dimension; ++n) {
                  s += A(p,i) * A(q,j) * A(m,k) * A(n,l) * B(p,q,m,n);
                }
              }
            }
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// 4th-order input
// \param A 4th-order tensor
// \param is input stream
// \return is input stream
//
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Tensor4<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          is >> A(i,j,k,l);
        }
      }
    }
  }

  return is;
}

//
// 4th-order output
// \param A 4th-order tensor
// \param os output stream
// \return os output stream
//
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  if (dimension == 0) {
    return os;
  }

  os << std::scientific << std::setprecision(17);

  for (Index i = 0; i < dimension; ++i) {

    for (Index j = 0; j < dimension; ++j) {

      for (Index k = 0; k < dimension; ++k) {

        os << std::setw(24) << A(i,j,k,0);

        for (Index l = 1; l < dimension; ++l) {

          os << "," << std::setw(24) << A(i,j,k,l);
        }

        os << std::endl;

      }

      os << std::endl;
      os << std::endl;

    }

    os << std::endl;

  }

  return os;
}

/// @}
} // namespace minitensor

#endif //MiniTensor_Tensor4_h
