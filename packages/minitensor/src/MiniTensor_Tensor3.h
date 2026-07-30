// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Tensor3_h)
#define MiniTensor_Tensor3_h

#include <iomanip>

#include "MiniTensor_Tensor.h"

namespace minitensor {

template<typename T, Index N>
using tensor3_store = Storage<T, dimension_power<N, 3>::value>;

///
/// Third-order tensor.
///
template<typename T, Index N = DYNAMIC>
class Tensor3 : public TensorBase<T, tensor3_store<T, N>>
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
  using Store = tensor3_store<T, N>;

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
  Tensor3(Filler const value);

  explicit
  KOKKOS_INLINE_FUNCTION
  Tensor3(Index const dimension, Filler const value);

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
  Tensor3(Tensor3<T, N> const & A);

  ///
  /// 3rd-order tensor simple destructor
  ///
  virtual
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
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
operator+(Tensor3<S, N> const & A, Tensor3<T, N> const & B);

///
/// 3rd-order tensor substraction
/// \param A 3rd-order tensor
/// \param B 3rd-order tensor
/// \return \f$ A - B \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
operator-(Tensor3<S, N> const & A, Tensor3<T, N> const & B);

///
/// 3rd-order tensor minus
/// \return \f$ -A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>
operator-(Tensor3<T, N> const & A);

///
/// 3rd-order tensor equality
/// Tested by components
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor3<T, N> const & A, Tensor3<T, N> const & B);

///
/// 3rd-order tensor inequality
/// Tested by components
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor3<T, N> const & A, Tensor3<T, N> const & B);

///
/// Scalar 3rd-order tensor product
/// \param s scalar
/// \param A 3rd-order tensor
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N>>::type
operator*(S const & s, Tensor3<T, N> const & A);

///
/// 3rd-order tensor scalar product
/// \param A 3rd-order tensor
/// \param s scalar
/// \return \f$ s A \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N>>::type
operator*(Tensor3<T, N> const & A, S const & s);

///
/// 3rd-order tensor scalar division
/// \param A 3rd-order tensor
/// \param s scalar
/// \return \f$ A / s \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
operator/(Tensor3<T, N> const & A, S const & s);

///
/// 3rd-order scalar tensor division
/// \param s scalar
/// \param A 3rd-order tensor
/// \return \f$ s / A \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
operator/(S const & s, Tensor3<T, N> const & A);

///
/// 3rd-order tensor 2nd-order tensor double dot product
/// \param A 3rd-order tensor
/// \param u 2nd-order tensor
/// \return \f$ B = A : u := B_i = A_{ijk} u_{jk} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
dotdot(Tensor3<T, N> const & A, Tensor<S, N> const & u);

///
/// 3rd-order tensor vector product
/// \param A 3rd-order tensor
/// \param u vector
/// \return \f$ B = A \cdot u := B_{ij} = A_{ijp} u_p \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot(Tensor3<T, N> const & A, Vector<S, N> const & u);

///
/// vector 3rd-order tensor product
/// \param A 3rd-order tensor
/// \param u vector
/// \return \f$ B = u \cdot A := B_{ij} = u_p A{pij} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot(Vector<S, N> const & u, Tensor3<T, N> const & A);

///
/// 3rd-order tensor vector product
/// \param A 3rd-order tensor
/// \param u vector
/// \return \f$ B = A \cdot u := B_{ij} = A_{ipj} u_p \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot2(Tensor3<T, N> const & A, Vector<S> const & u);

///
/// vector 3rd-order tensor product
/// \param u vector
/// \param A 3rd-order tensor
/// \return \f$ B = u \cdot A := B_{ij} = u_p A_{ipj} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot2(Vector<S, N> const & u, Tensor3<T, N> const & A);

///
/// 3rd-order tensor 2nd-order tensor product
/// \param A 3rd-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ijp} B_{pk} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Tensor3<T, N> const & A, Tensor<S, N> const & B);

///
/// 2nd-order tensor 3rd-order tensor product
/// \param A 2nd-order tensor
/// \param B 3rd-order tensor
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{pjk} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Tensor<S, N> const & A, Tensor3<T, N> const & B);

///
/// 3rd-order tensor 2nd-order tensor product
/// \param A 3rd-order tensor
/// \param B 2nd-order tensor
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ipj} B_{pk} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Tensor3<T, N> const & A, Tensor<S, N> const & B);

///
/// 2nd-order tensor 3rd-order tensor product
/// \param A 2nd-order tensor
/// \param B 3rd-order tensor
/// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{jpk} \f$
///
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Tensor<S, N> const & A, Tensor3<T, N> const & B);

///
/// Levi-Civita symbol
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
levi_civita_3();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC> const
levi_civita_3(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
levi_civita_3(Index const dimension);

///
/// Permutation symbol
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
permutation_3();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC> const
permutation_3(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
permutation_3(Index const dimension);

///
/// Alternating symbol
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
alternator_3();

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC> const
alternator_3(Index const dimension);

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
alternator_3(Index const dimension);

///
/// 3rd-order tensor input
/// \param A 3rd-order tensor
/// \param is input stream
/// \return is input stream
///
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Tensor3<T, N> & A);

///
/// 3rd-order tensor output
/// \param A 3rd-order tensor
/// \param os output stream
/// \return os output stream
///
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Tensor3<T, N> const & A);

} // namespace minitensor

namespace minitensor {

//
// 3rd-order tensor constructor with NaNs
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3() :
TensorBase<T, Store>::TensorBase()
{
  set_dimension(N);
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Index const dimension) :
TensorBase<T, Store>::TensorBase(dimension, ORDER)
{
  return;
}

//
// 3rd-order tensor constructor with a specified value
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Filler const value) :
TensorBase<T, Store>::TensorBase(N, ORDER, value)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Index const dimension, Filler const value) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  return;
}

//
//  Create 3rd-order tensor from array
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(T const * data_ptr) :
TensorBase<T, Store>::TensorBase(N, ORDER, data_ptr)
{
  return;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Index const dimension, T const * data_ptr) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  return;
}

//
// Copy constructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::Tensor3(Tensor3<T, N> const & A) :
TensorBase<T, Store>::TensorBase(A)
{
  return;
}

//
// 3rd-order tensor simple destructor
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N>::~Tensor3()
{
  return;
}

//
// Get dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Index
Tensor3<T, N>::get_dimension() const
{
  return TensorBase<T, Store>::get_dimension(ORDER);
}

//
// Set dimension
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
Tensor3<T, N>::set_dimension(Index const dimension)
{
  TensorBase<T, Store>::set_dimension(dimension, ORDER);
  return;
}

//
// 3rd-order tensor addition
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
operator+(Tensor3<S, N>const & A, Tensor3<T, N> const & B)
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
bool
operator==(Tensor3<T, N> const & A, Tensor3<T, N> const & B)
{
  return equal(A, B);
}

//
// 3rd-order tensor inequality
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Tensor3<T, N> const & A, Tensor3<T, N> const & B)
{
  return not_equal(A, B);
}

//
// Scalar 3rd-order tensor product
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N>>::type
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
KOKKOS_INLINE_FUNCTION
typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T>, N>>::type
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
levi_civita_3()
{
  Tensor3<T, N>
  A(N, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC> const
levi_civita_3(Index const dimension)
{
  Tensor3<T, DYNAMIC>
  A(dimension, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
levi_civita_3(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor3<T, DYNAMIC>
  A(dimension, Filler::ZEROS);

  fill_levi_civita(A);

  return A;
}
// Permutation symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
permutation_3()
{
  return levi_civita_3<T, N>();
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC> const
permutation_3(Index const dimension)
{
  return levi_civita_3<T>(dimension);
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
permutation_3(Index const dimension)
{
  return levi_civita_3<T, N>(dimension);
}

//
// Alternating symbol
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
alternator_3()
{
  return levi_civita_3<T, N>();
}

template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor3<T, DYNAMIC> const
alternator_3(Index const dimension)
{
  return levi_civita_3<T>(dimension);
}

template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<T, N> const
alternator_3(Index const dimension)
{
  return levi_civita_3<T, N>(dimension);
}

} // namespace minitensor
namespace minitensor {

//
// \return \f$ B = A : u := B_i = A_{ijk} u_{jk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
dotdot(Tensor3<T, N> const & A, Tensor<S, N> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Vector<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {

    typename Promote<S, T>::type
    s = 0.0;

    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        s += A(i,j,k) * u(j,k);
      }
    }
    B(i) = s;
  }

  return B;
}

//
// \return \f$ B = A \cdot u := B_{ij} = A_{ijp} u_p \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot(Tensor3<T, N> const & A, Vector<S, N> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  B(N);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(i,j,p) * u(p);
      }
      B(i,j) = s;
    }
  }

  return B;
}

//
// \return \f$ B = u \cdot A := B_{ij} = u_p A{pij} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot(Vector<S, N> const & u, Tensor3<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += u(p) * A(p,i,j);
      }
      B(i,j) = s;
    }
  }

  return B;
}


//
// \return \f$ B = A \cdot u := B_{ij} = A_{ipj} u_p \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot2(Tensor3<T, N> const & A, Vector<S, N> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(i,p,j) * u(p);
      }
      B(i,j) = s;
    }
  }

  return B;
}

//
// \return \f$ B = u \cdot A := B_{ij} = u_p A_{ipj} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot2(Vector<S, N> const & u, Tensor3<T, N> const & A)
{
  return dot2(A, u);
}

//
// \return \f$ C = A \cdot B := C_{ijk} = A_{ijp} B_{pk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Tensor3<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index k = 0; k < dimension; ++k) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,j,p) * B(p,k);
        }
        C(i,j,k) = s;
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{pjk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Tensor<S, N> const & A, Tensor3<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index k = 0; k < dimension; ++k) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,p) * B(p,j,k);
        }
        C(i,j,k) = s;
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijk} = A_{ipj} B_{pk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Tensor3<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index k = 0; k < dimension; ++k) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,p,j) * B(p,k);
        }
        C(i,j,k) = s;
      }
    }
  }

  return C;
}


//
// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{jpk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Tensor<S, N> const & A, Tensor3<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index k = 0; k < dimension; ++k) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,p) * B(j,p,k);
        }
        C(i,j,k) = s;
      }
    }
  }

  return C;
}


//
// 3rd-order tensor input
// \param A 3rd-order tensor
// \param is input stream
// \return is input stream
//
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Tensor3<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
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
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Tensor3<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  if (dimension == 0) {
    return os;
  }

  os << std::scientific << std::setprecision(17);

  for (Index i = 0; i < dimension; ++i) {

    for (Index j = 0; j < dimension; ++j) {

      os << std::setw(24) << A(i,j,0);

      for (Index k = 1; k < dimension; ++k) {
        os << "," << std::setw(24) << A(i,j,k);
      }

      os << std::endl;

    }

    os << std::endl;
    os << std::endl;

  }

  return os;
}

} // namespace minitensor

#endif //MiniTensor_Tensor3_h
