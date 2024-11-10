// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Mechanics_h)
#define MiniTensor_Mechanics_h

#include "MiniTensor_Tensor.h"

namespace minitensor {

///
/// Volumetric part of 2nd-order tensor
/// \param A tensor
/// \return \f$ \frac{1}{3} \mathrm{tr}\:A I \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
vol(Tensor<T, N> const & A);

///
/// Deviatoric part of 2nd-order tensor
/// \param A tensor
/// \return \f$ A - vol(A) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
dev(Tensor<T, N> const & A);

///
/// Push forward covariant vector
/// \param \f$ F, u \f$
/// \return \f$ F^{-T} u \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
push_forward_covariant(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Pull back covariant vector
/// \param \f$ F, u \f$
/// \return \f$ F^T u \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
pull_back_covariant(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Push forward contravariant vector
/// \param \f$ F, u \f$
/// \return \f$ F u \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
push_forward_contravariant(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Pull back contravariant vector
/// \param \f$ F, u \f$
/// \return \f$ F^{-1} u \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
pull_back_contravariant(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Push forward covariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F^{-T} A F^{-1} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
push_forward_covariant(Tensor<T, N> const & F, Tensor<T, N> const & A);

///
/// Pull back covariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F^T A F\f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
pull_back_covariant(Tensor<T, N> const & F, Tensor<T, N> const & A);

///
/// Push forward contravariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F A F^T \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
push_forward_contravariant(Tensor<T, N> const & F, Tensor<T, N> const & A);

///
/// Pull back contravariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F^{-1} A F^{-T} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
pull_back_contravariant(Tensor<T, N> const & F, Tensor<T, N> const & A);

///
/// Piola transformation for vector
/// \param \f$ F, u \f$
/// \return \f$ \det F F^{-1} u \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
piola(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Inverse Piola transformation for vector
/// \param \f$ F, u \f$
/// \return \f$ (\det F)^{-1} F u \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
piola_inverse(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Piola transformation for tensor, applied on second
/// index. Useful for transforming Cauchy stress to 1PK stress.
/// \param \f$ F, sigma \f$
/// \return \f$ \det F sigma F^{-T} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
piola(Tensor<T, N> const & F, Tensor<T, N> const & sigma);

///
/// Inverse Piola transformation for tensor, applied on second
/// index. Useful for transforming 1PK stress to Cauchy stress.
/// \param \f$ F, P \f$
/// \return \f$ (\det F)^{-1} P F^T \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
piola_inverse(Tensor<T, N> const & F, Tensor<T, N> const & P);

///
/// Smallest eigenvalue by inverse iteration.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
smallest_eigenvalue(Tensor<T, N> const & A);

///
/// Check strict ellipticity condition for 4th-order tensor.
/// Assume A has major symmetries.
/// \param A 4th-order tensor is transformed into 2nd-order
/// tensor.
/// \return whether the smallest eigenvalue of the 2nd-order
/// tensor is less or equal than zero.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
check_strict_ellipticity(Tensor4<T, N> const & A);

///
/// Check strong ellipticity condition for 4th-order tensor.
/// Assume A has major and minor symmetries.
/// \param A 4th-order tensor.
/// \return whether \f$ (m\odot n):A:(m\odot n) > 0 \forall m,n \neq 0 \f$.
///
template<typename T, Index N>
std::pair<bool, Vector<T, N>>
check_strong_ellipticity(Tensor4<T, N> const & A);

} // namespace minitensor

#include "MiniTensor_Mechanics.i.h"
#include "MiniTensor_Mechanics.t.h"

#endif // MiniTensor_Mechanics_h
