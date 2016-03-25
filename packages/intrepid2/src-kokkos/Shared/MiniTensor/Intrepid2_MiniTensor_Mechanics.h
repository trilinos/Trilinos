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

#if !defined(Intrepid2_MiniTensor_Mechanics_h)
#define Intrepid2_MiniTensor_Mechanics_h

#include "Intrepid2_MiniTensor_Tensor.h"

namespace Intrepid2 {

///
/// Volumetric part of 2nd-order tensor
/// \param A tensor
/// \return \f$ \frac{1}{3} \mathrm{tr}\:A I \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
vol(Tensor<T, N, ES> const & A);

///
/// Deviatoric part of 2nd-order tensor
/// \param A tensor
/// \return \f$ A - vol(A) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
dev(Tensor<T, N, ES> const & A);

///
/// Push forward covariant vector
/// \param \f$ F, u \f$
/// \return \f$ F^{-T} u \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
push_forward_covariant(Tensor<T, N, ES> const & F, Vector<T, N, ES> const & u);

///
/// Pull back covariant vector
/// \param \f$ F, u \f$
/// \return \f$ F^T u \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
pull_back_covariant(Tensor<T, N, ES> const & F, Vector<T, N, ES> const & u);

///
/// Push forward contravariant vector
/// \param \f$ F, u \f$
/// \return \f$ F u \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
push_forward_contravariant(Tensor<T, N, ES> const & F, Vector<T, N, ES> const & u);

///
/// Pull back contravariant vector
/// \param \f$ F, u \f$
/// \return \f$ F^{-1} u \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
pull_back_contravariant(Tensor<T, N, ES> const & F, Vector<T, N, ES> const & u);

///
/// Push forward covariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F^{-T} A F^{-1} \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
push_forward_covariant(Tensor<T, N, ES> const & F, Tensor<T, N, ES> const & A);

///
/// Pull back covariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F^T A F\f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
pull_back_covariant(Tensor<T, N, ES> const & F, Tensor<T, N, ES> const & A);

///
/// Push forward contravariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F A F^T \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
push_forward_contravariant(Tensor<T, N, ES> const & F, Tensor<T, N, ES> const & A);

///
/// Pull back contravariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F^{-1} A F^{-T} \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
pull_back_contravariant(Tensor<T, N, ES> const & F, Tensor<T, N, ES> const & A);

///
/// Piola transformation for vector
/// \param \f$ F, u \f$
/// \return \f$ \det F F^{-1} u \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
piola(Tensor<T, N, ES> const & F, Vector<T, N, ES> const & u);

///
/// Inverse Piola transformation for vector
/// \param \f$ F, u \f$
/// \return \f$ (\det F)^{-1} F u \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
piola_inverse(Tensor<T, N, ES> const & F, Vector<T, N, ES> const & u);

///
/// Piola transformation for tensor, applied on second
/// index. Useful for transforming Cauchy stress to 1PK stress.
/// \param \f$ F, sigma \f$
/// \return \f$ \det F sigma F^{-T} \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
piola(Tensor<T, N, ES> const & F, Tensor<T, N, ES> const & sigma);

///
/// Inverse Piola transformation for tensor, applied on second
/// index. Useful for transforming 1PK stress to Cauchy stress.
/// \param \f$ F, P \f$
/// \return \f$ (\det F)^{-1} P F^T \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
piola_inverse(Tensor<T, N, ES> const & F, Tensor<T, N, ES> const & P);

///
/// Smallest eigenvalue by inverse iteration.
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
smallest_eigenvalue(Tensor<T, N, ES> const & A);

///
/// Check strict ellipticity condition for 4th-order tensor.
/// Assume A has major symmetries.
/// \param A 4th-order tensor is transformed into 2nd-order
/// tensor.
/// \return whether the smallest eigenvalue of the 2nd-order
/// tensor is less or equal than zero.
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
check_strict_ellipticity(Tensor4<T, N, ES> const & A);

///
/// Check strong ellipticity condition for 4th-order tensor.
/// Assume A has major and minor symmetries.
/// \param A 4th-order tensor.
/// \return whether \f$ (m\odot n):A:(m\odot n) > 0 \forall m,n \neq 0 \f$.
///
template<typename T, Index N,  typename ES>
std::pair<bool, Vector<T, N, ES> >
check_strong_ellipticity(Tensor4<T, N, ES> const & A);

} // namespace Intrepid

#include "Intrepid2_MiniTensor_Mechanics.i.h"
#include "Intrepid2_MiniTensor_Mechanics.t.h"

#endif // Intrepid2_MiniTensor_Mechanics_h
