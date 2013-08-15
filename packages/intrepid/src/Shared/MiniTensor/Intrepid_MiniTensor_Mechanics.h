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

#if !defined(Intrepid_MiniTensor_Mechanics_h)
#define Intrepid_MiniTensor_Mechanics_h

#include "Intrepid_MiniTensor_Tensor.h"

namespace Intrepid {

///
/// Volumetric part of 2nd-order tensor
/// \param A tensor
/// \return \f$ \frac{1}{3} \mathrm{tr}\:A I \f$
///
template<typename T, Index N>
Tensor<T, N>
vol(Tensor<T, N> const & A);

///
/// Deviatoric part of 2nd-order tensor
/// \param A tensor
/// \return \f$ A - vol(A) \f$
///
template<typename T, Index N>
Tensor<T, N>
dev(Tensor<T, N> const & A);

///
/// Push forward covariant vector
/// \param \f$ F, u \f$
/// \return \f$ F^{-T} u \f$
///
template<typename T, Index N>
Vector<T, N>
push_forward_covariant(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Pull back covariant vector
/// \param \f$ F, u \f$
/// \return \f$ F^T u \f$
///
template<typename T, Index N>
Vector<T, N>
pull_back_covariant(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Push forward contravariant vector
/// \param \f$ F, u \f$
/// \return \f$ F u \f$
///
template<typename T, Index N>
Vector<T, N>
push_forward_contravariant(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Pull back contravariant vector
/// \param \f$ F, u \f$
/// \return \f$ F^{-1} u \f$
///
template<typename T, Index N>
Vector<T, N>
pull_back_contravariant(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Push forward covariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F^{-T} A F^{-1} \f$
///
template<typename T, Index N>
Tensor<T, N>
push_forward_covariant(Tensor<T, N> const & F, Tensor<T, N> const & A);

///
/// Pull back covariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F^T A F\f$
///
template<typename T, Index N>
Tensor<T, N>
pull_back_covariant(Tensor<T, N> const & F, Tensor<T, N> const & A);

///
/// Push forward contravariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F A F^T \f$
///
template<typename T, Index N>
Tensor<T, N>
push_forward_contravariant(Tensor<T, N> const & F, Tensor<T, N> const & A);

///
/// Pull back contravariant tensor
/// \param \f$ F, A \f$
/// \return \f$ F^{-1} A F^{-T} \f$
///
template<typename T, Index N>
Tensor<T, N>
pull_back_contravariant(Tensor<T, N> const & F, Tensor<T, N> const & A);

///
/// Piola transformation for vector
/// \param \f$ F, u \f$
/// \return \f$ \det F F^{-1} u \f$
///
template<typename T, Index N>
Vector<T, N>
piola(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Inverse Piola transformation for vector
/// \param \f$ F, u \f$
/// \return \f$ (\det F)^{-1} F u \f$
///
template<typename T, Index N>
Vector<T, N>
piola_inverse(Tensor<T, N> const & F, Vector<T, N> const & u);

///
/// Piola transformation for tensor, applied on second
/// index. Useful for transforming Cauchy stress to 1PK stress.
/// \param \f$ F, sigma \f$
/// \return \f$ \det F sigma F^{-T} \f$
///
template<typename T, Index N>
Tensor<T, N>
piola(Tensor<T, N> const & F, Tensor<T, N> const & sigma);

///
/// Inverse Piola transformation for tensor, applied on second
/// index. Useful for transforming 1PK stress to Cauchy stress.
/// \param \f$ F, P \f$
/// \return \f$ (\det F)^{-1} P F^T \f$
///
template<typename T, Index N>
Tensor<T, N>
piola_inverse(Tensor<T, N> const & F, Tensor<T, N> const & P);

///
/// Check strict ellipticity condition for 4th-order tensor.
/// Assume A has major symmetries.
/// \param A 4th-order tensor is transformed into 2nd-order
/// tensor.
/// \return whether the smallest eigenvalue of the 2nd-order
/// tensor is less or equal than zero.
///
template<typename T, Index N>
bool
check_strict_ellipticity(Tensor4<T, N> const & A);

///
/// Check strong ellipticity condition for 4th-order tensor.
/// Assume A has major symmetries.
/// \param A 4th-order tensor.
/// \return whether \f$ (m\odot n):A:(m\odot n) > 0 \forall m,n \neq 0 \f$.
///
template<typename T, Index N>
std::pair<bool, Vector<T, N> >
check_strong_ellipticity(Tensor4<T, N> const & A);

} // namespace Intrepid

#include "Intrepid_MiniTensor_Mechanics.i.h"
#include "Intrepid_MiniTensor_Mechanics.t.h"

#endif // Intrepid_MiniTensor_Mechanics_h
