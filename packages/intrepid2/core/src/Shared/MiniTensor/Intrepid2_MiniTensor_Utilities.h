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

#if !defined(Intrepid2_MiniTensor_Utilities_h)
#define Intrepid2_MiniTensor_Utilities_h

#include "Intrepid2_ConfigDefs.hpp"
#include "Sacado.hpp"

namespace Intrepid2 {

//
//swap function
//
template<typename T>
KOKKOS_INLINE_FUNCTION
void
swap(T & a, T & b);

//
// max function
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
max(const T & a,const  T & b);

//
// max function
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
min(const T & a,const  T & b);

///
/// Sign function
///
template<typename T>
KOKKOS_INLINE_FUNCTION
int
sgn(T const & s);

///
/// Copysign function
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
copysign(T const & a, T const & b);

///
/// NaN function. Necessary to choose the proper underlying NaN
/// for non-floating-point types.
/// Assumption: non-floating-point types have a typedef that
/// determines the underlying floating-point type.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
not_a_number();

///
/// Machine epsilon function. Necessary to choose the proper underlying
/// machine epsilon for non-floating-point types.
/// Assumption: non-floating-point types have a typedef that
/// determines the underlying floating-point type.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
machine_epsilon();

///
/// The circle constant
///
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
tau();

///
/// Random number generation. Uniform distribution U(-1,1)
/// which is the Teuchos default (!).
///
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
random();

///
/// Random number generation. Uniform distribution U(0,1).
///
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
random_uniform();

///
/// Random number generation. Normal distribution N(0,1).
///
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
random_normal();

///
/// Compute a non-negative integer power by binary manipulation.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
integer_power(T const & X, Index const exponent);

///
/// Utility for Kronecker delta in 2D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
kronecker_delta(Index const i, Index const j);

///
/// Utility for Kronecker delta in 3D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
kronecker_delta(Index const i, Index const j, Index const k);

///
/// Utility for Kronecker delta in 4D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
kronecker_delta(Index const i, Index const j, Index const k, Index const l);

///
/// Utility for Levi-Civita/permutation/alternating symbol in 2D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
levi_civita(Index const i, Index const j);

///
/// Utility for Levi-Civita/permutation/alternating symbol in 3D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
levi_civita(Index const i, Index const j, Index const k);

///
/// Utility for Levi-Civita/permutation/alternating symbol in 4D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
levi_civita(Index const i, Index const j, Index const k, Index const l);

} // namespace Intrepid

#include "Intrepid2_MiniTensor_Utilities.i.h"

#endif // Intrepid2_MiniTensor_Utilities_h
