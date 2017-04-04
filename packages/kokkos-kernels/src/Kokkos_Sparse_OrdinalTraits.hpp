/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_SPARSE_ORDINALTRAITS_HPP_
#define KOKKOS_SPARSE_ORDINALTRAITS_HPP_

/// \file Kokkos_Sparse_OrdinalTraits.hpp
/// \brief Declaration and definition of KokkosSparse::OrdinalTraits,
///   a traits class for "invalid" (flag) values of integer types that
///   KokkosKernels uses as local ordinals or global ordinals.

#include "KokkosKernels_config.h"
#include "Kokkos_Macros.hpp"
#include <climits>

namespace KokkosSparse {

/// \brief Traits class for "invalid" (flag) values of integer types
///   that Tpetra uses as local ordinals or global ordinals.
///
/// \tparam T Built-in integer type.
///
/// If T is signed, invalid() returns T(-1).  If T is unsigned,
/// invalid() returns the maximum representable value of T.  Do NOT
/// rely on these values!
///
/// I didn't choose the values this class calls "invalid."  For
/// backwards compatibility, they are the same as the values found in
/// Teuchos::OrdinalTraits<T>::invalid().  I can't call
/// Teuchos::OrdinalTraits<T>::invalid() because it is not marked as a
/// Kokkos device function.  I also can't use std::numeric_limits for
/// the same reason.  That's why this traits class needs to exist.
template<class T>
struct OrdinalTraits {
  static KOKKOS_INLINE_FUNCTION T invalid () { return -1; }
};

// template<>
// struct OrdinalTraits<char> {
//   static KOKKOS_INLINE_FUNCTION char invalid () { return CHAR_MAX; }
// };

template<>
struct OrdinalTraits<short int> {
  static KOKKOS_INLINE_FUNCTION short int invalid () { return -1; }
};

template<>
struct OrdinalTraits<unsigned short int> {
  static KOKKOS_INLINE_FUNCTION unsigned short int invalid () { return USHRT_MAX; }
};

template<>
struct OrdinalTraits<int> {
  static KOKKOS_INLINE_FUNCTION int invalid () { return -1; }
};

template<>
struct OrdinalTraits<unsigned int> {
  static KOKKOS_INLINE_FUNCTION unsigned int invalid () { return UINT_MAX; }
};

template<>
struct OrdinalTraits<long> {
  static KOKKOS_INLINE_FUNCTION long invalid () { return -1; }
};

template<>
struct OrdinalTraits<unsigned long> {
  static KOKKOS_INLINE_FUNCTION unsigned long invalid () { return ULONG_MAX; }
};

template<>
struct OrdinalTraits<long long> {
  static KOKKOS_INLINE_FUNCTION long long invalid () { return -1; }
};

template<>
struct OrdinalTraits<unsigned long long> {
  static KOKKOS_INLINE_FUNCTION unsigned long long invalid () { return ULLONG_MAX; }
};

} // namespace KokkosSparse

#endif // KOKKOS_SPARSE_ORDINALTRAITS_HPP_
