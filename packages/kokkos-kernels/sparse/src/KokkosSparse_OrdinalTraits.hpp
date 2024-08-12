//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_SPARSE_ORDINALTRAITS_HPP_
#define KOKKOS_SPARSE_ORDINALTRAITS_HPP_

/// \file KokkosSparse_OrdinalTraits.hpp
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
template <class T>
struct OrdinalTraits {
  static KOKKOS_INLINE_FUNCTION T invalid() { return -1; }
};

// template<>
// struct OrdinalTraits<char> {
//   static KOKKOS_INLINE_FUNCTION char invalid () { return CHAR_MAX; }
// };

template <>
struct OrdinalTraits<short int> {
  static constexpr KOKKOS_INLINE_FUNCTION short int invalid() { return -1; }
};

template <>
struct OrdinalTraits<unsigned short int> {
  static constexpr KOKKOS_INLINE_FUNCTION unsigned short int invalid() { return USHRT_MAX; }
};

template <>
struct OrdinalTraits<int> {
  static constexpr KOKKOS_INLINE_FUNCTION int invalid() { return -1; }
};

template <>
struct OrdinalTraits<unsigned int> {
  static constexpr KOKKOS_INLINE_FUNCTION unsigned int invalid() { return UINT_MAX; }
};

template <>
struct OrdinalTraits<long> {
  static constexpr KOKKOS_INLINE_FUNCTION long invalid() { return -1; }
};

template <>
struct OrdinalTraits<unsigned long> {
  static constexpr KOKKOS_INLINE_FUNCTION unsigned long invalid() { return ULONG_MAX; }
};

template <>
struct OrdinalTraits<long long> {
  static constexpr KOKKOS_INLINE_FUNCTION long long invalid() { return -1; }
};

template <>
struct OrdinalTraits<unsigned long long> {
  static constexpr KOKKOS_INLINE_FUNCTION unsigned long long invalid() { return ULLONG_MAX; }
};

}  // namespace KokkosSparse

#endif  // KOKKOS_SPARSE_ORDINALTRAITS_HPP_
