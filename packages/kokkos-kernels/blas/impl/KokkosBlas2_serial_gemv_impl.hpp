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
#ifndef __KOKKOSBLAS_GEMV_SERIAL_IMPL_HPP__
#define __KOKKOSBLAS_GEMV_SERIAL_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBlas_util.hpp"
#include "KokkosBlas2_serial_gemv_internal.hpp"

namespace KokkosBlas {

template <typename ArgTrans, typename ArgAlgo>
struct SerialGemv {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType /*alpha*/, const AViewType & /*A*/, const xViewType & /*x*/,
                                           const ScalarType /*beta*/, const yViewType & /*y*/);
};

}  // namespace KokkosBlas

#include "KokkosBlas2_serial_gemv_tpl_spec_decl.hpp"

namespace KokkosBlas {

///
/// Serial Impl
/// ===========

///
/// NT
///

template <>
template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
KOKKOS_INLINE_FUNCTION int SerialGemv<Trans::NoTranspose, Algo::Gemv::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const xViewType &x, const ScalarType beta, const yViewType &y) {
  return Impl::SerialGemvInternal<Algo::Gemv::Unblocked>::invoke(A.extent(0), A.extent(1), alpha, A.data(),
                                                                 A.stride_0(), A.stride_1(), x.data(), x.stride_0(),
                                                                 beta, y.data(), y.stride_0());
}

template <>
template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
KOKKOS_INLINE_FUNCTION int SerialGemv<Trans::NoTranspose, Algo::Gemv::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const xViewType &x, const ScalarType beta, const yViewType &y) {
  return Impl::SerialGemvInternal<Algo::Gemv::Blocked>::invoke(A.extent(0), A.extent(1), alpha, A.data(), A.stride_0(),
                                                               A.stride_1(), x.data(), x.stride_0(), beta, y.data(),
                                                               y.stride_0());
}

///
/// T
///

template <>
template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
KOKKOS_INLINE_FUNCTION int SerialGemv<Trans::Transpose, Algo::Gemv::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const xViewType &x, const ScalarType beta, const yViewType &y) {
  return Impl::SerialGemvInternal<Algo::Gemv::Unblocked>::invoke(A.extent(1), A.extent(0), alpha, A.data(),
                                                                 A.stride_1(), A.stride_0(), x.data(), x.stride_0(),
                                                                 beta, y.data(), y.stride_0());
}

template <>
template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
KOKKOS_INLINE_FUNCTION int SerialGemv<Trans::Transpose, Algo::Gemv::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const xViewType &x, const ScalarType beta, const yViewType &y) {
  return Impl::SerialGemvInternal<Algo::Gemv::Blocked>::invoke(A.extent(1), A.extent(0), alpha, A.data(), A.stride_1(),
                                                               A.stride_0(), x.data(), x.stride_0(), beta, y.data(),
                                                               y.stride_0());
}

///
/// CT
///

template <>
template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
KOKKOS_INLINE_FUNCTION int SerialGemv<Trans::ConjTranspose, Algo::Gemv::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const xViewType &x, const ScalarType beta, const yViewType &y) {
  return Impl::SerialGemvInternal<Algo::Gemv::Unblocked>::invoke(Impl::OpConj(), A.extent(1), A.extent(0), alpha,
                                                                 A.data(), A.stride_1(), A.stride_0(), x.data(),
                                                                 x.stride_0(), beta, y.data(), y.stride_0());
}

template <>
template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
KOKKOS_INLINE_FUNCTION int SerialGemv<Trans::ConjTranspose, Algo::Gemv::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const xViewType &x, const ScalarType beta, const yViewType &y) {
  return Impl::SerialGemvInternal<Algo::Gemv::Blocked>::invoke(Impl::OpConj(), A.extent(1), A.extent(0), alpha,
                                                               A.data(), A.stride_1(), A.stride_0(), x.data(),
                                                               x.stride_0(), beta, y.data(), y.stride_0());
}

}  // namespace KokkosBlas

#endif
