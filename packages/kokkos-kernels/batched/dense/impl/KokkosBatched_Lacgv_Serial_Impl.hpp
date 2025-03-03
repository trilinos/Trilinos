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

#ifndef KOKKOSBATCHED_LACGV_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_LACGV_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Lacgv_Serial_Internal.hpp"

namespace KokkosBatched {
template <typename XViewType>
KOKKOS_INLINE_FUNCTION int SerialLacgv::invoke(const XViewType &x) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::lacgv: XViewType is not a Kokkos::View.");
  static_assert(XViewType::rank == 1, "KokkosBatched::lacgv: XViewType must have rank 1.");

  return Impl::SerialLacgvInternal::invoke(x.extent_int(0), x.data(), x.stride(0));
}
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_LACGV_SERIAL_IMPL_HPP_
