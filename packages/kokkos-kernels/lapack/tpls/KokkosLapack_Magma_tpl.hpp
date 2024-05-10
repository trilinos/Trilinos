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
#ifndef KOKKOSLAPACK_MAGMA_TPL_HPP_
#define KOKKOSLAPACK_MAGMA_TPL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
#include <KokkosLapack_magma.hpp>

namespace KokkosLapack {
namespace Impl {

MagmaSingleton::MagmaSingleton() {
  magma_int_t stat = magma_init();
  if (stat != MAGMA_SUCCESS) Kokkos::abort("MAGMA initialization failed\n");

  Kokkos::push_finalize_hook([&]() { magma_finalize(); });
}

MagmaSingleton& MagmaSingleton::singleton() {
  static MagmaSingleton s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)

#endif  // KOKKOSLAPACK_MAGMA_TPL_HPP_
